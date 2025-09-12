
using DataFrames
using Images
using SparseArrays
using LsqFit
using Statistics
using SparseArrays

#using ..JPELT: pelt_fit, steps_from_pelt, getsteps_pelt, pelt_find_optimal_noise_region
#using ..SimpleLogger: debug
#using ..StepAnalysis: fit_traces, getsteps


"""
Exponential fit results: A₀*exp(-t/τ) + C₀

Fields: dataX, fitX, A0, tau, C0, χ2_red
"""
struct ExponentialFit{T<:Real}
    dataX::Vector{T}
    fitX::Vector{T}
    A0::Float64
    tau::Float64
    C0::Float64
    χ2_red::Float64
end

"""
Automatic Step Finder (ASF) fit results.

Fields: dataX, fitX, S_curve, best_shot, iter, steps, stepHeights, stepTimes, χ2_red
"""
struct AsfFit{T<:Real}
    dataX::Vector{T}
    fitX::Vector{T}
    S_curve::Vector{T}
    best_shot::Real
    iter::Int
    steps::Vector{Int}
    stepHeights::Vector{T}
    stepTimes::Vector{Int}
    χ2_red::Float64
end

"""
PELT (Pruned Exact Linear Time) changepoint detection results.

Fields: dataX, fitX, cost, penalty, nstep, stepHeight, stepTime, stepLength, si, sf, χ2_red
"""
struct PeltFit{T<:Real}
	dataX::Vector{T}
    fitX::Vector{T}
	cost::Float64
	penalty::Float64
    nstep::Int
    stepHeight::Vector{Float64}
    stepTime::Vector{Int}
    stepLength::Vector{Int}
    si::Int
    sf::Int
    χ2_red::Float64
end


"""
    chi2_dof(trace, pfit, n_params)

Calculate reduced chi-squared from residuals.
"""
function chi2_dof(trace::Vector{T}, pfit::Vector{T}, n_params::Int) where {T<:Real}
    residuals = (trace .- pfit) #/ sqrt.(trace)
    χ2 = sum(residuals .^ 2)
    dof = length(trace) - n_params
    
    if dof <= 0 
        return nothing
    end
    
    return χ2 / dof 
end

"""
    fit_traces_exponential(trace; initial_params=nothing, max_iterations=1000)

Fit single trace to exponential decay A₀*exp(-t/τ) + C₀.
Returns ExponentialFit or nothing if fit fails.
"""
function fit_traces_exponential(trace::Vector{T}; 
                                initial_params=nothing,
                                max_iterations::Int=1000) where {T<:Real}
    
   
    # Check for empty trace
    if isempty(trace)
        return nothing
    end
    
    # Define exponential decay model
    exponential_model(t, p) = p[1] .* exp.(-t ./ p[2]) .+ p[3]
    
    t = collect(1:length(trace))
    y_float = Float64.(trace)
    
    # Determine initial parameters
    if initial_params !== nothing
        p0 = collect(Float64.(initial_params))
    else
        # Automatic parameter estimation
        y_max = maximum(y_float)
        y_min = minimum(y_float)
        p0 = [y_max - y_min, 10.0, y_min]
    end
    
    # Ensure positive initial parameters where needed
    p0[1] = max(p0[1], 1e-6)  # Amplitude must be positive
    p0[2] = max(p0[2], 1e-6)  # Tau must be positive
    
    # Attempt curve fitting
    fit_result = try
        curve_fit(exponential_model, Float64.(t), y_float, p0; 
                 maxIter=max_iterations)
    catch err
         @warn "Empty trace provided or fit failed"
        return nothing
    end

    # Check convergence
    if !fit_result.converged
        @warn "Fit did not converge"
        return nothing
    end
    
    # Calculate fitted values and statistics
    fitted_vals = exponential_model(Float64.(t), fit_result.param)
    χ2_red = chi2_dof(y_float, fitted_vals, length(fit_result.param))
    
    # Check dof
    if χ2_red == nothing
        @warn "no degrees of freedom"
        return nothing
    end
    
    # Return struct with results
    return ExponentialFit(
        trace,
        fitted_vals,
        fit_result.param[1],      # A0
        fit_result.param[2],      # tau
        fit_result.param[3],      # C0
        χ2_red
    )
end

"""
    fit_traces_exponential(TRZS, df; initial_params=nothing, max_iterations=1000)

Fit multiple traces from sparse matrix to exponential decay.
Returns DataFrame with fit results and trace data arrays.
"""
function fit_traces_exponential(TRZS::SparseMatrixCSC{Vector{T}, Int}, df::DataFrame;
                                initial_params=nothing,max_iterations::Int=1000) where {T<:Real}


    # Pre-allocate result arrays
    trace_indices_i = Int[]
    trace_indices_j = Int[]  
    original_data = Vector{Vector{T}}()
    fitted_data = Vector{Vector{Float64}}()
	
    n_traces = nrow(df)
    fit_results = DataFrame(
        i = Int[], j = Int[],
        A0 = Float64[], tau = Float64[], C0 = Float64[],
        chi2 = Float64[]
    )
    
    for row in eachrow(df)
        i, j = row.i, row.j
        y = TRZS[i, j]
        
        # Use helper function for fitting
        expfit = fit_traces_exponential(y; initial_params=initial_params, max_iterations=max_iterations)
        
        # Handle failed fits
        if expfit === nothing
            @warn "Empty trace or fit failed at ($i, $j), skipping"
            continue
        end
        
        # Store successful fit
        push!(fit_results, (
            i, j,  expfit.A0, expfit.tau,expfit.C0, expfit.χ2_red  
        ))
        
        push!(trace_indices_i, i)
        push!(trace_indices_j, j)
        push!(original_data, y)
        push!(fitted_data, expfit.fitX)
    end

    return fit_results, trace_indices_i, trace_indices_j, original_data, fitted_data
end


"""
    fit_pelt(trace, cro, ic; si=1, sf=5)

Fit single trace using PELT changepoint detection.
Returns PeltFit or nothing if fit fails.
"""
function fit_pelt(trace::Vector{T}, cro::Dict, ic::Int; si::Int=1, sf::Int=5) where {T<:Real}
    # get steps from PELT for this ic
    pfit, pen, cost, ncp = steps_from_pelt(trace, cro, ic)
    
    if pen <= 0
        return nothing
    end
    
    if pen == 0 && cost == 0 && ncp == 0
        return nothing
    end
    
    # Protection: handle small length differences between trace and fitted steps
    len_diff = length(trace) - length(pfit)
    if abs(len_diff) == 1
        if len_diff > 0
            # trace is 1 element longer, extend pfit
            push!(pfit, pfit[end])
        else
            # pfit is 1 element longer, extend trace (shouldn't normally happen, but for safety)
            trace = vcat(trace, trace[end])
        end
    elseif len_diff != 0
        @warn "Length mismatch between trace ($(length(trace))) and fitted steps ($(length(pfit))): difference = $len_diff"
    end
    
    c1 = cro["changepoints"][ic]
    nsteps, sth, stt, stl = getsteps_pelt(trace, cro, ic)
    
    # Calculate reduced chi-squared
    χ2_red = chi2_dof(trace, pfit, nsteps)
    
    # Check dof
    if χ2_red === nothing
        @warn "no degrees of freedom"
        return nothing
    end
    
    return PeltFit(trace,      # dataX
                   pfit,       # fitX
                   cost,       # cost
                   pen,        # penalty
                   nsteps,     # nstep
                   sth,        # stepHeight
                   stt,        # stepTime
                   stl,        # stepLength
                   si,         # si
                   sf,         # sf
                   χ2_red)     # χ2_red
end


"""
    fit_pelt(trzs, df; ic=1, si=1, sf=5, auto_noise_region=true)

Fit multiple traces from sparse matrix using PELT.
Returns DataFrame with results and trace data arrays.
"""
function fit_pelt(trzs::SparseMatrixCSC{Vector{T}, Int}, df::DataFrame; 
                  ic=1, si=1, sf=5, auto_noise_region=true) where {T<:Real}
    I = Int[]
    J = Int[]
    DX = Vector{Vector{T}}()
    FX = Vector{Vector{T}}()
    
    ng = 0
    
    df2 = DataFrame(i=Int[], j=Int[], nmol=Int[], cost=Float32[], penalty=Float32[], 
                    nstep=Int[], stepHeight=Float32[], 
                    stepTime=Int[], stepLength=Int[], chi2=Float32[])
    
    for row in eachrow(df)
        debug("fitting trace ($(row.i),$(row.j))")
        
        tz = trzs[row.i, row.j] 
        
        
            cro, actual_i, actual_f = try
                pelt_fit(tz; i=si, f=sf, auto_noise_region=auto_noise_region)
        catch e
            debug("PELT fitting failed with error: $e")
            continue
        end

        if length(cro["changepoints"]) < ic
            ic = length(cro["changepoints"])
        end

        #println("number = ", ng, "cro = ", length(cro["changepoints"][ic]), " i = ", actual_i, " sf = ", actual_f)
        #println("trz = ", length(tz))
        pfit = fit_pelt(tz, cro, ic; si=actual_i, sf=actual_f)
        
        if pfit !== nothing
            ng += 1
            push!(I, row.i)
            push!(J, row.j)
            push!(DX, pfit.dataX)
            push!(FX, pfit.fitX)
            
            for k in 1:pfit.nstep
                push!(df2, (row.i, row.j, ng, pfit.cost, pfit.penalty, pfit.nstep, 
                           pfit.stepHeight[k], pfit.stepTime[k], pfit.stepLength[k], pfit.χ2_red))
            end
        end
    end
    
    return df2, I, J, DX, FX
end


"""
    fitpelt(tz; max_ic=5, si=1, sf=5, auto_noise_region=true, tolerance=1e-4)

Optimize PELT fitting across multiple ic values.
Returns number of changepoint sets and vector of PeltFit results.
"""
function fitpelt(tz::Vector{T}; max_ic::Int=5, si::Int=1, sf::Int=5,
                 auto_noise_region=true, tolerance::Float64=1e-4) where {T<:Real}
    
    PFIT = Vector{PeltFit{T}}()
    
    cro, actual_i, actual_f = try
        pelt_fit(tz; i=si, f=sf, auto_noise_region=auto_noise_region)
    catch e
        debug("PELT fitting failed with error: $e")
        return nothing, Vector{PeltFit{T}}()
    end
    
    if cro === nothing
        return nothing, Vector{PeltFit{T}}()
    end
    
    icx = length(cro["changepoints"])
    icm = min(icx, max_ic)
    
    for ic in 1:icm
        pfit = fit_pelt(tz, cro, ic; si=actual_i, sf=actual_f)
        if pfit !== nothing
            # Protection: handle small length differences between trace and fitted steps  
            len_diff = length(tz) - length(pfit.fitX)
            if abs(len_diff) == 1
                if len_diff > 0
                    # tz is 1 element longer, extend pfit.fitX
                    pfit = PeltFit(pfit.dataX, vcat(pfit.fitX, pfit.fitX[end]), pfit.cost, pfit.penalty, 
                                  pfit.nstep, pfit.stepHeight, pfit.stepTime, pfit.stepLength, 
                                  pfit.si, pfit.sf, pfit.χ2_red)
                else
                    # pfit.fitX is 1 element longer, extend tz (shouldn't normally happen, but for safety)
                    tz = vcat(tz, tz[end])
                end
            elseif len_diff != 0
                @warn "Length mismatch between trace ($(length(tz))) and fitted steps ($(length(pfit.fitX))): difference = $len_diff"
            end
            push!(PFIT, pfit)
        end
    end
    
    return icx, PFIT
end


"""
    fit_asf(tz; niter=10, tresH=0.15)

Fit single trace using Automatic Step Finder (ASF).
Returns AsfFit or nothing if fit fails.
"""
function fit_asf(tz::Vector{T}; niter::Int=10, tresH::Real=0.15) where {T<:Real}
    # Call fit_traces with sel="core" (fixed)
    dataX, FitX, S_curve, best_shot, iter, _ = fit_traces(tz; niter=niter, tresH=tresH, sel="core")
    
    # Check if fit was successful
    if best_shot <= 0
        return nothing
    end
    
    # Extract step information using getsteps
    sth, stt, stl = getsteps(FitX)
    
    # Calculate reduced chi-squared
    χ2_red = chi2_dof(dataX, FitX, length(sth))
    
    # Check dof
    if χ2_red === nothing
        @warn "no degrees of freedom"
        return nothing
    end

    # Return structured result (match AsfFit field order)
    return AsfFit(
        dataX,      # dataX: original trace  
        FitX,       # fitX: fitted values
        S_curve,    # S_curve: optimization curve
        best_shot,  # best_shot: quality measure
        iter,       # iter: iterations performed
        stl,        # steps: step lengths (durations)
        sth,        # stepHeights: step heights
        stt,        # stepTimes: step start times
        χ2_red      # χ2_red: reduced chi-squared
    )
end

"""
    elbow(x::AbstractVector, y::AbstractVector; method=:kneedle,
          smooth_window::Int=5, smooth_passes::Int=1, return_index::Bool=false)

Finds the elbow (knee) point of a concave curve defined by (x, y).

Methods:
- `:kneedle`     — point of maximum perpendicular distance to the chord (first→last).
- `:curvature`   — point of maximum discrete curvature after smoothing.

Options:
- `smooth_window` (odd ≥ 3): moving-average window used for optional smoothing (curvature only by default).
- `smooth_passes`: how many smoothing passes to apply.

Returns: (x_elbow, y_elbow) by default, or index if `return_index=true`.

Notes:
- Input is sorted by x; NaNs are dropped; duplicates in x are collapsed by last occurrence.
- For noisy data, prefer `method=:curvature` or increase smoothing.
"""
function elbow(x::AbstractVector, y::AbstractVector; method=:kneedle,
               smooth_window::Int=5, smooth_passes::Int=1, return_index::Bool=true)
    @assert length(x) == length(y) "x and y must have same length"
    n = length(x)
    @assert n ≥ 3 "need at least 3 points"

    # Copy & clean
    xs = collect(x); ys = collect(y)

    # Drop NaNs
    mask = .!(isnan.(xs) .| isnan.(ys))
    xs, ys = xs[mask], ys[mask]
    @assert length(xs) ≥ 3 "need at least 3 valid points"

    # Sort by x, keep last when x duplicates
    p = sortperm(xs, alg=QuickSort)
    xs, ys = xs[p], ys[p]
    # Collapse duplicate x by keeping last
    keep = trues(length(xs))
    for i in 1:length(xs)-1
        if xs[i] == xs[i+1]
            keep[i] = false
        end
    end
    xs, ys = xs[keep], ys[keep]
    n = length(xs)
    @assert n ≥ 3 "need ≥3 distinct x-values"

    if method === :kneedle
        # Normalize to [0,1] for scale invariance
        x0, x1 = xs[1], xs[end]
        y0, y1 = ys[1], ys[end]
        xr = (xs .- x0) ./ (x1 - x0 + eps())  # avoid div0 if degenerate
        yr = (ys .- y0) ./ (y1 - y0 + eps())

        # Line coefficients for chord in (xr,yr): ax + by + c = 0
        dx = xr[end] - xr[1]; dy = yr[end] - yr[1]
        # Distance from point to line: |dy*x - dx*y + (xr[end]*yr[1] - yr[end]*xr[1])| / √(dx^2+dy^2)
        denom = hypot(dx, dy) + eps()
        c = xr[end]*yr[1] - yr[end]*xr[1]
        d = @. abs(dy*xr - dx*yr + c) / denom

        # exclude endpoints (they're zero by construction)
        i = argmax(d)
        (i == 1 || i == n) && (i = argmax(@view d[2:n-1]) + 1)

        return return_index ? i : (xs[i], ys[i])

    elseif method === :curvature
        # Light smoothing to stabilize numerical derivatives
        xsm, ysm = xs, ys
        if smooth_window ≥ 3
            @assert isodd(smooth_window) "smooth_window must be odd"
            for _ in 1:max(1, smooth_passes)
                xsm = movavg(xsm, smooth_window)
                ysm = movavg(ysm, smooth_window)
            end
        end

        # First and second finite differences (non-uniform spacing aware)
        # Curvature κ ≈ |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
        κ = zeros(Float64, n)
        @inbounds for k in 2:n-1
            x₋, x₀, x₊ = xsm[k-1], xsm[k], xsm[k+1]
            y₋, y₀, y₊ = ysm[k-1], ysm[k], ysm[k+1]

            h₁ = x₀ - x₋
            h₂ = x₊ - x₀
            # derivatives via three-point unequal spacing formulas
            x′ = 1.0            # dx/dx
            y′ = (-(h₂/(h₁*(h₁+h₂)))*y₋ + ((h₂ - h₁)/(h₁*h₂))*y₀ + (h₁/(h₂*(h₁+h₂)))*y₊)
            # second derivatives
            x″ = 0.0
            y″ = (2/(h₁*(h₁+h₂)))*y₋ - (2/(h₁*h₂))*y₀ + (2/(h₂*(h₁+h₂)))*y₊

            denom = (x′^2 + y′^2)^(3/2) + eps()
            κ[k] = abs(x′*y″ - y′*x″) / denom
        end

        i = argmax(@view κ[2:n-1]) + 1
        return return_index ? i : (xs[i], ys[i])

    else
        Base.error("unknown method = $method (use :kneedle or :curvature)")
    end
end

# simple centered moving average with edge replication
function movavg(v::AbstractVector, w::Int)
    n = length(v); r = (w - 1) >> 1
    out = similar(v, Float64)
    @inbounds for i in 1:n
        lo = max(1, i - r); hi = min(n, i + r)
        out[i] = sum(@view v[lo:hi]) / (hi - lo + 1)
    end
    out
end