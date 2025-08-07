"""
Analysis Functions for Image Processing and Peak Detection

This module contains functions for analyzing fluorescence microscopy data,
including peak detection and local maxima identification.
"""

using DataFrames
using Images
using SparseArrays
using LsqFit
using Statistics

export detect_local_maxima, fit_traces_exponential

"""
    detect_local_maxima(frame::AbstractMatrix{<:Real}; 
                        threshold=0.0, dx=0, dy=0, 
                        window_size=(3,3)) -> DataFrame

Detect local maxima in a 2D image frame using a sliding window approach.

# Arguments
- `frame`: 2D image matrix
- `threshold`: Minimum intensity a peak must exceed (default: 0.0)
- `dx`: Margin to exclude on left and right edges (default: 0)
- `dy`: Margin to exclude on top and bottom edges (default: 0)
- `window_size`: Size of sliding window for maxima detection (default: (3,3))

# Returns
DataFrame with columns: i (row), j (column), intensity

# Notes
- Uses replicate padding to handle borders during window operations
- Excludes peaks within dx/dy pixels of frame edges
- Window size must be odd in both dimensions

# Examples
```julia
# Create a test image with peaks
frame = [1.0 2.0 1.0; 2.0 5.0 2.0; 1.0 2.0 1.0]
peaks = detect_local_maxima(frame; threshold=3.0)

# Use custom window size
peaks = detect_local_maxima(frame; threshold=2.0, window_size=(5,5))
```
"""
function detect_local_maxima(frame::AbstractMatrix{<:Real}; 
                            threshold=0.0, dx=0, dy=0, 
                            window_size=(3,3))
    
    # Validate window size
    if any(iseven.(window_size))
        throw(ArgumentError("Window size must be odd in both dimensions"))
    end
    
    # Calculate center index for the window
    center_idx = div(prod(window_size) + 1, 2)
    
    # Find local maxima using sliding window
    is_max = mapwindow(x -> x[center_idx] == maximum(x), frame, window_size; 
                       border=Pad(:replicate))
    
    # Apply threshold and find candidates
    candidates = findall(is_max .& (frame .> threshold))
    
    # Pre-filter candidates based on edge constraints
    height, width = size(frame)
    valid_candidates = filter(candidates) do I
        i, j = Tuple(I)
        return i > dy && j > dx && i <= height - dy && j <= width - dx
    end
    
    # Extract coordinates and intensities
    n_peaks = length(valid_candidates)
    i_vals = Vector{Int}(undef, n_peaks)
    j_vals = Vector{Int}(undef, n_peaks)
    intensities = Vector{Float64}(undef, n_peaks)
    
    for (idx, I) in enumerate(valid_candidates)
        i, j = Tuple(I)
        i_vals[idx] = i
        j_vals[idx] = j
        intensities[idx] = frame[I]
    end
    
    return DataFrame(i = i_vals, j = j_vals, intensity = intensities)
end

"""
    fit_traces_exponential(TRZS::SparseMatrixCSC{Vector{T}, Int}, 
                          df::DataFrame; 
                          initial_params=nothing, 
                          max_iterations=1000) -> Tuple

Fit exponential decay models to fluorescence traces extracted from image data.

# Arguments
- `TRZS`: Sparse matrix containing trace data as vectors, indexed by (i,j) coordinates
- `df`: DataFrame containing trace coordinates with columns `i` and `j`
- `initial_params`: Optional tuple (A0, tau, C0) for initial parameter guesses. If `nothing`, 
  parameters are estimated automatically from each trace
- `max_iterations`: Maximum number of iterations for curve fitting (default: 1000)

# Returns
A tuple containing:
- `fit_results::DataFrame`: Fitting results with columns:
  - `i`, `j`: Trace coordinates  
  - `A0`: Amplitude parameter
  - `tau`: Decay time constant
  - `C0`: Offset parameter
  - `chi2`: Chi-squared value
  - `chi2_red`: Reduced chi-squared value
  - `n`: Number of data points in trace
  - `converged`: Whether fit converged successfully
- `trace_indices_i::Vector{Int}`: Row indices for successful fits
- `trace_indices_j::Vector{Int}`: Column indices for successful fits  
- `original_data::Vector{Vector{T}}`: Original trace data for successful fits
- `fitted_data::Vector{Vector{T}}`: Fitted exponential curves

# Model
Fits the exponential decay model: `f(t) = A0 * exp(-t/tau) + C0`

# Notes
- Traces that fail to converge are excluded from the results
- Automatic parameter estimation uses: A0 = max(y), tau = 10.0, C0 = min(y)
- Chi-squared values help assess goodness of fit
- All successful fits are guaranteed to have converged

# Examples
```julia
# Basic usage with automatic parameter estimation
results, i_vals, j_vals, orig_data, fit_data = fit_traces_exponential(traces, peaks)

# With custom initial parameters
results, i_vals, j_vals, orig_data, fit_data = fit_traces_exponential(
    traces, peaks; initial_params=(100.0, 15.0, 5.0)
)
```
"""
function fit_traces_exponential(
    TRZS::SparseMatrixCSC{Vector{T}, Int}, 
    df::DataFrame;
    initial_params=nothing,
    max_iterations::Int=1000
) where {T<:Real}

    # Pre-allocate result arrays
    n_traces = nrow(df)
    fit_results = DataFrame(
        i = Int[], j = Int[],
        A0 = Float64[], tau = Float64[], C0 = Float64[],
        chi2 = Float64[], chi2_red = Float64[], 
        n = Int[], converged = Bool[]
    )
    
    trace_indices_i = Int[]
    trace_indices_j = Int[]  
    original_data = Vector{Vector{T}}()
    fitted_data = Vector{Vector{T}}()

    # Define exponential decay model
    exponential_model(t, p) = p[1] .* exp.(-t ./ p[2]) .+ p[3]

    for row in eachrow(df)
        i, j = row.i, row.j
        y = TRZS[i, j]
        
        # Skip empty traces
        if isempty(y)
            @warn "Empty trace at ($i, $j), skipping"
            continue
        end
        
        t = collect(1:length(y))

        # Determine initial parameters
        if initial_params !== nothing
            p0 = collect(Float64.(initial_params))
        else
            # Automatic parameter estimation with bounds checking
            y_max = maximum(y)
            y_min = minimum(y) 
            p0 = [Float64(y_max - y_min), 10.0, Float64(y_min)]
        end
        
        # Ensure positive initial parameters where needed
        p0[1] = max(p0[1], 1e-6)  # Amplitude must be positive
        p0[2] = max(p0[2], 1e-6)  # Tau must be positive
        
        # Attempt curve fitting with error handling
        fit_result = try
            curve_fit(exponential_model, Float64.(t), Float64.(y), p0; 
                     maxIter=max_iterations)
        catch err
            @warn "Fit failed for trace ($i, $j): $(typeof(err))" 
            continue
        end

        # Check convergence
        if !fit_result.converged
            @warn "Fit did not converge for trace ($i, $j)"
            continue
        end

        # Calculate fitted values and statistics
        fitted_vals = exponential_model(Float64.(t), fit_result.param)
        residuals = Float64.(y) .- fitted_vals
        χ2 = sum(residuals .^ 2)
        dof = length(y) - length(fit_result.param)
        χ2_red = dof > 0 ? χ2 / dof : Inf

        # Store successful fit
        push!(fit_results, (
            i, j, 
            fit_result.param[1], fit_result.param[2], fit_result.param[3],
            χ2, χ2_red, length(y), true
        ))
        
        push!(trace_indices_i, i)
        push!(trace_indices_j, j)
        push!(original_data, y)
        push!(fitted_data, fitted_vals)
    end

    return fit_results, trace_indices_i, trace_indices_j, original_data, fitted_data
end