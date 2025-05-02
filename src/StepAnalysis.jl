module StepAnalysis
using Revise
using StatsPlots
using Statistics
using StatsBase
import Measures
using SparseArrays
using DataFrames

using NPZ

export get_stats, pixel_trace, traces_stats, traces_h2d, traces_thr
# get_traces
export plot_traces_stats, plot_stats, plot_traces_h2d
export plot_steps, plot_trx, fit_data, fit_data3, fit_traces, getsteps, plotfit, plotfit2
export traces_above_thr, build_traces
export plot_trace, plot_traces, plot_sdf, plot_frames
export renumber_nmol!

# plt_ft_trace

include("JStepFinder.jl")
using .JStepFinder

include("histos.jl")
using .histos

include("SimpleLogger.jl")
using .SimpleLogger
global_log_level = ERROR


"""
    fit_traces(trz::AbstractVector{<:AbstractVector{T}}; 
               n::Int=1, niter::Int=10) where {T<:Real}

Pick out the n-th trace from any 1-D collection of 1-D traces and delegate
to the single-trace `fit_traces(::AbstractVector{T},…)` method.

# Arguments
- `trz`: any vector-like container of vector-like traces, e.g. `Vector{Vector{Float64}}`,
  `Vector{SubArray{Float32,1,…}}`, etc.
- `n`: index of the trace to fit (default 1)
- `niter`: number of iterations for the step finder (default 10)

# Returns
whatever your existing `fit_traces(dataX::AbstractVector{T},…)` returns.
"""
function fit_traces(trz::AbstractVector{<:AbstractVector{T}}; 
                    n::Int=1, niter=10, thr=0.0) where {T<:Real}
    dataX = trz[n]
    return fit_traces(dataX; niter=niter, thr=thr)
end


"""
    fit_traces(trz::AbstractMatrix{<:AbstractVector{T}}, i::Int, j::Int; niter::Int=10) where {T<:Real}

Extract the 1-D trace stored at `(i,j)` in a 2-D array of vectors and run
the same step-fitting routine on it.

# Arguments
- `trz`: A 2-D array whose each entry `trz[i,j]` is itself a 1-D vector of type `T<:Real`.
- `i, j`: row and column indices selecting which trace to fit.
- `niter`: number of iterations for the underlying fit (default 10).

# Returns
Whatever your `fit_traces(dataX::AbstractVector{T}, niter)` returns:
`(dataX, FitX0, FitX, S_curve, thr, best_shot)`.
"""
function fit_traces(trz::AbstractMatrix{<:AbstractVector{T}}, 
	                i::Int, j::Int; niter=10, thr=0.0, sel="auto") where {T<:Real}
    # pick out the vector at (i,j)
    dataX = trz[i, j]
    # delegate to the single‐vector fit_traces
    return fit_traces(dataX; niter=niter, tresH=thr, sel=sel)
end


"""
    fit_traces(dataX::AbstractVector{T}, niter::Int) where {T<:Real}

Perform a step‐fit analysis on a single 1D signal, exploring a range of thresholds
to select the optimal step detection parameters.

# Arguments
- `dataX::AbstractVector{T}`  
  The input trace (time series) of type `T<:Real`, e.g. `Vector{Float64}`.
- `niter::Int`  
  Number of iterations the underlying `stepfindcore` algorithm should perform.

# Returns
A 6‐tuple `(dataX, FitX0, FitX, S_curve, thr, best_shot)` where:
- `dataX`    — the original input trace.
- `FitX0`    — the step fit result at zero threshold (baseline fit).
- `FitX`     — the step fit result at the chosen optimal threshold.
- `S_curve`  — the S‐curve values recorded during the final fit.
- `thr`      — the threshold value determined to be best.
- `best_shot`— an integer flag indicating fit quality (0 means fit failed).

# Behavior
1. Runs a zero‐threshold fit (`tresH=0.0`) to produce `FitX0`.  
2. Iterates thresholds in `range(0.05, 0.5; length=30)`, stopping when the fit fails (`best_shot==0`),  
   then uses the previous threshold as `thr`.  
3. If a nonzero `thr` was found, reruns the fit at `thr` to produce final `FitX` and `S_curve`.  
"""
function fit_traces(dataX::AbstractVector{T}; niter::Int=10, tresH::Real=0.15,
	                tol::Real=1e-3, max_iter::Int=5, sel="auto") where {T<:Real}
    
    # baseline fit at zero threshold



	if sel == "auto"
		S_curve, best_shot, FitX, iter, cc = auto_step_main(dataX; tresH=tresH, N_iter=niter, 
                        tol=tol, max_iter=max_iter)
	else
		FitX, _, _, S_curve, best_shot = stepfindcore(dataX; tresH=tresH, N_iter=niter)
		iter = -10
		cc = -10.0
	end

    return dataX, FitX, S_curve, best_shot, iter, cc
end


"""
    fit_data(trz::AbstractVector{<:AbstractVector{T}}; niter=10, thr=0.0, rmax=1000, fplot=50) where {T<:Real}

Fits step functions to a vector-like collection of time traces using the `stepfindcore` algorithm.

# Arguments
- `trz`: Any vector-like container of vector-like traces, e.g. `Vector{Vector{Float64}}`, `Vector{SubArray{Float32,1,…}}`, etc.
- `niter`: Number of iterations used in the step-fitting algorithm (default: 10).
- `thr`: Threshold value for step detection (default: 0.0).
- `rmax`: Maximum number of traces to process (default: 1000).
- `fplot`: Frequency of traces to store for plotting (every `fplot`-th trace is stored, default: 50).

# Returns
- `DX`: Vector of raw traces selected for plotting.
- `FX`: Vector of fitted traces corresponding to `DX`.
- `TL`: Vector containing the number of steps found in each trace.
- `StepHeight`: Concatenated step heights from all traces (only if more than one step is found).
- `StepTime`: Concatenated step times from all traces.
- `StepLength`: Concatenated lengths of the segments between steps.
"""
function fit_data(trz::AbstractVector{<:AbstractVector{T}}; 
	niter=10, thr=0.0, rmax=1000, fplot=50) where {T<:Real}

	DX = Vector{Vector{T}}()
	FX = Vector{Vector{T}}()
	StepHeight = T[]
	StepTime = Int[]
	StepLength = Int[]
	TL = Int[]

	nmax = min(rmax, length(trz))

	for n in 1:nmax
		dataX = trz[n]
		FitX, _, _, S_curve, best_shot0 = stepfindcore(dataX; 
												tresH=thr, N_iter=niter, 
												demo=false)

		if mod(n, fplot) == 0
			push!(DX, dataX)
			push!(FX, FitX)
		end

		sth, stt, stl = getsteps(FitX)
		push!(TL, length(sth))

		if length(sth) > 1 
			append!(StepHeight, sth)
			append!(StepTime, stt)
			append!(StepLength, stl)
		end
	end

	return DX, FX, TL, StepHeight, StepTime, StepLength
end


"""
    fit_data3(trz::AbstractMatrix{<:AbstractVector{T}}; niter=10, thr=0.0, fplot=50) where {T<:Real}

Analyze a 2D matrix of time traces by performing step-fitting on each pixel trace and collecting trace and step statistics.

# Arguments
- `trz`: A 2D matrix where each element is a 1D time trace (i.e., a vector of type `T<:Real`).
- `niter`: Number of iterations used in the step-fitting algorithm (default: 10).
- `thr`: Threshold value for step detection (default: 0.0).
- `fplot`: Interval for saving fitted traces for plotting (default: every 50th row).

# Returns
- `DX`: Vector of raw traces sampled every `fplot` rows.
- `FX`: Vector of fitted traces corresponding to `DX`.
- `TL`: Vector containing number of steps for each plotted trace.
- `df`: A `DataFrame` with columns:
    - `nmol`: a unique integer index per trace,
    - `nstep`: number of steps in the trace,
    - `stepHeight`: height of each step,
    - `stepTime`: starting index of each step,
    - `stepLength`: duration of each step.

- ped is the pedestal (1600 for a 4 x 4 binnin)
- tt is the time treshold (tt = 1 eliminates the first bin, typically exp decay)
"""
function fit_data3(trz::AbstractMatrix{<:AbstractVector{T}}; 
                   niter=10, thr=0.0, ped=1600.0, tt=1.0, stdf=0.1, sel="auto") where {T<:Real}

	I = Int[]
	J = Int[]
    DX = Vector{Vector{T}}()
    FX = Vector{Vector{T}}()
    ITER = Int[]
	CC = Float64[]

    df = DataFrame(i=Int[], j=Int[], nmol=Int[], nstep=Int[],
                   stepHeight=T[], stepTime=Int[], stepLength=Int[])

    n, m = size(trz)
    nc = 0  # unique molecule index
	nf = 0; # failed fit
	ng = 0; # good fit
    for i in 1:n
        for j in 1:m
            nc += 1
            dataX, FitX, S_curve, best_shot, iter, cc = fit_traces(trz, i, j; 
                                                         niter=niter,
														 thr=thr, sel=sel)

            if best_shot > 0 && std(FitX) > stdf
                sth, stt, stl = getsteps(FitX)
                nsteps = length(sth)
				ng+=1

				push!(I,i)
				push!(J,j)
				push!(DX, dataX)
				push!(FX, FitX)
				push!(ITER, iter)
				push!(CC, cc)

                #if mod(nc, fplot) == 0
                    
                #end

                for k in 1:nsteps
					if stt[k] > tt
                    	push!(df, (i, j, ng, nsteps, sth[k] - ped, stt[k], stl[k]))
					end
                end
            else
				nf += 1
                warn("Fit to trace in pixel ($(i), $(j)) failed bshot = $(best_shot), thr = $(thr)")
            end
        end
    end

	MDX = SparseArrays.sparse(I, J, DX, n, m)
	MFX = SparseArrays.sparse(I, J, FX, n, m)

	df, nc, ng, nf, I, J, MDX, MFX, ITER, CC

    #return DX, FX, TL, df, nc, nf
end


"""
    getsteps(arr::AbstractVector{T}; atol=1e-8) where {T<:Real}

Identifies constant segments in a piecewise constant signal and extracts their values, start indices, and lengths.

# Arguments
- `arr`: Input signal as a 1D array of real values.
- `atol`: Absolute tolerance for detecting change points in the signal (default: `1e-8`).

# Returns
- `constant_values`: Vector containing the constant value of each segment.
- `segment_starts`: Vector of indices marking the beginning of each constant segment.
- `segment_lengths`: Vector of segment lengths (number of elements in each constant region).

# Description
This function is typically used on step-like signals, where it detects transitions by comparing adjacent values within a specified tolerance. It is useful for analyzing step fits in time traces.
"""
function getsteps(arr::AbstractVector{T}; atol=1e-8) where {T<:Real}
    # Find where the difference between consecutive elements is not close to zero
    change_mask = .!isapprox.(diff(arr), 0.0, atol=atol)
    change_indices = findall(identity, change_mask) .+ 1

    # Get the constant values of each segment
    segment_values = arr[change_indices]
    constant_values = vcat(arr[1], segment_values)

    # Compute segment start and end indices
    segment_starts = vcat(1, change_indices)
    segment_ends = vcat(change_indices, length(arr) + 1)
    segment_lengths = segment_ends .- segment_starts

    return constant_values, segment_starts, segment_lengths
end


"""
Computes the gloabl statistics of the stack. Sum, mean and stds. Also creates
a vector of means and stds.
imst is the stack  
"""
function get_stats(imst::Array{T, 3}; bck=200.0) where {T<:AbstractFloat}
	
	n_frames = size(imst, 3)
	
	# Preallocate vectors
	total_intensity = zeros(T, n_frames)
	mean_intensity = zeros(T, n_frames)
	std_intensity = zeros(T, n_frames)

	
	# Loop over each image in the stack
	for i in 1:n_frames
	    frame = imst[:, :, i] .- bck
	    total_intensity[i] = sum(frame)
	    mean_intensity[i] = mean(frame)
	    std_intensity[i] = std(frame)
	end
	total_intensity, mean_intensity, std_intensity 
end


"""
Get the time trace for a given pixel 
"""
function pixel_trace(imstack::Array{T, 3}, i::Int, j::Int; bck=200.0) where {T<:AbstractFloat}
    # Get the number of frames (3rd dimension)
    n_frames = size(imstack, 3)
    
    # Extract intensity at (i,j) across all frames
    trace = [imstack[i, j, t] .- bck for t in 1:n_frames]
    
    return trace
end

"""
Trace statistics
"""
function traces_stats(imst::Array{T, 3}) where {T<:AbstractFloat}
	imx = size(imst)[1]
	jmx = size(imst)[2]
	xsum = zeros(imx,jmx)
	xmean = zeros(imx,jmx)
	xstd = zeros(imx,jmx)
	for i in  1:imx
		for j in  1:jmx
			trace = pixel_trace(imst, i, j)
			xsum[i,j] = sum(trace)
			xmean[i,j] = mean(trace)
			xstd[i,j] =  std(trace)
		end
	end
	xsum, xmean, xstd
end


function traces_h2d(xsum, xmean, xstd; bins=50)
	function h2d(x1, x2)
		# Ensure matrices have the same shape
		@assert size(x1) == size(x2)
	
		# Flatten both to vectors
		x = vec(x1)
		y = vec(x2)

		# Compute the 2D histogram with 50×50 bins
		h = fit(Histogram, (x, y), nbins=(bins, bins))

	end
	h2d(xsum, xstd), h2d(xmean, xstd)
end

"""
    build_traces(imst::AbstractArray{T,3}; nsigma::Real=10) where {T<:Real}

Given a 3D image stack `imst` of size (n × m × t), this function:

1. Computes the mean `μ` and standard deviation `σ` of the **first** frame `imst[:, :, 1]`.
2. Marks as “noisy” any pixel whose first‐frame intensity exceeds `μ + nsigma*σ`.
3. Builds and returns a matrix `TRZ` of size (n × m), where each entry is a `Vector{T}` of length `t`:
   - For a **good** pixel, the vector is that pixel’s full time‑trace `imst[i, j, 1:t]`.
   - For a **noisy** pixel, the vector is all zeros.

# Arguments
- `imst::AbstractArray{T,3}`  
  The input stack, with dimensions (rows, cols, frames).
- `nsigma::Real=10`  
  Threshold multiplier: pixels above `μ + nsigma*σ` in frame 1 are masked.

# Returns
- `TRZ::Matrix{Vector{T}}`  
  An n×m matrix of vectors, each of length t.

"""
# function build_traces(imst::AbstractArray{T,3}; nsigma::Real=10) where {T<:Real}

	
# 	n, m, t = size(imst)
# 	# 1) Compute threshold from first frame
	
# 	frame1 = imst[:, :, 1]
# 	μ = mean(frame1)
# 	σ = std(frame1)
# 	threshold = μ + nsigma * σ
	
# 	# 2) Build a boolean mask of "noisy" pixels
# 	#    mask[i,j] == true  ⇔ pixel (i,j) is to be zeroed
# 	mask = frame1 .> threshold
	
# 	# 3) Prepare a zero‐trace template
# 	zero_trace = zeros(T, t)
	
# 	# 4) Allocate output matrix of vectors
# 	TRZ = Matrix{Vector{T}}(undef, n, m)
	
# 	# 5) Fill in each pixel's trace or zeros
# 	for i in 1:n, j in 1:m
# 	    if mask[i, j]
# 	        # masked pixel → all‐zero trace
# 	        TRZ[i, j] = copy(zero_trace)
# 	    else
# 	        # good pixel → extract its time series
# 	        TRZ[i, j] = vec(imst[i, j, :])
# 	    end
# 	end

# 	return mask, TRZ
# end



function build_traces(imst::AbstractArray{T,3}) where {T<:Real}
	n, m, t = size(imst)
	
	# Allocate output matrix of vectors
	TRZ = Matrix{Vector{T}}(undef, n, m)
	
	# Fill in each pixel's trace 
	for i in 1:n, j in 1:m
	    TRZ[i, j] = vec(imst[i, j, :])
	end
	
	return TRZ
end


"""
Obtain traces above threshold 
"""
function traces_thr(imst::Array{T, 3}; xmthr=0.0, xsthr=10.0) where {T<:AbstractFloat}
	imx = size(imst)[1]
	jmx = size(imst)[2]

	M = Matrix{Vector{T}}(undef, imx, jmx) 
	for i in  1:imx
		for j in  1:jmx
			trace = pixel_trace(imst, i, j)
			xmean = mean(trace)
			xstd =  std(trace)
			if xmean > xmthr && xstd > xsthr
				M[i,j] = trace
			else
				M[i,j] = [0.0]
			end
		end
	end
	[v for v in M if v != [0.0]]
end


"""
    traces_above_thr(imst::Array{T, 3}; nsigma=3.0) where {T<:AbstractFloat}

Extracts pixel traces from a 3D image stack (time-series of images) for which the pixel value 
at the first frame exceeds a threshold, and stores them in a sparse matrix.

# Arguments
- `imst`: A 3D array (`imx × jmx × t`) of image data over time.
- `nsigma`: Threshold in units of standard deviations above the mean of the first frame (default: 3.0).

# Returns
- A sparse matrix of size `imx × jmx`, where each nonzero element is a `Vector{T}` representing the 
  time trace of a pixel that exceeds the intensity threshold at time `t = 1`.

# Description
The function computes the mean and standard deviation of the first frame of the image stack.
It selects all pixels where the intensity is greater than `mean + nsigma * std` and stores 
their full time traces in a sparse matrix. This is useful for efficiently handling sparse 
events in large image sequences, such as bright spots in microscopy images.

"""
function traces_above_thr(imst::Array{T, 3}; nsigma=3.0) where {T<:AbstractFloat}
	imx = size(imst)[1]
	jmx = size(imst)[2]

	I = Int[]
	J = Int[]
	V =Vector{T}[]
	
	xm = mean(imst[:,:,1])
	xstd = std(imst[:,:,1])
	xtrh = xm + nsigma * xstd
	for i in  1:imx
		for j in  1:jmx
			if imst[i,j,1] >= xtrh
				push!(I,i)
				push!(J,j)
				push!(V,pixel_trace(imst, i, j))
			end
		end
	end
	sparse(I, J, V, imx, jmx)
	
end


#####
## Plots
#####

function plot_frames(imst::AbstractArray{T,3}; nscale::Int=20) where {T<:Real}

	FF = []

	
    for i in 1:9
        fn = (i-1) * nscale + i
		if fn > size(imst)[3]
			warn("requested frame = $(fn) is to large, set smaller nscale")
			fn = size(imst)[3] -i
			warn("set fn = $(fn)")
		end
        push!(FF, heatmap(imst[:, :, fn],
            colorbar=false,  # Optional: removes extra space
            title="Frame $fn",
            titlefontsize=7,
            tickfontsize=6,
            guidefontsize=6,
            titlelocation=:left,
            aspect_ratio=:equal))
    end

    plot(FF...;
        layout=(3, 3),
        size=(900, 900),
        margin=1.0*Measures.mm,
        top_margin=1.0*Measures.mm,
        bottom_margin=1.0*Measures.mm,
        left_margin=1.0*Measures.mm,
        right_margin=1.0*Measures.mm,
        plot_titlefontsize=7,
        legendfontsize=6)
end


function plot_frames_with_centroids(imst::AbstractArray{T,3}, region_stats::Vector{Vector{Dict{Symbol, Any}}}; 
	nscale::Int=20) where {T<:Real}

	FF = []

	
    for i in 1:9
        fn = (i-1) * nscale + i
		if fn > size(imst)[3]
			warn("requested frame = $(fn) is to large, set smaller nscale")
			fn = size(imst)[3] -i
			warn("set fn = $(fn)")
		end
		frame = imst[:, :, fn]
		centroids = [r[:centroid] for r in region_stats[fn]]
	
		p = heatmap(frame, color=:grays, aspect_ratio=1, title="Frame $fn",
		titlefontsize=7,
		tickfontsize=6,
		guidefontsize=6,
		titlelocation=:left,
		aspect_ratio=:equal)
		p = scatter!(p, [c[1] for c in centroids], [c[2] for c in centroids], markersize=4, color=:red)
        push!(FF, pp)
    end

    plot(FF...;
        layout=(3, 3),
        size=(900, 900),
        margin=1.0*Measures.mm,
        top_margin=1.0*Measures.mm,
        bottom_margin=1.0*Measures.mm,
        left_margin=1.0*Measures.mm,
        right_margin=1.0*Measures.mm,
        plot_titlefontsize=7,
        legendfontsize=6)
end



function plotfit2(dataX::AbstractVector{<:Real}, 
	FitX1::AbstractVector{<:Real}, 
	FitX2::AbstractVector{<:Real}; 
	figsize=(1000, 800))  # Size in pixels, similar to figsize in matplotlib

	plt1 = plot(1:length(dataX), dataX, 
	label="Noisy Signal", color=:gray, lw=1,
	xlabel="time steps", ylabel="Intensity", title="Trace", 
	legend=:topright, grid=true)

	plot!(plt1, 1:length(FitX1), FitX1, 
	label="Fit1", color=:blue, lw=2)
	plot!(plt1, 1:length(FitX2), FitX2, 
	label="FitX2 ", color=:red, lw=2)
	plt1
end


function plotfit(dataX::AbstractVector{<:Real}, 
	             FitX1::AbstractVector{<:Real}, 
	             S_curve::AbstractVector{<:Real}, i::Int, j::Int; 
	             figsize=(1000, 800))  # Size in pixels, similar to figsize in matplotlib

	plt1 = plot(1:length(dataX), dataX, 
	label="Noisy Signal", color=:gray, lw=1,
	xlabel="time steps", ylabel="Intensity", title="Trace in pixel ($(i), $(j))", 
	legend=:topright, grid=true)

	plot!(plt1, 1:length(FitX1), FitX1, 
	label="Step Fit", color=:red, lw=2)

	plt2 = plot(1:length(S_curve), S_curve, 
	marker=:circle, label="S-curve",
	xlabel="Iteration", ylabel="S-value (MSE ratio)", 
	title="Goodness of Fit (S-curve)", 
	grid=true, legend=:topright)

	plot(plt1, plt2, layout=(2, 1), size=figsize)
end


function plot_traces_stats(xsum, xmean, xstd; meanmx=25.0, summx=1e+3, stdmx=50.0)
	xmin = 0.0
	xmax = meanmx
	data = vec(xmean)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p1 = histogram(filtered, bins=collect(edges), title="Mean of traces", xlabel="trace index", ylabel="Frequency")

	xmin = 0.0
	xmax = summx
	data = vec(xsum)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p2 = histogram(filtered, bins=collect(edges), title="Sum of traces", xlabel="trace index", ylabel="Frequency")

	xmin = 0.0
	xmax = stdmx
	data = vec(xstd)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p3 = histogram(filtered, bins=collect(edges), title="std of traces", xlabel="trace index", ylabel="Frequency")

	plot(p1, p2, p3;
     layout=(3, 1),
     size=(700, 900),
     margin=5 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=10,
     guidefontsize=9,
     tickfontsize=8,
     legendfontsize=8)
end


function plot_stats(total_intensity, mean_intensity, std_intensity)
	# Plot mean and std over time
	n_frames = length(total_intensity)
	# Create a 3-row, 1-column plot layout
	plot_layout = @layout [a; b; c]

	p1 = plot(1:n_frames, total_intensity,
          label="Total Intensity", xlabel="Frame", ylabel="Total",
          title="Total Intensity Over Time", lw=1)

	p2 = plot(1:n_frames, mean_intensity,
          label="Mean Intensity", xlabel="Frame", ylabel="Mean",
          title="Mean Intensity Over Time", lw=1)

	p3 = plot(1:n_frames, std_intensity,
          label="Std Dev", xlabel="Frame", ylabel="Std",
          title="Standard Deviation Over Time", lw=1)

	# Combine them in a single figure
	#plot(p1, p2, p3, layout=plot_layout, size=(700, 900))
	plot(p1, p2, p3;
     layout=(3, 1),
     size=(700, 900),
     margin=5 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=10,
     guidefontsize=9,
     tickfontsize=8,
     legendfontsize=8)
end


function plot_traces_h2d(xsum, xmean, xstd; bins=50)
	function h2d(x1, x2, x1l, x2l, xt)
		# Ensure matrices have the same shape
		@assert size(x1) == size(x2)
	
		# Flatten both to vectors
		x = vec(x1)
		y = vec(x2)
	
		# 2D histogram
		histogram2d(x, y;
		       bins=bins,
		       xlabel=x1l,
		       ylabel=x2l,
		       title=xt)
	end
	p1 = h2d(xsum, xstd, "xsum", "xstd", "Sum vs Std")
	p2 = h2d(xmean, xstd, "mean", "xstd", "Mean vs Std")
	plot(p1, p2;
     layout=(1, 2),
     margin=1 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=8,
     guidefontsize=9,
     tickfontsize=6,
     legendfontsize=8)
end


function plot_steps(TL::Vector{Int}, 
	StepHeight::Vector{Float64},
	StepTime::Vector{Int}, 
	StepLength::Vector{Int};
	nbins = 50,
	shlims = (-750.0, 100.0),
	figsize=(1200, 1000))

	p1 = step_hist(TL, xlabel="# steps", ylabel="# of occurrences")
	p2 = step_hist(-1.0 .* StepHeight, 
	xlim=shlims, nbins=nbins,
	xlabel="Step height", 
	ylabel="# of occurrences", 
	logy=true)
	p3 = step_hist(StepTime, xlabel="Step times", ylabel="# of occurrences", logy=true)
	p4 = step_hist(StepLength, xlabel="Step lengths", ylabel="# of occurrences")

	plot(p1, p2, p3, p4, layout=(2, 2), size=figsize)
end


function plot_trx(DX::AbstractVector{<:AbstractVector{T}}, FX::AbstractVector{<:AbstractVector{T}}; nx=4, ny=4, figsize=(1500, 500)) where {T<:Real}

	nplots = min(length(DX), nx * ny)
	plots_array = Vector{Any}(undef, nx * ny)
	
	for i in 1:nplots
		plt = plot(1:length(DX[i]), DX[i],
			label="Noisy Signal", color=:gray, lw=1,
			xlabel="time steps", ylabel="Intensity", title="Trace", 
			legend=:topright, grid=true)

		plot!(plt, 1:length(FX[i]), FX[i], label="Step Fit thr=0", 
		color=:blue, lw=2)
		plots_array[i] = plt
	end

	# Fill unused plots if any
	for i in nplots+1:nx*ny
		plots_array[i] = plot(title="(empty)")
	end

	plot(plots_array...; layout=(nx, ny), size=figsize)
end

"""
    plot_trx(DX, FX; ni=1, nf=9, layout=(3,3), figsize=(1500, 500))

Plot a subset of raw vs fitted traces stored in two matrices of vectors.

# Arguments
- `DX::AbstractMatrix{<:AbstractVector{T}}`  
  Matrix (dense or sparse) whose nonzero entries are the raw 1-D traces.
- `FX::AbstractMatrix{<:AbstractVector{T}}`  
  Matrix (same size & sparsity pattern as `DX`) of fitted traces.
- `ni::Int=1`  
  Index (into the list of nonzero traces) of the first trace to plot.
- `nf::Int=9`  
  Index of the last trace to plot.
- `layout::Tuple{Int,Int}=(3,3)`  
  Number of rows and columns in the subplot grid.
- `figsize=(1500, 500)`  
  Figure size in pixels `(width, height)`.

# Returns
A `Plots.Plot` object with `(nf−ni+1)` subplots showing raw vs fit traces.

# Throws
- If `DX` and `FX` do not have identical dimensions or nonzero positions.  
- If `ni`/`nf` are out of the valid range.  
- If `layout` is too small to hold `(nf−ni+1)` subplots.
"""
function plot_trx(DX::AbstractMatrix{<:AbstractVector{T}},
                  FX::AbstractMatrix{<:AbstractVector{T}};
                  ni::Int=1, nf::Int=9,
                  layout::Tuple{Int,Int}=(3,3),
                  figsize=(1500, 500)) where {T<:Real}

    @assert size(DX) == size(FX) "DX and FX must have the same dimensions"

    # Extract nonzero indices/values (works for sparse or dense)
    rowsD, colsD, dxs = findnz(DX)
    rowsF, colsF, fxs = findnz(FX)
    @assert rowsD == rowsF && colsD == colsF "DX and FX must have the same nonzero positions"

    total = length(dxs)
    @assert 1 ≤ ni ≤ nf ≤ total "Need 1 ≤ ni ≤ nf ≤ number of traces ($total)"

    nx, ny = layout
    needed = nf - ni + 1
    @assert needed ≤ nx * ny "Layout $(layout) too small for $needed plots"

    # Build each subplot
    plots = Vector{Any}(undef, nx * ny)
    idx = 1
    for k in ni:nf
        i, j = rowsD[k], colsD[k]
        raw = dxs[k]
        fit = fxs[k]

        p = plot(1:length(raw), raw;
                 color=:gray, lw=1,
                 label="Raw @($i,$j)",
                 xlabel="Time step", ylabel="Intensity",
                 title="Trace ($i,$j)",
                 legend=:topright, grid=true)

        plot!(p, 1:length(fit), fit;
              color=:red, lw=2, label="Fit")

        plots[idx] = p
        idx += 1
    end

    # Fill remaining slots with empty plots
    for k in idx:(nx*ny)
        plots[k] = plot()  # blank
    end

    # Combine into a single figure
    return plot(plots...; layout=layout, size=figsize)
end


"""
    plot_trace(df::DataFrame, mol_id::Int)

Plot the trace of a given molecule from a dataframe.

# Arguments
- `df`: DataFrame with columns `nmol`, `nstep`, `stepHeight`, `stepTime`, `stepLength`.
- `mol_id`: Integer representing the renumbered molecule index to plot.

# Description
Assumes that `nmol` values have been renumbered to start from 1. Plots a step function with heights from `stepHeight`,
transitions at `stepTime`, and final segment length taken from `stepLength`.
"""
function plot_trace(df::DataFrame, mol_id::Int)
    # Filter rows for the specified molecule
    mol_df = filter(:nmol => ==(mol_id), df)
    sort!(mol_df, :stepTime)

    xs = Float64[]
    ys = Float64[]

    for i in 1:nrow(mol_df)
        h = mol_df.stepHeight[i]
        t_start = mol_df.stepTime[i]
        if i < nrow(mol_df)
            t_end = mol_df.stepTime[i+1]
        else
            t_end = t_start + mol_df.stepLength[i]
        end
        append!(xs, [t_start, t_end])
        append!(ys, [h, h])
    end

	if length(ys) >0
		ymax = maximum(ys) + 10
		xmax = maximum(xs) 
    	#plot(xs, ys, ylims=(0, ymax), xlims=(0, xmax), label="Trace $mol_id", lw=2, 
		# xlabel="Time", ylabel="Step Height", title="")
		plot(xs, ys, xlims=(0, xmax), label="Trace $mol_id", lw=2, 
		 xlabel="Time", ylabel="Step Height", title="")
	else
		plot(label="Trace $mol_id", xlabel="Time", ylabel="Step Height", title="")
	end
end


function plot_traces(df::DataFrame, mol_ids::Vector{Int}; lyt=(10,5), size=(1200, 800))
	P = []
	
	for mol in mol_ids
		push!(P, plot_trace(df, mol))
	end
	plot(P..., layout=lyt, size=size)
end


function plot_sdf(sdf::DataFrame; nbins=(10,50,50,50), 
	              xlim=((0.0,10.0),(-100.0, 800.0), (0.0, 200.0), (0.0, 200.0)), 
				  log=(false,true,true,true))
	pnst = step_hist(sdf.nstep .-1;
	                   nbins= nbins[1],
	                   xlim = xlim[1],
					   logy=log[1],
	                   xlabel = "number of steps",
	                   ylabel= "Frequency")
	
	
	psth = step_hist(sdf.stepHeight;
	                   nbins= nbins[2],
	                   xlim = xlim[2],
					   logy=log[2],
	                   xlabel = "step heigth",
	                   ylabel= "Frequency")
	pstt = step_hist(sdf.stepTime;
	                   nbins= nbins[3],
	                   xlim = xlim[3],
					   logy=log[3],
	                   xlabel = "step time",
	                   ylabel= "Frequency")
	pstl = step_hist(sdf.stepLength;
	                   nbins= nbins[4],
	                   xlim = xlim[4],
					   logy=log[4],
	                   xlabel = "step length",
	                   ylabel= "Frequency")

	plot(pnst, psth, pstt, pstl, layout=(2, 2), size=(1000, 800))
	
end



"""
    plot_traces_by_nmolx(df::DataFrame, molx_id::Int; lyt=(3, 3), size=(1200, 800))

Plot traces for all molecules with a given `nmolx` group in the provided DataFrame.

# Arguments
- `df`: DataFrame containing columns `nmolx` and `nmol`.
- `molx_id`: The `nmolx` group to select.
- `lyt`: Tuple specifying the subplot layout.
- `size`: Tuple specifying the plot size.
"""
function plot_traces_by_nmolx(df::DataFrame, molx_id::Int; lyt=(3, 3), size=(1200, 800))
    # Filter rows for molecules with the same nmolx
    mol_df = filter(:nmolx => ==(molx_id), df)
    # Get unique nmol values from this subset
    nmol_ids = unique(mol_df.nmol)
    # Plot the corresponding traces
    plot_traces(df, nmol_ids; lyt=lyt, size=size)
end


###
# Data Frame


function renumber_nmol!(df::DataFrame)
    unique_ids = unique(df.nmol)
    id_map = Dict(id => i for (i, id) in enumerate(sort(unique_ids)))
    df.nmol .= [id_map[id] for id in df.nmol]
    return df
end



end #module
