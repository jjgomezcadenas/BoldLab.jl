module JPELT
using Statistics
using StatsBase
using SparseArrays
using DataFrames
using LinearAlgebra
using Changepoints
using Distributions


export pelt_fit, getsteps_pelt, steps_from_pelt, pelt_find_optimal_noise_region


include("SimpleLogger.jl")
using .SimpleLogger
global_log_level = WARN



"""
    pelt_fit(tz::AbstractVector{T}; i::Int=1, f::Int=5) where {T<:Real}

Apply the CROPS (Changepoints for a Range of PenaltieS) version of the PELT algorithm to a single trace `tz`, using a Gaussian cost model with known standard deviation and unknown mean.

# Arguments
- `tz`: A 1D real-valued vector representing the time trace to analyze.
- `i`: Starting index of the initial segment used to estimate the noise standard deviation `σ` (default: 1).
- `f`: Final index of the segment used to estimate `σ` (default: 5). The range `[i:f]` should represent a flat region of the trace.

# Returns
- A dictionary `crops_output` containing the output of the `CROPS` routine, including:


  * `out["penalty"]` : vector of penalties between minimum and maximum pen
  * `out["number"]` : vector of number of changepoints for each penalty
  * `out["constrained"]` : vector of optimal segmentation costs for each penalty
  * `out["changepoints"]` : vector of changepoints sets for each penalty
  

# Notes
- The cost function used is `NormalMeanSegment`, which assumes normally distributed data with known variance and unknown mean.
- The penalty range is set between `log(n)` (likely to overfit) and `100 × log(n)` (encouraging fewer changepoints).
- The optimal number of changepoints for each penalty is selected internally by the CROPS algorithm.

"""

function pelt_fit(tz::Vector{T}; i::Int=1, f::Int=5, auto_noise_region::Bool=true) where {T<:Real}

	 # Automatically find optimal noise region if requested
    if auto_noise_region
        i, f = pelt_find_optimal_noise_region(tz)
        debug("Auto-selected noise region: [$i:$f]")
    end

	# Check if trace is long enough for the requested noise estimation range
    if length(tz) < f
        debug("Trace too short ($(length(tz))) for noise estimation range si=$i, sf=$f")
        return nothing, i, f
    end

	# Estimate the std of the distribution for range i:f
	σi = std(tz[i:f])

	# compute the penalty functions
	pn1 =log(length(tz))  # minimum penalty (tends to overfit)
	pn2 = 100 * pn1 # defined the maximum range (PELT will mange optimal)
	
	# Cost distribution is normal with known std and unknown mean
	costfi = NormalMeanSegment(tz, σi)

	# run CROPS
	crops_output = CROPS(costfi , length(tz), pn1, pn2)
	#crops_output = @PELT tz Normal(:?, σi) pn1 pn2

	return crops_output, i, f
end

"""
    steps_from_pelt(tz, cro, ic)

Generate a piecewise constant fit from a PELT changepoint analysis result.

# Arguments
- `tz`: A real-valued vector representing the time trace to be segmented.
- `cro`: A dictionary output from the `CROPS` routine, containing keys like `"changepoints"`, `"penalty"`, and `"constrained"`.
- `ic`: Index of the changepoint set to use from `cro["changepoints"]` (default: 1).

# Returns
A tuple of:
1. A vector representing the stepwise fit to the input trace (`tz`), where each segment has a constant value equal to the segment mean.
2. The penalty value associated with the selected changepoint set.
3. The cost value associated with the selected changepoint set.

# Notes
- If the changepoint list has fewer than 2 elements, the function returns an empty fit and zeros for penalty and cost.
- If `ic` exceeds the number of available changepoint sets, it defaults to 1 with a warning.
- The returned fit has the same length as `tz` and is constructed by repeating the mean of each segment over its duration.
- The final segment is extended explicitly from the last changepoint to the end of the trace.
"""
function steps_from_pelt(tz::Vector{T}, cro, ic::Int) where {T<:Real}

	debug("Calling steps_from_pelt, ic =$(ic)")
	debug("cro, ic =$(cro)")
	
	if length(cro["changepoints"][ic] ) <2
		return [], 0.0, 0.0, 0
	end

	if length(cro["changepoints"] )< ic
		warn("ic =$(ic) is larger than fit set, take ic = 1")
		ic = 1
	end
	
	c1 = cro["changepoints"][ic]
	pen = cro["penalty"][ic]
	cost = cro["constrained"][ic]
	nc = cro["number"][ic]

	values = Float64[]
	lengths = Int64[]
	
	# Handle changepoints properly - they mark boundaries between segments
	# If first changepoint is 0, start from index 1
	start_idx = c1[1] == 0 ? 1 : c1[1]
	
	for i in 2:length(c1)
		# Each segment goes from previous changepoint to current changepoint
		t0 = start_idx
		t1 = c1[i]
		if t1 > t0  # Only add segment if it has positive length
			push!(values, mean(tz[t0:t1]))
			push!(lengths, t1 - t0 + 1)
		end
		start_idx = t1 + 1  # Next segment starts after current changepoint
	end
	
	# Add the final segment from the last changepoint to the end
	if start_idx <= length(tz)
		push!(values, mean(tz[start_idx:end]))
		push!(lengths, length(tz) - start_idx + 1)
	end

	vcat([fill(val, len) for (val, len) in zip(values, lengths)]...), pen, cost, nc
			  	
end


"""
    getsteps_pelt(tz::AbstractVector{<:Real}, c1)

Extract step statistics from a time trace `tz` given a list of changepoint positions `c1`.

# Arguments
- `tz`: A real-valued vector representing the time trace.
- `c1`: A vector of changepoint indices, including the starting (0 or 1) and ending positions.

# Returns
A tuple of:
1. `nsteps`: Number of steps (equal to `length(c1) - 1`)
2. `hh`: Vector of mean values for each segment (step height)
3. `ht`: Vector of changepoint indices (step times)
4. `hl`: Vector of step durations (step lengths)

# Notes
- Step durations are computed as `t1 - t0 + 1`, where `t0` and `t1` are consecutive changepoint positions.
- If a changepoint value is 0, it is replaced with 1 to avoid indexing errors.
"""
function getsteps_pelt(tz::Vector{T}, cro, ic::Int) where {T<:Real}
	c1 = cro["changepoints"][ic]
    nsteps = length(c1) -1
	hh = Float64[]
	ht = Int[]
	hl = Int[]

	for i in 2:length(c1) 
		t0 = ifelse(c1[i-1] == 0, 1, c1[i-1])
		t1 = c1[i]
		push!(ht, t1)
		push!(hl, t1-t0+1)
		push!(hh, mean(tz[t0:t1]))
	end
	
	nsteps, hh, ht, hl 
end


"""
    pelt_find_optimal_noise_region(tz::Vector{T}; n_chunks=5, min_chunk_size=10) -> (Int, Int)

Find the optimal region for noise estimation by analyzing trace in chunks.
Returns the start and end indices of the chunk with the lowest coefficient of variation.

# Arguments
- `tz`: Trace vector
- `n_chunks`: Number of chunks to divide the trace into (default: 5)  
- `min_chunk_size`: Minimum size for each chunk (default: 10)

# Returns
- Tuple (si, sf) with optimal start and end indices for noise estimation
"""
function pelt_find_optimal_noise_region(tz::Vector{T}; n_chunks=10, min_chunk_size=10, low_mean=1e-5) where {T<:Real}

    n = length(tz)
    
    # Ensure we have enough data for chunking
    if n < n_chunks * min_chunk_size
        # Fall back to using last 20% of trace or minimum chunk size
        chunk_size = max(min_chunk_size, div(n, 5))
        return max(1, n - chunk_size + 1), n
    end
    
    chunk_size = div(n, n_chunks)
    best_cv = Inf
    best_si, best_sf = 1, min(chunk_size, n)
    
    for i in 1:n_chunks
        si = (i - 1) * chunk_size + 1
        sf = min(i * chunk_size, n)
        
        if sf - si + 1 >= min_chunk_size
            chunk = tz[si:sf]
            chunk_mean = mean(chunk)
            chunk_std = std(chunk)
            
            # Use coefficient of variation (CV = std/mean) as noise metric
            # For low/near-zero signals, use std directly
            cv = chunk_mean > low_mean ? chunk_std / abs(chunk_mean) : chunk_std
            
            if cv < best_cv
                best_cv = cv
                best_si, best_sf = si, sf
            end
        end
    end
    
    debug("Optimal noise region: [$best_si:$best_sf], CV = $best_cv")
    return best_si, best_sf
end





end #end module 