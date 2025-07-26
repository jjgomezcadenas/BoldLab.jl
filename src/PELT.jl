module JPELT
using Revise
using StatsPlots
using Statistics
using StatsBase
import Measures
using SparseArrays
using DataFrames
using LinearAlgebra
using Changepoints
using Distributions


export stepfindcore, auto_step_main

include("SimpleLogger.jl")
using .SimpleLogger
global_log_level = WARN



"""
    find_fit_pelt(trzs, df; ic=1, si=1, sf=5, ped=0.0)

Apply changepoint detection using the PELT algorithm to a set of time traces extracted from a trace matrix `trzs`, using the coordinate list from the dataframe `df`.

# Arguments
- `trzs`: A 2D matrix where each element is a vector representing a time trace (`trzs[i, j]` is the trace for pixel `(i, j)`).
- `df`: A DataFrame with columns `i` and `j` indexing the pixels of interest in `trzs`.
- `ic`: Index of the cost function used in the `changepoints` dictionary (default: 1).
- `si`: Start index of the trace to fit (default: 1).
- `sf`: End index of the trace to fit (default: 5).
- `ped`: Pedestal value to subtract from each trace before storing the results (default: 0.0).

# Returns
A tuple of six items:
1. `df2`: A DataFrame containing one row per detected step, with columns:
   - `i`, `j`: pixel indices
   - `nmol`: molecule counter
   - `bestShot`: value of the cost for the fit
   - `nstep`: number of steps found
   - `stepHeight`: height of each individual step
   - `stepHeightMin`: minimum step height among all steps
   - `stepTime`: time of the step (index)
   - `stepLength`: duration of each step
2. `I`: Vector of `i` indices for successful fits.
3. `J`: Vector of `j` indices for successful fits.
4. `DX`: Vector of raw traces after pedestal subtraction.
5. `FX`: Vector of fitted traces after pedestal subtraction.
6. `SC`: Vector of changepoint detection result dictionaries (from `pelt_fit`).

# Notes
- The function filters out fits where either the penalty or cost returned is not strictly positive.
- The `debug(...)` statements help trace progress and data quality.
"""


function find_fit_pelt(trzs::AbstractVector{<:AbstractVector{T}}, df::DataFrame; 
                       ic =1, si=1, sf=5, ped=0.0) where {T<:Real}
	I = Int[]
	J = Int[]
    DX = Vector{Vector{T}}()
    FX = Vector{Vector{T}}()
	SC = Vector{Dict{String, Array}}()
	
	ng = 0

	df2 = DataFrame(i=Int[], j=Int[], nmol=Int[], bestShot=Float32[], 
					nstep=Int[], stepHeight=Float32[], stepHeightMin=Float32[], 
					stepTime=Int[], stepLength=Int[])

	for row in eachrow(df)
		debug("fitting trace ($(row.i),$(row.i))")
		
		tz = trzs[row.i, row.j] 

		debug("mean between $(si) and $(sf) $(mean(tz[si:sf]))")
		
		cro = pelt_fit(tz; i=si, f=sf)
		c1 = cro["changepoints"][ic]
		pfit, pen, cost =  steps_from_pelt(tz, cro; ic=ic)
		
		if pen >0 && cost > 0 # fit OK
			ng+=1
			push!(I,row.i)
			push!(J,row.j)
			push!(DX, tz .-ped)
			push!(FX, pfit .- ped)
			push!(SC, cro)

			nsteps, sth, stt, stl = getsteps_pelt(tz,c1)  
			sthmx = minimum(sth)
			for k in 1:nsteps
                push!(df2, (row.i, row.j, ng, cost, nsteps, 
							sth[k] - ped, sthmx - ped, stt[k], stl[k]))
			end
            
		end
	end
	
	df2, I, J, DX, FX, SC
end


"""
    pelt_fit(tz::AbstractVector{T}; i::Int=1, f::Int=5) where {T<:Real}

Apply the CROPS (Changepoints for a Range of PenaltieS) version of the PELT algorithm to a single trace `tz`, using a Gaussian cost model with known standard deviation and unknown mean.

# Arguments
- `tz`: A 1D real-valued vector representing the time trace to analyze.
- `i`: Starting index of the initial segment used to estimate the noise standard deviation `σ` (default: 1).
- `f`: Final index of the segment used to estimate `σ` (default: 5). The range `[i:f]` should represent a flat region of the trace.

# Returns
- A dictionary `crops_output` containing the output of the `CROPS` routine, including:
  - `"penalties"`: a vector of penalty values explored
  - `"costs"`: corresponding cost values
  - `"changepoints"`: a vector of changepoint sets, one for each penalty

# Notes
- The cost function used is `NormalMeanSegment`, which assumes normally distributed data with known variance and unknown mean.
- The penalty range is set between `log(n)` (likely to overfit) and `100 × log(n)` (encouraging fewer changepoints).
- The optimal number of changepoints for each penalty is selected internally by the CROPS algorithm.
"""
function pelt_fit(tz::AbstractVector{T}; i::Int=1, f::Int=5) where {T<:Real}
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
function getsteps_pelt(tz::AbstractVector{T}, c1) where {T<:Real}
    nsteps = length(c1) -1
	hh = []
	ht =[]
	hl = []

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
    steps_from_pelt(tz, cro; ic=1)

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
function steps_from_pelt(tz, cro; ic=1)

	debug("Calling steps_from_pelt, ic =$(ic)")
	debug("cro, ic =$(cro)")
	
	if length(cro["changepoints"][1] ) <2
		return [], 0.0, 0.0
	end

	if length(cro["changepoints"] )< ic
		warn("ic =$(ic) is larger than fit set, take ic = 1")
		ic = 1
	end
	
	c1 = cro["changepoints"][ic]
	pen = cro["penalty"][ic]
	cost = cro["constrained"][ic]

	values = Float64[]
	lengths = Int64[]
	for i in 1:length(c1) -1 
		t0 = ifelse(c1[i] == 0, 1, c1[i])
		t1 = c1[i+1]
		push!(values, mean(tz[t0:t1]))
		push!(lengths, t1-t0+1)
	end
	
	push!(values, mean(tz[length(c1)+1:length(tz)]))
	push!(lengths, 1)
	
	c1e = c1[end]
	tze = length(tz)

	debug("c1 = $(c1)")
	debug("mean from i = $(c1e) to j = $(tze)")
	debug("mean f= $(mean(tz[c1e:tze]))")
	
	push!(values, mean(tz[c1e:tze]))
	push!(lengths, length(tz) -length(c1)+1)

	vcat([fill(val, len) for (val, len) in zip(values, lengths)]...), pen, cost
			  	
		
end

end #end module 