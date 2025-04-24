using DataFrames, Distributions, Random
using ImageFiltering
using .SimpleLogger
"""Compute the radiative rate"""
function rad_rate(fm::FBMolecule, exc::LaserExcitation, row)
    B = cross_section(fm, exc.λ)
    I = intensity(exc, row[:x], row[:y])
    B * sin(row[:theta])^2 * sin(row[:phi])^2 * I
end

"""Compute the MO efficiency"""
function mo_eff(NA, row)
    delta = asin(NA)
    tetha = row[:theta]
    sincosm = sin(tetha - delta) * cos(tetha - delta)
    sincosp = sin(tetha + delta) * cos(tetha + delta)
    #sincos0 = sin(delta) * cos(delta)  # computed but not used further
    return (2 * delta - sincosp + sincosm) / π / 2
end

""" Compute a time trajector based on exponential waiting times"""
function trajectory(molecule::FBMolecule, r::typeof(1.0Hz))
    t = 0.0s
    hist = []
    dpb = r / molecule.M             # rate/PB-rate 

    debug("typeof(molecule.M) = $(typeof(molecule.M))")
    ddark = r / molecule.Darks[1]    # rate/PB-rate

    debug("typeof(ddark) = $(typeof(ddark))")
    
    argtdark = (1.0/ddark*Hz)

    debug("typeof(argtdark) = $(typeof(argtdark))")

    argtpb = (1.0/dpb*Hz)
    argtgr = (molecule.Darks[2])/1.0s

    debug("dpb = $(dpb), ddark = $(ddark)")
    tdark = rand(Exponential(argtdark)) * 1.0s
    tpb = rand(Exponential(argtpb)) * 1.0s

    debug("tdark = $(tdark), tpb = $(tpb)")
    debug("typof (dark) = $(typeof(tdark)), tpb = $(tpb)")
    while tdark < tpb
        t += tdark
        push!(hist, t)
        debug("time in dark = $(t)")
       
        tground = rand(Exponential(argtgr)) * 1.0s
        debug("tground = $(tground) ")
        t += tground
        push!(hist, t)
        debug("time ground = $(t) ")
        tdark = rand(Exponential(argtdark)) * 1.0s
        tpb = rand(Exponential(argtpb)) * 1.0s

        debug("new dice: tdark = $(tdark), tpb = $(tpb)")
    end
    t += tpb
    push!(hist, t)
    debug("time to photobleach = $(t)")
    hist
end

function spectral_fraction(fb::FBMolecule, λmin::typeof(1.0nm), λmax::typeof(1.0nm))
	quadgk(fb.fm.pdf, λmin/nm, λmax/nm)[1]
end


"""Generate data based on the molecule, excitation, and a DataFrame sample"""
function generate_data(molecule::FBMolecule, exc::LaserExcitation,
					   sample::DataFrame, NA::Float64, 
					   λmin::typeof(1.0nm), λmax::typeof(1.0nm))
    data = deepcopy(sample)
    rs = []
    moeffs = Float64[]
	spx  = Float64[]
    trajs = []
    names = String[]
	#molecule = fmb.fm
    
    for row in eachrow(sample)
		rate = rad_rate(molecule, exc, row)
        push!(rs, rate)
		push!(spx, spectral_fraction(molecule, λmin, λmax))
        push!(moeffs, mo_eff(NA, row))
        push!(trajs, trajectory(molecule, rate))
        push!(names, molecule.name)
    end
    
    data.R = rs
    data.MOeff = moeffs
    data.Trajectories = trajs
    data.Name = names
	data.Spx = spx
    
    return data
end


"""real (MC generated, no semaring) traces"""
function real_trace(row::DataFrameRow, t_meas::Unitful.Time, t_resol::Unitful.Time)
    
    # the time interval goes from 0 to the t_mean (-t_resol) in steps of t_resol
    tmax = (t_meas - t_resol)/s
    tstep = t_resol/s
    t = collect(0:tstep:tmax)

    debug("Compute real trace from 0 to $(tmax) s in steps of $(tstep) s")
    
    # Initialize the trace vector to ones (same length as t)
    trace = ones(length(t))
    # load the trajectories
    times = row[:Trajectories] ./s
    
    # Calculate the scaling intensity factor  ( R x Moesff x SPX)
    # MOeff is the molecule efficiency, SPX the spectral factor, both computes in the 
    # input DF
    int0 = row[:R]/Hz * row[:MOeff] * row[:Spx]
    
    #   for i=1: (-1)^1 = -1, i=2: (-1)^2 = 1, etc.
    for (i, ti) in enumerate(times)
        factor = (-1)^i
        # Create a heaviside step function: heaviside(x) = 1 if x>=0, else 0.
        # We compute (t - ti) elementwise, and then use ifelse broadcasting.
        heaviside = ifelse.((t .- ti) .>= 0, 1.0, 0.0)
        trace .+= factor .* heaviside
    end
    
    # Multiply the final trace by int0
    trace .*= int0
    # return with proper units
    return t*1.0s, trace*1.0Hz
end


"""average the trace based on exposure time (t_exp)"""
function measured_trace(t::Vector{typeof(1.0s)}, arr::Vector{typeof(1.0Hz)}, texp::Unitful.Time)
    tres = (t[2] - t[1])                   # Calculate time resolution
	
    binning = Int(round(texp / tres))     # Compute number of points per bin)

    # Determine the number of elements that fit evenly into bins
    n = (length(arr) ÷ binning) * binning

	#println("mesured trace: binning =$(binning), n=$(n)")

    # Truncate arrays to fit exactly into bins
    arr = arr[1:n]
    t = t[1:n]

    # Reshape arrays into matrices with 'binning' rows and average each column
    # then, multiply by time of expostion the array.

    arr_b = mean(reshape(arr, binning, :), dims=1)[:] .* texp
    t_b = mean(reshape(t, binning, :), dims=1)[:]

	#println("size of arrb = $(length(arr_b))")
	#println("size of tb = $(length(t_b))")

    return t_b, arr_b
end


"""compute traces for each row of a DataFrame"""
function traces(data::DataFrame, rinfo::Dict{Any, Any})
    
    traces = []
    
    # Iterate over each row of the DataFrame
    for row in eachrow(data)
        # Compute the real trace 
        t, trace = real_trace(row, rinfo["t_meas"], rinfo["t_resol"])
		
        # Compute the measured trace from the real trace 
        t, mtrace = measured_trace(t, trace, rinfo["t_exp"])
		

        # Append the measured trace to the array of traces
        push!(traces, mtrace)
    end

    return traces
end


"""create a DataFrame with binned trace coordinates and corresponding traces"""
function df_traces(data::DataFrame, rinfo::Dict{Any, Any})
    
    # Obtain traces 
    trs = traces(data, rinfo)
    
    # Compute spatial resolution based on pixel scale and division factor
    res = rinfo["pixelscale"] / rinfo["pixeldiv"]
    info("Generate traces with a resolution $(res)")
    
    # Calculate binned indices i and j from positions x and y (convert microns to indices)
    # Julia uses integer division operator `÷` (equivalent to Python's `//`)
    i_indices = Int.(data.x  .÷ rinfo["res"])
    j_indices = Int.(data.y  .÷ rinfo["res"])

    # Construct DataFrame with indices and corresponding traces
    TracesDf = DataFrame(i = i_indices, j = j_indices, traces = trs)

    # Filter out rows where indices exceed the provided limits (edges lost due to binning)
    TracesDf = filter(row -> row.i < rinfo["NHR"][1] && row.j < rinfo["NHR"][2], TracesDf)
    
    return TracesDf
end


""" 
Prepare a high resolution image in the CMOS

	- returns an image im[hrdim x hrdim] such that each pixel im[i,j] contains a trace.

	- NB, HR image comes in counts
"""
function hr_image(n::Int64, tracesDf::DataFrame, rinfo::Dict{Any, Any})
	# Initialize HR image array with zeros
	im = zeros(Float64, rinfo["NHR"][1], rinfo["NHR"][2])  

	# Fill in the image from the trace data for molecule number "i"
	
	for row in eachrow(tracesDf)
		im[row.i + 1, row.j + 1] = row.traces[n + 1]
	end

	im
end


"""
- Compute the low resolution obtained from the HR image 
- Pass a gaussian filter to reflect diffractive limit

"""
function frame2D(n::Int64, tracesDf::DataFrame, rinfo::Dict{Any, Any})
	function rebin(M, binning)
		NHR = size(M)
		rows = (NHR[1] ÷ binning) * binning
		cols = (NHR[2] ÷ binning) * binning

		#trim the matrix, to fit potential mismatches due to integer division
		Mt = M[1:rows, 1:cols]
		# Reshape into (rows_per_bin, num_bins_rows, cols_per_bin, num_bins_cols)
		M_rebinned = reshape(Mt, binning, 
							 size(Mt,1) ÷ binning, 
							 binning, 
							 size(Mt,2) ÷ binning)
		
		# Group the bins by their row and column position
		#We reorder the dimensions to group the data by bin position:
	    # 1st dim (2) → row bins
	    # 2nd dim (4) → column bins
	    # 3rd dim (1) → rows within the bin
	   # 4th dim  (3) → cols within the bin
		# so now we have (num_bins_rows, num_bins_cols, bin_h, bin_w)
		M_rebinned = permutedims(M_rebinned, (2, 4, 1, 3))
		
		# Sum over each bin×bin block (dims 3 and 4)
		M_binned = sum(M_rebinned, dims=(3, 4))
		
		# Remove singleton dimensions to get a rows×col result
		M_binned = dropdims(M_binned, dims=(3, 4))
	end
	
    # Get HR image
    im = hr_image(n, tracesDf, rinfo)
    
    dlr = rinfo["dl"] /2.0 # diffraction limit radius
    res = rinfo["res"]     # resolution in the HR image = pixelscale/pixeldiv

    # Pass a gaussian filter to simulate the effect of difractive limit. 
    sigma = uconvert(nm, dlr)/uconvert(nm, res)

    if rinfo["gaussf"]==true
        imgf = imfilter(im, reflect(Kernel.gaussian(sigma)))
    else
        imgf = im
    end 

    # Reshape and sum to bin the image
    binning = rinfo["pixeldiv"]
    debug(" BINNING = $(binning)")
	im_binned = rebin(imgf, binning)
	
    debug("size of im_binned = $(size(im_binned))")
    
    debug(" constant background = $(rinfo["bck"])")
    debug("t_exp = $(rinfo["t_exp"])")

    im_binned .+ rinfo["bck"] .* rinfo["t_exp"]
end



"""Add noise shot noise, electronics noise and dark current noise to the frame"""
function frame2Dn(f2D::Matrix{Float64}, rinfo::Dict{Any, Any}, cam::CMOS;  λ=550nm)

    ns = noise(cam, λ, f2D, rinfo["t_exp"])
    debug("adding noise: size of noise matrix =$(size(ns))")
    im_n = [rand(Normal(f2D[i, j], 
                        ns[i, j])) for i in axes(f2D, 1), j in axes(f2D, 2)]
    im_n
end


"""Generates the temporal evolution of a simulated image sequence over time, 
using trace data (tracesDf).

Each time point corresponds to a frame/image where intensities are computed using the CPS function.

This is done like this:

- im0 = CPS(0, TracesDf, rinfo) 

Calls the function CPS to compute the image at time index 0 (first time frame.)

- data = [zeros(size(im0)), im0]

Initializes a list of 2D arrays (images).

	•	The first element is a black (zeroed) image of the same size as im0.

	•	The second is im0 itself.

	•	This is a workaround to build an array of consistent shapes — later the black frame will be removed.

- N = length(TracesDf[1, :Traces])

	• This gets the number of time points by checking the length of the trace vector in the first molecule.

- for i in 2:N-1
    push!(data, CPS(i, TracesDf, rinfo))
end

	•	Iterates over time indices starting from 2 to N-1.

	•	For each time index i, it computes the corresponding image using CPS(...) and appends it to data.

- reduce(vcat, [reshape(d, 1, size(d,1), size(d,2)) for d in data])

Each 2D image d is reshaped to a 3D array with shape (1, height, width).

	•	This adds a time dimension in the first axis.

	•	vcat (vertical concatenation) stacks all these time frames along the first dimension, yielding a final 3D array:

	•	Shape: (T, H, W) — T time steps, image of height H and width W.
"""
function frame3D(tracesDf::DataFrame, rinfo::Dict{Any, Any})

    # generate first frame and prepre 3D matrix
	im0 = frame2D(0, tracesDf, rinfo)
	#data = [zeros(size(im0)), im0]
    data = [im0]
	
	N = length(tracesDf[1, :traces])

    # fill all other frames
	for i in 1:N-1
		push!(data, frame2D(i, tracesDf, rinfo))
	end
	return reduce(vcat, [reshape(d, 1, size(d,1), size(d,2)) for d in data])
end


"""Add noise to frame3D"""
function frame3Dn(tracesDf::DataFrame, rinfo::Dict{Any, Any}, cam::CMOS)
	im0 = frame2Dn(frame2D(0, tracesDf, rinfo), rinfo, cam)
	data = [im0]
	
	N = length(tracesDf[1, :traces])
	for i in 2:N-1
		d2 = frame2Dn(frame2D(i, tracesDf, rinfo), rinfo, cam)
		push!(data, d2)
	end
	return reduce(vcat, [reshape(d, 1, size(d,1), size(d,2)) for d in data])
end