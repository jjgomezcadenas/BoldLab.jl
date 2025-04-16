### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 34233420-116a-11f0-0e9c-75681c7ddb3f
using Pkg; Pkg.activate(ENV["JBoldLab"])

# ╔═╡ 3738ca4c-c6e2-4a4d-90f2-1ee3284b7c89
push!(LOAD_PATH, ENV["JBoldLab"] * "/src")

# ╔═╡ bba8e102-398c-4191-9d88-5cbc7879d6cf
begin
	using Revise
	using BoldLab
	using SimpleLogger
	#using dffunctions
end

# ╔═╡ bf3dd6ae-7cf8-45b4-b407-f480cb9e7f30
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Plots 
	#using Measures
	using Printf
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using DelimitedFiles
	using Images, FileIO
	using Interpolations
	using QuadGK
	using Distributions, Random
end

# ╔═╡ 1b562d29-0419-4c3b-b699-271888b547a5
begin
	using Unitful
	using UnitfulEquivalences
	using PhysicalConstants.CODATA2018
	
end

# ╔═╡ e6f2b315-b5b0-4983-a1c0-4b396b337b97
using ImageFiltering

# ╔═╡ 5b04c22a-1a29-4d17-be3b-3a42fa1ecdf9
import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W

# ╔═╡ c7301bb6-e147-4998-883e-2fa1c542c14e
barn = 1E-24 * cm^2

# ╔═╡ 934b4ebc-a4f4-4900-9d19-f8b04726a635
ENV["JBoldLab"]

# ╔═╡ 21345bd9-5cd7-4775-805c-e1407072a616
names(BoldLab)

# ╔═╡ a650449d-6f90-4655-ae83-64fe4eb67b22
names(SimpleLogger)

# ╔═╡ 41053e81-3b76-41f3-9f54-076518a64e98
#names(dffunctions)

# ╔═╡ 4badd96d-2461-4507-b771-1b4882cd239d
function plot_molecule(n3f, n3c)
	p1 = plot(n3f.em_spectrum[1], n3f.em_spectrum[2],
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "a.u.",
			  title = "NAPH3 measured emission spectrum",
		      legend = :topright)
	p1 = plot!(p1, n3c.em_spectrum[1], n3c.em_spectrum[2],
			  label = "NAPH3-Ba")
	
	ls = 400:10:800
	p2 = plot(ls, n3f.pdf.(ls),
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "PDF",
			  title = "NAPH3 emission PDF",
		      legend = :topright)
			  
	p2 = plot!(p2, ls, n3c.pdf.(ls),
			  label = "NAPH3-Ba")

	p3 = plot(n3f.abs_spectrum[1], n3f.abs_spectrum[2],
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "\$\\epsilon (M^{-1} cm^{-1})\$",
			  title = "NAPH3 \$\\epsilon \$",
		      legend = :topright)
	p3 = plot!(p3, n3c.abs_spectrum[1], n3c.abs_spectrum[2],
			  label = "NAPH3-Ba")
	
	ls = 200:1:500
	p4 = plot(ls, n3f.xs.(ls),
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "\$\\sigma\$ (barn)",
			  title = "NAPH3 \$\\sigma\$ (barn)",
		      legend = :topright)
	p4 = plot!(p4, ls, n3c.xs.(ls),
			  label = "NAPH3-Ba")
	plot(p1,p2, p3, p4, layout=(2,2), size=(700, 500))
end

# ╔═╡ 8d8cc6ff-3006-4cd1-9308-15b4c449a998
md"""
## NAPH3 molecules
"""

# ╔═╡ 1d59335e-7c76-49c2-88f6-aef6e4312941
begin
	jbl = ENV["JBoldLab"]
	naph3f, naph3c = fm_from_csv(joinpath(jbl, "data"), 
									  "naph3_emission.csv",
									  "naph3_abs.csv")
end

# ╔═╡ a2076949-49f3-4964-a2f8-6b3d6c6303fb
plot_molecule(naph3f, naph3c)

# ╔═╡ f163b606-952f-4d14-9e9c-e9d9ce8722fd
md"""
## Define Setup
"""

# ╔═╡ cbfd9e9c-2faa-4334-8324-dd7351e681b6
lsr = Laser(450nm, 1mW)

# ╔═╡ c3440909-a0c4-4e73-8bc2-953b2d26bbc3
md"""
- photon density for 1 mW
- photon density for 3.115 mW
"""

# ╔═╡ e67eb86d-ea64-4a5c-97c3-15fa67725b7a
I = uconvert(Hz*cm^-2, photon_density(450nm, 1mW, 15*μm, 15*μm))

# ╔═╡ f31cd8db-1d37-4921-9c60-971ef17a3f22
md"""
### Laser: 
- Wavelength = $(lsr.λ)
- Power = $(lsr.P)
- Photon energy = $(photon_energy(lsr.λ))
- Intensity = $(I)
"""

# ╔═╡ e34083c1-526f-46bb-b75b-4717ac600fb7
I2 = uconvert(Hz*cm^-2, photon_density(450nm, 3.115mW, 15*μm, 15*μm))

# ╔═╡ 10be9bea-b310-432a-a14e-f79a0d76b796
obj = Objective("highNA", 0.95, 100.0)

# ╔═╡ 58fe32b8-d516-40f0-a875-5cd4df780355
obj.M

# ╔═╡ 90e48289-e186-4118-a7d7-f911cce4b50e
struct BandFilter
	λmin::Unitful.Length
	λmax::Unitful.Length
end

# ╔═╡ 53768fe9-bfe7-423f-bd5f-8631e11c38c2
bf = BandFilter(450.0nm, 800.0nm)

# ╔═╡ 5ae75425-b152-4099-b2fd-d212a75cca51
md"""
### Objective and filter: 
- NA = $(obj.NA)
- Magnification = $(obj.M))
- Filter λmin = $(bf.λmin)
- Filter λmax = $(bf.λmax)
"""

# ╔═╡ a850b3f4-03f7-4978-b245-54b2460f8825
function oflash4_eff(λ::Unitful.Length)
	l = λ/nm
		if l < 400.0 || l > 900.0
			return 0.
		else
			wl = 400.:50.:900.
			ϵ = [0.4,0.7,0.78,0.82,0.82,0.8,0.72,0.62,0.5,0.37,
			  0.24]
			e = CubicSplineInterpolation(wl, ϵ)
			return e(l)
		end
end

# ╔═╡ 845833df-5a13-489c-b7a3-9250abf92c36
struct CMOS
	name::String
	pixelsize::Unitful.Length  # pixel size is pixelsize x pixelsize
	npixel::Integer            # pixel array is npixel x npixel 
    readout_noise::Float64     # in electrons/pixel
    dark_current::typeof(1.0Hz)      # in electrons/pixel/second
 	binning::Integer           # binning size
	qe::Function               # as a function of λ
	sensorsize::Unitful.Length
	binPixelsize::Unitful.Length 
	binNpixel::Integer

	function CMOS(name, pixelsize, npixel, readout_noise, dark_current, binning, qe)
		sensorsize = pixelsize * npixel
    	binPixelsize = pixelsize * binning
		biNpixel = npixel ÷ binning
		new(name, pixelsize, npixel, readout_noise, dark_current, binning, qe,
		   sensorsize, binPixelsize, biNpixel)
	end
end

# ╔═╡ 8ce9609f-efb4-4760-a4c4-308551669558
begin

	nphe(cam::CMOS, λ, photons) = cam.qe(λ) * photons
	
	# Compute noise (https://camera.hamamatsu.com/eu/en/learn/technical_information/camera_articles/qcmos_vs_emccd.html)

	# Check noise for readout and DC

	function noise(cam::CMOS, λ, photons, texp)
		#println("nphe = $(nphe(cam, λ, photons))")
		#println("noise1 = $((cam.readout_noise * cam.binning)^2)")
		#println("noise2 = $(cam.dark_current * cam.binning^2)")
		#println("texp = $(texp)")
	    return sqrt.(
	        nphe(cam, λ, photons) .+
	        cam.readout_noise * cam.binning .+
	        cam.dark_current * texp/(Hz*s) * cam.binning * cam.binning # area!
	    )
	end
end

# ╔═╡ 872458a6-89ae-4fa2-85b0-220bad7503f4
begin
	of = CMOS("ORCA_Flash4", 6.5μm, 2049, 1.0, 0.06Hz, 4, oflash4_eff)
end

# ╔═╡ 67c79039-c97e-4efb-a1d0-acb6b3667c49
md"""
### CMOS: 
- name = $(of.name)
- pixel size = $(of.pixelsize)
- number of pixels = $(of.npixel)
- readout noise = $(of.readout_noise) p.e.
- dark current = $(of.dark_current)  p.e./pixel/s
- binning = $(of.binning)
- QE @ 550 nm = $(of.qe(550nm))
- sensor size = $(of.sensorsize)
- logical pixel size = $(of.binPixelsize)
- effective number of pixels = $(of.binNpixel)
"""

# ╔═╡ c7d45989-5a77-4b05-805e-1ded0d9beb76
noise(of, 550nm, 100, 0.5s)

# ╔═╡ 9543ff4c-7e48-4009-9a79-812fc3ffffcb
md"""
- Define miscelanous run parameters in dictionary RunInfo
"""

# ╔═╡ 9b37c1b4-8528-4b5e-8a30-ca7bfebda9ad
RunInfo = Dict("N" => 10, # number of molecules to generate
			   "sigma" => 15μm, # sigma of gaussian beam in sample
			   "center" =>(15μm,15μm), # center of laser beam
			   "dimensions" => [30μm, 30μm], # effectise "lasar light spot"
			   "pbcycles"  => 5e+4, # photobleaching cycles (nγ abosorbed before pB)
			   "dkcycles"  => 1e+5, # dark cycles (nγ abosorbed before DT)
			   "tdark"  => 2s, # time in dark states before decaying to ground
			   "t_div" => 5, # number of divisions on time-exp for HR calculation 
			   "t_meas" => 100s, #Measurement time (length of run)
			   "t_exp" => 0.5s, #time of exposition per step
			  )

# ╔═╡ 8ad1d34e-7125-4cc9-934c-11522d6a74df
md"""
- Define the physical molecule
"""

# ╔═╡ a8869a63-7cf4-4d68-86e6-963c0a3e74c9
naph3fB = FBMolecule(naph3f, "naph3-free", RunInfo["pbcycles"], 
					 (RunInfo["dkcycles"], RunInfo["tdark"]))

# ╔═╡ bbca2d6f-a5dc-458e-9488-8ff01e0d92ec
cross_section(naph3fB, lsr.λ)/barn

# ╔═╡ 3ae994ca-44ec-45bf-b23a-9a578b35c5e2
naph3f.xs(lsr.λ/nm)

# ╔═╡ 55d8bfdf-8d6a-4d74-89d6-862fed0d2260
md"""
- Define Laser Excitation 
"""

# ╔═╡ 7d2d9c04-7237-4177-8711-a67df460aec7
begin
	sigma=15μm 
	center=(sigma,sigma)
	lexc = LaserExcitation(lsr, RunInfo["sigma"], RunInfo["center"] )
end

# ╔═╡ 5c24b381-88ce-469e-8311-8136725261c4
md"""
- Generate a sample of $(RunInfo["N"]) molecules with dimensions $(RunInfo["dimensions"])
"""

# ╔═╡ be574290-a0c3-4cb1-990e-f1b8f66e1034
begin
	sample=Sample_3D(10000,RunInfo["dimensions"])
end

# ╔═╡ 18b12b7b-db7d-40e0-be9b-cbf5a6f77438
md"""
- Generate "true data" for this sample, given the laser and the objective
"""

# ╔═╡ c1cd36d1-e43a-476c-ac1c-833b11c576f9
begin
	data=generate_data(naph3fB, lexc, sample, obj.NA, bf.λmin, bf.λmax)
end

# ╔═╡ 523009ad-d221-4f60-a7ed-010fb20c0a09
RunInfo

# ╔═╡ f76cb6da-3b6b-43d2-b020-bd5aff6acab7
function real_trace(row, t_meas, t_resol)
    
    t = collect(0:t_resol/s:(t_meas - t_resol)/s)
    
    # Initialize the trace vector to ones (same length as t)
    trace = ones(length(t))
    
    
    times = row[:Trajectories] ./s
    
    # Calculate the scaling intensity factor 
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
    
    return t*1.0s, trace*1.0Hz
end

# ╔═╡ 8bcbc2c2-fb3e-458d-a643-d55d7e15173c
function plot_real_traces(data, t_meas, t_resol)
	P=[]
	for i in 1:9
		tx, ty = real_trace(data[i,:], t_meas, t_resol)
		#println("tx =$(tx)")
		#println("ty =$(ty)")
		push!(P, plot(tx, ty,
			  label = "NAPH3 steps",
		      xlabel = "t (s)",
		      ylabel = "Rate",
			  title = "NAPH3 steps",
		      legend = false))
	end
	plot(P...,  layout=(3,3), size=(700, 500))
end
	

# ╔═╡ f683cfc7-9bf3-45e7-ab28-fb1532e829dc
begin
	t_div=RunInfo["t_div"] 
	t_meas= RunInfo["t_meas"]
	t_exp = RunInfo["t_exp"]
	t_resol=t_exp/t_div
	RunInfo["t_resol"] = t_resol
	#datax=generate_data(naph3fB, lexc, sample, obj.NA, bf.λmin, bf.λmax)
	plot_real_traces(data, t_meas, t_resol)
end

# ╔═╡ a95a294d-f91f-4ef6-b2df-ee6da0b69155
begin

	dl = ( (bf.λmin + bf.λmax) / 2 ) / (2 * obj.NA) #diffractive limit"
	pixelscale=of.binPixelsize/obj.M 
	pixeldiv=4
	res=pixelscale/pixeldiv
	RunInfo["pixelscale"] = pixelscale
	RunInfo["pixeldiv"] = pixeldiv
	RunInfo["res"] = res
	RunInfo["dl"] = dl
	RunInfo["bck"] = 400Hz # background
	md"""
	- Time of measurement = $(t_meas) 
	- Time of exposition = $(t_exp)
	- Divisions = $(t_div)
	- Time (resolution) = $(t_resol)
	- Diffractive limit = $(dl)
	- pixel scale = $(uconvert(nm, pixelscale))
	- pixel division = $(pixeldiv)
	- resolution =$(uconvert(nm, res))
	"""

end

# ╔═╡ 448f82e4-95d5-4bbb-97c7-9d4269aad25a
begin
	NHR = Int.(RunInfo["dimensions"]   .÷ res) .- 1
	RunInfo["NHR"] = NHR
end

# ╔═╡ 82f10b3c-9cdd-4a0a-81da-610e2160588c
"""average the trace based on exposure time (t_exp)"""
function measured_trace(t, arr, texp)
    res = (t[2] - t[1])                   # Calculate time resolution
	
    binning = Int(round(texp / res))     # Compute number of points per bin)

    # Determine the number of elements that fit evenly into bins
    n = (length(arr) ÷ binning) * binning

	#println("mesured trace: binning =$(binning), n=$(n)")

    # Truncate arrays to fit exactly into bins
    arr = arr[1:n]
    t = t[1:n]

    # Reshape arrays into matrices with 'binning' rows and average each column
    arr_b = mean(reshape(arr, binning, :), dims=1)[:] .* texp
    t_b = mean(reshape(t, binning, :), dims=1)[:]

	#println("size of arrb = $(length(arr_b))")
	#println("size of tb = $(length(t_b))")

    return t_b, arr_b
end

# ╔═╡ 69964f32-021d-4f99-af05-b68d41dc89d0
"""computes traces for each row of a DataFrame"""
function traces(data, info)
    
    traces = []
    
    # Iterate over each row of the DataFrame
    for row in eachrow(data)
        # Compute the real trace 
        t, trace = real_trace(row, info["t_meas"], info["t_resol"])
		#println("real trace")
		#println("t size =$(length(t))")
		#println("trace size: t =$(length(trace))")

        # Compute the measured trace from the real trace 
        t, mtrace = measured_trace(t, trace, info["t_exp"])
		#println("measured trace")
		#println("t size =$(length(t))")
		#println("trace size: t =$(length(mtrace))")

        # Append the measured trace to the array of traces
        push!(traces, mtrace)
    end

    return traces
end

# ╔═╡ 9200aad7-9aa6-4dc6-930b-6a0e7552fb90
function plot_traces(data, info)
	tr = traces(data, info)
	P=[]
	#println("tr =$(tr[1])")
	texp = info["t_exp"]
	t = collect(1:length(tr[1]))*texp
	for i in 1:9
		push!(P, plot(t, tr[i],
			  label = "NAPH3 steps",
		      xlabel = "t",
		      ylabel = "Rate",
			  title = "Measured Traces",
		      legend = false))
	end
	plot(P...,  layout=(3,3), size=(700, 500))
end

# ╔═╡ c54c8180-4f47-4c7c-b6a9-949c1685502a
plot_traces(data, RunInfo)

# ╔═╡ e44059a0-311c-45e2-92b2-85d652f690e3
"""create a DataFrame with binned trace coordinates and corresponding traces"""
function df_traces(data, info)
    
    # Obtain traces 
    trs = traces(data, info)
    
    # Compute spatial resolution based on pixel scale and division factor
    res = info["pixelscale"] / info["pixeldiv"]
    
    # Calculate binned indices i and j from positions x and y (convert microns to indices)
    # Julia uses integer division operator `÷` (equivalent to Python's `//`)
    i_indices = Int.(data.x  .÷ info["res"])
    j_indices = Int.(data.y  .÷ info["res"])

    # Construct DataFrame with indices and corresponding traces
    TracesDf = DataFrame(i = i_indices, j = j_indices, traces = trs)

    # Filter out rows where indices exceed the provided limits (edges lost due to binning)
    TracesDf = filter(row -> row.i < info["NHR"][1] && row.j < info["NHR"][2], TracesDf)
    
    return TracesDf
end

# ╔═╡ ae4c38d8-0bd8-45a4-916f-78289754dd1a
dft = df_traces(data, RunInfo)

# ╔═╡ c7451c25-e798-47cf-99c0-90de0be61adc
""" Prepare a high resolution image in the CMOS

	- returns an image im[hrdim x hrdim] such that each pixel im[i,j] contains a trace.

	-NB, HR image comes in counts
"""
function hr_image(i, tracesDf, info)
	# Initialize HR image array with zeros
	im = zeros(Float64, info["NHR"][1], info["NHR"][2])  

	# Fill in the image from the trace data for molecule number "i"
	
	for row in eachrow(tracesDf)
		im[row.i + 1, row.j + 1] = row.traces[i + 1]
	end

	im
end

# ╔═╡ 9f383740-1721-4f86-8a21-bcc8bbf4ee28
	"""
	- Compute the low resolution obtained from the HR image 
	- Pass a gaussian filter to reflect diffractive limit

	To compute the low resolution image we do 2D block binning by summation. For each binning × binning block of the original image, we compute the sum and place it in a downsampled 2D result.

	Then, if the original image was, say 100×100 (this corresponds to HR, below the diffractive limit) and binning = 10 (for example), the result will be a 10×10 image where each pixel is the sum of a 10×10 block from the original image.

	- In julia this is done like this: 

	- im is a 2D array (image) of size NHR[1] × NHR[2]. The goal is to reduce its resolution by aggregating each binning × binning block into a single pixel.

	- The instruction: 
	reshaped = reshape(im, (NHR[1] ÷ binning, binning, NHR[2] ÷ binning, binning))

	Does the following:
	reshaping the 2D image into a 4D array:
	•	Dimension 1: number of vertical bins (rows ÷ binning)
	•	Dimension 2: rows inside each bin
	•	Dimension 3: number of horizontal bins (cols ÷ binning)
	•	Dimension 4: columns inside each bin

	- The instruction:
	im_binned = sum(sum(reshaped, dims=(2, 4)), dims=(2, 4))[:,:,1,1]

	This sums over:
	•	dims=2: within each vertical bin
	•	dims=4: within each horizontal bin

	Then we do a second sum to reduce those 2D axes completely.

	The resulting array has dimensions:
		•	(NHR[1] ÷ binning, 1, NHR[2] ÷ binning, 1)

	Finally:
	[:,:,1,1]

	This collapses the 4D result into a 2D image by taking only the non-summed axes (the binned image dimensions).
	
	"""
	function CPS(n, tracesDf, info)
		# Get HR image
	    im = hr_image(n, tracesDf, info)
		println("in cps")
		
		dlr = info["dl"] /2.0 # diffraction limit radius
		res = info["res"]     # resolution in the HR image = pixelscale/pixeldiv

		# Pass a gaussian filter to simulate the effect of difractive limit. 
	    sigma = uconvert(nm, dlr)/uconvert(nm, res)

		println("sigma = $(sigma)")
	    imgf = imfilter(im, Kernel.gaussian(sigma)) 
		
		kk = 0
		for i in 1:size(im, 1), j in 1:size(im, 2)
		    if im[i, j] != 0 && kk < 6
				println("im[$i, $j] = ", im[i, j])
		        println("imgf[$i, $j] = ", imgf[i, j])
				kk+=1
		    end
		end

		 # Reshape and sum to bin the image
	    binning = info["pixeldiv"]
	    NHR = info["NHR"]
	
	    reshaped = reshape(imgf, 
						   (NHR[1] ÷ binning, 
							binning, 
							NHR[2] ÷ binning,
							binning))

	    im_binned = sum(sum(reshaped, dims=(2, 4)), dims=(2, 4))[:,1,:,1]

		println("im binned")
		kk = 0
		for i in 1:size(im_binned, 1), j in 1:size(im_binned, 2)
		    if im_binned[i, j] != 0 && kk < 6
		        println("im_binned[$i, $j] = ", im_binned[i, j])
				kk+=1
		    end
		end

		# add the constant background bck and multiply by time of exposition
		# to obtain count per second
		println("bck = $(info["bck"])")
		println("t_exp = $(info["t_exp"])")
	    cps = im_binned .+ info["bck"] .* info["t_exp"]

		println("cps")
		kk=0
		for i in 1:size(cps, 1), j in 1:size(cps, 2)
		    if cps[i, j] > 200.0 && kk < 6
		        println("cps[$i, $j] = ", cps[i, j])
				kk+=1
		    end
		end
	    cps
	end

# ╔═╡ 966054f8-c7c6-48d9-a820-da1bb813148c
cps = CPS(1, dft, RunInfo)


# ╔═╡ c6ffec58-094e-4e32-bd2f-259e14a62214
heatmap(cps, colorbar=true, title="CPS", aspect_ratio=1)

# ╔═╡ f46bbd84-a0f8-4373-9758-81837374c3c2
"""Add shot noise due to electronics and dark current to CPS"""
	function CPS_noisy(cps, info, cam;  λ=550nm)
		ns = noise(cam, λ, cps, info["t_exp"])
		println(size(ns))
	    im_n = [rand(Normal(cps[i, j], 
							ns[i, j])) for i in axes(cps, 1), j in axes(cps, 2)]
    	im_n
	end

# ╔═╡ 3f715111-af88-454f-bdc2-814515f9c79b
cps_n = CPS_noisy(cps, RunInfo, of)

# ╔═╡ b4f15620-d868-4827-963c-9334b52f4303
heatmap(cps_n, colorbar=true, title="CPS_n", aspect_ratio=1)

# ╔═╡ 0c95b794-9794-4f34-8ae4-a88e1994483d
"""Generates the temporal evolution of a simulated image sequence over time, 
using trace data (tracesDf).

Each time point corresponds to a frame/image where intensities are computed using the CPS function.

This is done like this:

- im0 = CPS(0, TracesDf, info) 

Calls the function CPS to compute the image at time index 0 (first time frame.)

- data = [zeros(size(im0)), im0]

Initializes a list of 2D arrays (images).

	•	The first element is a black (zeroed) image of the same size as im0.

	•	The second is im0 itself.

	•	This is a workaround to build an array of consistent shapes — later the black frame will be removed.

- N = length(TracesDf[1, :Traces])

	• This gets the number of time points by checking the length of the trace vector in the first molecule.

- for i in 2:N-1
    push!(data, CPS(i, TracesDf, info))
end

	•	Iterates over time indices starting from 2 to N-1.

	•	For each time index i, it computes the corresponding image using CPS(...) and appends it to data.

- reduce(vcat, [reshape(d, 1, size(d,1), size(d,2)) for d in data])

Each 2D image d is reshaped to a 3D array with shape (1, height, width).

	•	This adds a time dimension in the first axis.

	•	vcat (vertical concatenation) stacks all these time frames along the first dimension, yielding a final 3D array:

	•	Shape: (T, H, W) — T time steps, image of height H and width W.
"""
function temporal_evolution_CPS(TracesDf, info)
	im0 = CPS(0, TracesDf, info)
	data = [zeros(size(im0)), im0]
	
	N = length(TracesDf[1, :traces])
	for i in 2:N-1
		push!(data, CPS(i, TracesDf, info))
	end
	return reduce(vcat, [reshape(d, 1, size(d,1), size(d,2)) for d in data])
end

# ╔═╡ e8394f9d-b752-4d42-97e6-51488d5afe06
cpd3d = temporal_evolution_CPS(dft, RunInfo)

# ╔═╡ 50d3c818-e3f9-47f1-b1cd-84961c7e5d11
size(cpd3d)

# ╔═╡ 1ee8e8fb-a757-458a-bbf8-a6c8559f8c01
heatmap(cpd3d[150,:,:], colorbar=true, title="CPS_n", aspect_ratio=1)

# ╔═╡ 4ddf742f-1a89-4da3-b8cd-7f81ca9e1fd8
	"""Add noise to temporal evolution"""
function temporal_evolution_CPS_noisy(TracesDf, info, c)
	im0 = CPS_noisy(CPS(0, TracesDf, info), info, c)
	data = [zeros(size(im0)), im0]
	
	N = length(TracesDf[1, :traces])
	for i in 2:N-1
		d2 = CPS_noisy(CPS(i, TracesDf, info), info, c)
		push!(data, d2)
	end
	return reduce(vcat, [reshape(d, 1, size(d,1), size(d,2)) for d in data])
end


# ╔═╡ 993688d6-d02a-49ca-8ffe-56393c61e9b7
cpd3dn = temporal_evolution_CPS_noisy(dft, RunInfo, of)

# ╔═╡ 38cfe310-33cd-4172-b2fb-3538a9f51f10
heatmap(cpd3dn[150,:,:], colorbar=true, title="CPS_n", aspect_ratio=1)

# ╔═╡ d3fc191a-1f48-4cc0-8cea-52f031a05efc



# ╔═╡ Cell order:
# ╠═34233420-116a-11f0-0e9c-75681c7ddb3f
# ╠═bf3dd6ae-7cf8-45b4-b407-f480cb9e7f30
# ╠═1b562d29-0419-4c3b-b699-271888b547a5
# ╠═5b04c22a-1a29-4d17-be3b-3a42fa1ecdf9
# ╠═e6f2b315-b5b0-4983-a1c0-4b396b337b97
# ╠═c7301bb6-e147-4998-883e-2fa1c542c14e
# ╠═3738ca4c-c6e2-4a4d-90f2-1ee3284b7c89
# ╠═934b4ebc-a4f4-4900-9d19-f8b04726a635
# ╠═bba8e102-398c-4191-9d88-5cbc7879d6cf
# ╠═21345bd9-5cd7-4775-805c-e1407072a616
# ╠═a650449d-6f90-4655-ae83-64fe4eb67b22
# ╠═41053e81-3b76-41f3-9f54-076518a64e98
# ╠═4badd96d-2461-4507-b771-1b4882cd239d
# ╠═8d8cc6ff-3006-4cd1-9308-15b4c449a998
# ╠═1d59335e-7c76-49c2-88f6-aef6e4312941
# ╠═a2076949-49f3-4964-a2f8-6b3d6c6303fb
# ╠═f163b606-952f-4d14-9e9c-e9d9ce8722fd
# ╠═f31cd8db-1d37-4921-9c60-971ef17a3f22
# ╠═cbfd9e9c-2faa-4334-8324-dd7351e681b6
# ╠═c3440909-a0c4-4e73-8bc2-953b2d26bbc3
# ╠═e67eb86d-ea64-4a5c-97c3-15fa67725b7a
# ╠═e34083c1-526f-46bb-b75b-4717ac600fb7
# ╠═58fe32b8-d516-40f0-a875-5cd4df780355
# ╠═10be9bea-b310-432a-a14e-f79a0d76b796
# ╠═90e48289-e186-4118-a7d7-f911cce4b50e
# ╠═53768fe9-bfe7-423f-bd5f-8631e11c38c2
# ╠═5ae75425-b152-4099-b2fd-d212a75cca51
# ╠═67c79039-c97e-4efb-a1d0-acb6b3667c49
# ╠═a850b3f4-03f7-4978-b245-54b2460f8825
# ╠═845833df-5a13-489c-b7a3-9250abf92c36
# ╠═8ce9609f-efb4-4760-a4c4-308551669558
# ╠═872458a6-89ae-4fa2-85b0-220bad7503f4
# ╠═c7d45989-5a77-4b05-805e-1ded0d9beb76
# ╠═9543ff4c-7e48-4009-9a79-812fc3ffffcb
# ╠═9b37c1b4-8528-4b5e-8a30-ca7bfebda9ad
# ╠═8ad1d34e-7125-4cc9-934c-11522d6a74df
# ╠═a8869a63-7cf4-4d68-86e6-963c0a3e74c9
# ╠═bbca2d6f-a5dc-458e-9488-8ff01e0d92ec
# ╠═3ae994ca-44ec-45bf-b23a-9a578b35c5e2
# ╠═55d8bfdf-8d6a-4d74-89d6-862fed0d2260
# ╠═7d2d9c04-7237-4177-8711-a67df460aec7
# ╠═5c24b381-88ce-469e-8311-8136725261c4
# ╠═be574290-a0c3-4cb1-990e-f1b8f66e1034
# ╠═18b12b7b-db7d-40e0-be9b-cbf5a6f77438
# ╠═c1cd36d1-e43a-476c-ac1c-833b11c576f9
# ╠═f683cfc7-9bf3-45e7-ab28-fb1532e829dc
# ╠═c54c8180-4f47-4c7c-b6a9-949c1685502a
# ╠═a95a294d-f91f-4ef6-b2df-ee6da0b69155
# ╠═448f82e4-95d5-4bbb-97c7-9d4269aad25a
# ╠═ae4c38d8-0bd8-45a4-916f-78289754dd1a
# ╠═966054f8-c7c6-48d9-a820-da1bb813148c
# ╠═c6ffec58-094e-4e32-bd2f-259e14a62214
# ╠═3f715111-af88-454f-bdc2-814515f9c79b
# ╠═b4f15620-d868-4827-963c-9334b52f4303
# ╠═e8394f9d-b752-4d42-97e6-51488d5afe06
# ╠═50d3c818-e3f9-47f1-b1cd-84961c7e5d11
# ╠═1ee8e8fb-a757-458a-bbf8-a6c8559f8c01
# ╠═993688d6-d02a-49ca-8ffe-56393c61e9b7
# ╠═38cfe310-33cd-4172-b2fb-3538a9f51f10
# ╠═9200aad7-9aa6-4dc6-930b-6a0e7552fb90
# ╠═8bcbc2c2-fb3e-458d-a643-d55d7e15173c
# ╠═523009ad-d221-4f60-a7ed-010fb20c0a09
# ╠═e44059a0-311c-45e2-92b2-85d652f690e3
# ╠═69964f32-021d-4f99-af05-b68d41dc89d0
# ╠═f76cb6da-3b6b-43d2-b020-bd5aff6acab7
# ╠═82f10b3c-9cdd-4a0a-81da-610e2160588c
# ╠═c7451c25-e798-47cf-99c0-90de0be61adc
# ╠═9f383740-1721-4f86-8a21-bcc8bbf4ee28
# ╠═f46bbd84-a0f8-4373-9758-81837374c3c2
# ╠═0c95b794-9794-4f34-8ae4-a88e1994483d
# ╠═4ddf742f-1a89-4da3-b8cd-7f81ca9e1fd8
# ╠═d3fc191a-1f48-4cc0-8cea-52f031a05efc
