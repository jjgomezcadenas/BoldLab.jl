### A Pluto.jl notebook ###
# v0.20.19

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
	using NPZ
	import Measures 
	import TiffImages
end

# ╔═╡ e6f2b315-b5b0-4983-a1c0-4b396b337b97
using ImageFiltering, ImageView

# ╔═╡ 5b04c22a-1a29-4d17-be3b-3a42fa1ecdf9
import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W

# ╔═╡ 203c85b0-fe14-4c5b-9dc7-e2bea8642933
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
                :(eval(x) = $(Expr(:core, :eval))($name, x)),
                :(include(x) = $(Expr(:top, :include))($name, x)),
                :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
                :(include($path))))
    m
end

# ╔═╡ 388cf7a0-4ba1-4178-8b74-0650945ba28c
begin
	cmdir=joinpath(ENV["DATA"], "HD5t/precdr")
	pdir =joinpath(ENV["PROJECTS"], "BoldLab")
end

# ╔═╡ 8fe42754-595d-4448-b26e-805eb1dcca3e
begin
	jn = ingredients(string(pdir,"/src/BoldLab.jl"))
end

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

# ╔═╡ 139f39e1-176c-410a-91bc-62234b31153d
md"""
## Run parameters
"""

# ╔═╡ 7280a644-4a95-4e97-83a6-16c1a9332d7c
begin
	NM = 1.1e+3
	N=Int(NM)
	with_noise=true
	laser_power = 2W
	dimensions = [1e+3μm, 1e+3μm]  #This are the dimensions of the sample

	sigma  = 1e+3μm # sigma of gaussian beam in sample
	center = (1μm,1μm) # center of laser beam
	
	pbcycles = 1e+5 # photobleaching cycles (nγ abosorbed before pB)
	dkcycles = 5e+5 # dark cycles (nγ abosorbed before DT)
	tdark    =2s    # time in dark states before decaying to ground
	file = generate_filename(NM, laser_power, pbcycles, dkcycles, tdark; with_noise=with_noise)
end

# ╔═╡ 478baae3-f2e8-424d-8efc-ea532a8dbf95
md"""
- Number of molecules $(NM)
- Surface = $(dimensions[1] * dimensions[2])
- num of molecules/surface = $(NM/(dimensions[1] * dimensions[2]))
"""

# ╔═╡ f163b606-952f-4d14-9e9c-e9d9ce8722fd
md"""
## Define Setup
"""

# ╔═╡ cbfd9e9c-2faa-4334-8324-dd7351e681b6
begin
	lsr = Laser(450nm, laser_power)
	I = uconvert(Hz*cm^-2, photon_density(450nm, laser_power, sigma, sigma))
	#I2 = uconvert(Hz*cm^-2, photon_density(450nm, 3.115mW, 15*μm, 15*μm))

	md"""
	#### Laser: 
	- Wavelength = $(lsr.λ)
	- Power = $(lsr.P)
	- Photon energy = $(photon_energy(lsr.λ))
	- Intensity = $(I)
	"""
end

# ╔═╡ 5ae75425-b152-4099-b2fd-d212a75cca51
begin
	obj = Objective("highNA", 0.95, 100.0)
	bf = BandFilter(450.0nm, 800.0nm)
	md"""
	#### Objective and filter: 
	- NA = $(obj.NA)
	- Magnification = $(obj.M))
	- Filter λmin = $(bf.λmin)
	- Filter λmax = $(bf.λmax)
	"""
end

# ╔═╡ 67c79039-c97e-4efb-a1d0-acb6b3667c49
begin
	of = CMOS("ORCA_Flash4", 6.5μm, 2049, 1.0, 0.06Hz, 60, oflash4_eff)
	ofn = CMOS("ORCA_Flash4-Nominal", 6.5μm, 2049, 1.0, 0.06Hz, 1, oflash4_eff)
	md"""
	#### CMOS: 
	- name = $(of.name)
	- pixel size = $(ofn.pixelsize)
	- number of pixels = $(ofn.npixel)
	- readout noise = $(ofn.readout_noise) p.e.
	- dark current = $(ofn.dark_current)  p.e./pixel/s
	- binning = $(ofn.binning)
	- QE @ 550 nm = $(ofn.qe(550nm))
	- sensor size = $(ofn.sensorsize)
	- logical pixel size = $(ofn.binPixelsize)
	- effective number of pixels = $(ofn.binNpixel)

	- Notice that I also define a "of" camera to simulate the scaled binning of the small sample being considered
	"""
end

# ╔═╡ c7d45989-5a77-4b05-805e-1ded0d9beb76
md"""
Typical noise per photon and second exposition = $(noise(ofn, 550nm, 1, 1.0s))

"""

# ╔═╡ 9543ff4c-7e48-4009-9a79-812fc3ffffcb
md"""
#### RunInfo
- Define miscelanous run parameters in dictionary RunInfo
"""

# ╔═╡ 9b37c1b4-8528-4b5e-8a30-ca7bfebda9ad
begin
	RunInfo = Dict()
	RunInfo["t_div"] = 5 # number of divisions on time-exp for HR calculation 
	RunInfo["t_meas"] = 100s #Measurement time (length of run)
	RunInfo["t_exp"] = 0.5s #time of exposition per step
	RunInfo["t_resol"] = RunInfo["t_exp"] / RunInfo["t_div"] # time resolution
end

# ╔═╡ 8ad1d34e-7125-4cc9-934c-11522d6a74df
begin 
	
	naph3ff = FBMolecule(naph3f, "naph3-free", pbcycles, (dkcycles, tdark))
	md"""
	#### Define the physical molecule:
	- name = $(naph3ff.name)
	- QY =$(naph3ff.QY)
	- M = $(naph3ff.M)
	- Darks = $(naph3ff.Darks)
	- cross section at $(lsr.λ)) = $(cross_section(naph3ff, lsr.λ)/barn) barn
	"""
end

# ╔═╡ 945d715b-75a4-4e87-81a4-42eb74a0d008
begin 
	
	naph3fB = FBMolecule(naph3c, "naph3-ba", pbcycles, (dkcycles, tdark))
	md"""
	#### Define the physical molecule:
	- name = $(naph3fB.name)
	- QY =$(naph3fB.QY)
	- M = $(naph3fB.M)
	- Darks = $(naph3fB.Darks)
	- cross section at $(lsr.λ)) = $(cross_section(naph3fB, lsr.λ)/barn) barn
	"""
end

# ╔═╡ 55d8bfdf-8d6a-4d74-89d6-862fed0d2260
md"""
#### Define Laser Excitation 
"""

# ╔═╡ 7d2d9c04-7237-4177-8711-a67df460aec7
begin
	
	

	#sigma  = 1e+3μm # sigma of gaussian beam in sample
	#center = (1e+3μm,1e+3μm) # center of laser beam
	
	
	
	lexc = LaserExcitation(lsr, sigma, center)
end

# ╔═╡ 5c24b381-88ce-469e-8311-8136725261c4
md"""
### Generate Data
- Generate a sample of $(N) molecules with (x,y) coordinates distributed in a spot of  $(dimensions). Each molecule is oriented with orientation defined by angles (theta, phi)
"""

# ╔═╡ be574290-a0c3-4cb1-990e-f1b8f66e1034
begin
	sample=Sample_3D(N,dimensions)
end

# ╔═╡ 47f7dd48-5562-4dcd-9125-561526958b9d
begin 
	
	scatter(sample.x, sample.y)
end

# ╔═╡ 18b12b7b-db7d-40e0-be9b-cbf5a6f77438
md"""
- Generate integrated signal data for this sample.
- The sample is iluminated with a laser excitation of intensity $(lexc.I), through an objective of numerical aperture $(obj.NA), and a filter defined between $(bf.λmin) and $(bf.λmax)
- The data DF includes the emission rate R (for a Molecular efficiency MOeff) and spectral fraction Spx, without computing trajectories.
"""

# ╔═╡ c1cd36d1-e43a-476c-ac1c-833b11c576f9
begin
	data=generate_data_integrated_signal(naph3fB, lexc, sample, obj.NA, bf.λmin, bf.λmax)
end

# ╔═╡ cd1b2041-e291-4519-993d-df05c9ce3d26
begin
htp2, ptp2 = jn.BoldLab.step_hist(data.R/Hz .* data.MOeff;
                 nbins = 50,
                 xlim   = (0.0, 1500.0),
                 xlabel = " N photons",
                 ylabel = "Frequency",
                 title=" photons emitted ")
	ptp2
end

# ╔═╡ 1832df91-a3cf-427c-9c9e-d4f4e6496635
1e+5/2/1e+4

# ╔═╡ ee09e501-77bc-4b37-89e9-dabeedb4b743
htp2

# ╔═╡ 4eaeb791-ac84-47db-a1c7-9eb27713846a


# ╔═╡ fe2e3696-80ab-4eb2-8523-188ed3fd88bd
Lf = mean(data.R/Hz.* data.MOeff)

# ╔═╡ 3a72a34d-7403-4679-a235-f4ed8379fdbd
stn(Lf, m, eps) = Lf/sqrt(Lf * m * eps)

# ╔═╡ e188d97f-05e2-4887-bb7f-8bf0e05c06d5
stn(Lf*100.0, 1e+4, 1e-2)

# ╔═╡ 265b2789-70e9-4de7-93bd-bf1402e2171f
function total_signal_in_sample(Lf; nmol=1e+6)
	return Lf * nmol
end

# ╔═╡ 5d0b17d5-fe5c-48f4-a6de-a1b209c281b3
ts = total_signal_in_sample(Lf*100.0, nmol=1e+3) 

# ╔═╡ 81beea3a-efc4-4da3-b320-6eda461d09ff
stn(ts, 1e+4, 1e-2)

# ╔═╡ 8a15e991-6a2e-4f2e-a5e7-f248a7b5b36f


# ╔═╡ a9639e54-98a7-4293-a61c-ce4b93cbbdfc
ptrap(r, rho, F) = 1.0 - exp(-π*r^2*F)


# ╔═╡ 2f72900d-7714-4178-adb6-1289d642dc14
ptrap(1.0, 1.0/100.0, 5.0)

# ╔═╡ 839eb0ca-51cf-4889-ab19-e36661d1e7bc
md"""
- mean signal per molecule and second: $Lf
- total signal in sample per second : $ts
"""

# ╔═╡ aa0e0474-7d40-4261-aca7-65ee0e5ecbf0
2615-2457

# ╔═╡ 020341d2-9558-4c40-a993-2fff8b1bdffe
md"""
#### Create Integrated Signal Heatmap
"""

# ╔═╡ a95a294d-f91f-4ef6-b2df-ee6da0b69155
begin
	dl = ( (bf.λmin + bf.λmax) / 2 ) / (2 * obj.NA) #diffractive limit
	pixelscale=of.binPixelsize/obj.M
	pixeldiv=4
	res=pixelscale/pixeldiv
	RunInfo["pixelscale"] = pixelscale
	RunInfo["pixeldiv"] = pixeldiv
	RunInfo["res"] = res
	RunInfo["dl"] = dl
	RunInfo["NHR"] = Int.( dimensions .÷ res .-1)
	imaging_params_ready = true
	md"""
	- Diffractive limit = $(dl)
	- pixel size = $(of.pixelsize)
	- bin pixel size = $(of.binPixelsize)
	- pixel scale = $(uconvert(nm, pixelscale))
	- pixel division = $(pixeldiv)
	- resolution =$(uconvert(nm, res))
	- Image dimensions (NHR) = $(RunInfo["NHR"])
	"""
end

# ╔═╡ bfc31410-0571-4066-b80d-1acfbc0ef7ee
begin
	# Create integrated signal image
	if imaging_params_ready
		# Initialize image
		integrated_signal = zeros(Float64, RunInfo["NHR"][1], RunInfo["NHR"][2])

		# Calculate binned indices for each molecule
		for row in eachrow(data)
			i_idx = Int(row.x ÷ res)
			j_idx = Int(row.y ÷ res)

			# Check bounds
			if i_idx < RunInfo["NHR"][1] && j_idx < RunInfo["NHR"][2]
				# Integrated signal = R * MOeff * Spx
				signal = (row.R / Hz) * row.MOeff * row.Spx
				integrated_signal[i_idx + 1, j_idx + 1] += signal
			end
		end

		integrated_signal
	end
end

# ╔═╡ f61499cf-3f2f-489b-b826-135e1e116916
begin
	# Apply Gaussian filter to simulate diffraction limit
	dlr = RunInfo["dl"] / 2.0
	sigma_int = uconvert(nm, dlr) / uconvert(nm, res)

	integrated_signal_filtered = imfilter(integrated_signal, reflect(Kernel.gaussian(sigma_int)))

	heatmap(integrated_signal_filtered,
		colorbar=true,
		title="Integrated Signal (R × MOeff × Spx) with diffraction blur",
		aspect_ratio=1,
		xlabel="x (pixels)",
		ylabel="y (pixels)")
end

# ╔═╡ 314a2753-1248-4b10-94b9-7192c9b709cf
begin
	lms = mean(integrated_signal_filtered)
	lts = sum(integrated_signal_filtered)
	md"""
	- mean signal per pixel = $lms
	- integrated signal $lts
	"""
end

# ╔═╡ 651d02d2-e7bb-4316-b5c1-ac1001e3905b
begin
	epsilon = 0.01
	sn = 5
	lb = ls * epsilon
	nmb = ls^2/(lb*sn^2)
end

# ╔═╡ 09f89004-4d5f-4094-8d1a-85dbdcdd560d
begin
	# Bin the high-resolution integrated signal image
	function rebin_integrated(M, binning)
		NHR = size(M)
		rows = (NHR[1] ÷ binning) * binning
		cols = (NHR[2] ÷ binning) * binning

		Mt = M[1:rows, 1:cols]
		M_rebinned = reshape(Mt, binning, size(Mt,1) ÷ binning, binning, size(Mt,2) ÷ binning)
		M_rebinned = permutedims(M_rebinned, (2, 4, 1, 3))
		M_binned = sum(M_rebinned, dims=(3, 4))
		dropdims(M_binned, dims=(3, 4))
	end

	binned_signal = rebin_integrated(integrated_signal_filtered, pixeldiv)

	heatmap(binned_signal,
		colorbar=true,
		title="Binned Integrated Signal (binning=$(pixeldiv))",
		aspect_ratio=1,
		xlabel="x (binned pixels)",
		ylabel="y (binned pixels)")
end

# ╔═╡ 1cace43c-2f6a-4a50-ba65-23b78eb55faa
mean(binned_signal)

# ╔═╡ 743f3cc4-98e9-4bc8-858e-54ce328aa628
res

# ╔═╡ Cell order:
# ╠═34233420-116a-11f0-0e9c-75681c7ddb3f
# ╠═bf3dd6ae-7cf8-45b4-b407-f480cb9e7f30
# ╠═1b562d29-0419-4c3b-b699-271888b547a5
# ╠═5b04c22a-1a29-4d17-be3b-3a42fa1ecdf9
# ╠═e6f2b315-b5b0-4983-a1c0-4b396b337b97
# ╠═203c85b0-fe14-4c5b-9dc7-e2bea8642933
# ╠═388cf7a0-4ba1-4178-8b74-0650945ba28c
# ╠═8fe42754-595d-4448-b26e-805eb1dcca3e
# ╠═c7301bb6-e147-4998-883e-2fa1c542c14e
# ╠═3738ca4c-c6e2-4a4d-90f2-1ee3284b7c89
# ╠═934b4ebc-a4f4-4900-9d19-f8b04726a635
# ╠═bba8e102-398c-4191-9d88-5cbc7879d6cf
# ╠═21345bd9-5cd7-4775-805c-e1407072a616
# ╠═a650449d-6f90-4655-ae83-64fe4eb67b22
# ╠═41053e81-3b76-41f3-9f54-076518a64e98
# ╠═8d8cc6ff-3006-4cd1-9308-15b4c449a998
# ╠═1d59335e-7c76-49c2-88f6-aef6e4312941
# ╠═a2076949-49f3-4964-a2f8-6b3d6c6303fb
# ╠═139f39e1-176c-410a-91bc-62234b31153d
# ╠═7280a644-4a95-4e97-83a6-16c1a9332d7c
# ╠═478baae3-f2e8-424d-8efc-ea532a8dbf95
# ╠═f163b606-952f-4d14-9e9c-e9d9ce8722fd
# ╠═cbfd9e9c-2faa-4334-8324-dd7351e681b6
# ╠═5ae75425-b152-4099-b2fd-d212a75cca51
# ╠═67c79039-c97e-4efb-a1d0-acb6b3667c49
# ╠═c7d45989-5a77-4b05-805e-1ded0d9beb76
# ╠═9543ff4c-7e48-4009-9a79-812fc3ffffcb
# ╠═9b37c1b4-8528-4b5e-8a30-ca7bfebda9ad
# ╠═8ad1d34e-7125-4cc9-934c-11522d6a74df
# ╠═945d715b-75a4-4e87-81a4-42eb74a0d008
# ╠═55d8bfdf-8d6a-4d74-89d6-862fed0d2260
# ╠═7d2d9c04-7237-4177-8711-a67df460aec7
# ╠═5c24b381-88ce-469e-8311-8136725261c4
# ╠═be574290-a0c3-4cb1-990e-f1b8f66e1034
# ╠═47f7dd48-5562-4dcd-9125-561526958b9d
# ╠═18b12b7b-db7d-40e0-be9b-cbf5a6f77438
# ╠═c1cd36d1-e43a-476c-ac1c-833b11c576f9
# ╠═cd1b2041-e291-4519-993d-df05c9ce3d26
# ╠═1832df91-a3cf-427c-9c9e-d4f4e6496635
# ╠═ee09e501-77bc-4b37-89e9-dabeedb4b743
# ╠═4eaeb791-ac84-47db-a1c7-9eb27713846a
# ╠═fe2e3696-80ab-4eb2-8523-188ed3fd88bd
# ╠═3a72a34d-7403-4679-a235-f4ed8379fdbd
# ╠═e188d97f-05e2-4887-bb7f-8bf0e05c06d5
# ╠═265b2789-70e9-4de7-93bd-bf1402e2171f
# ╠═5d0b17d5-fe5c-48f4-a6de-a1b209c281b3
# ╠═81beea3a-efc4-4da3-b320-6eda461d09ff
# ╠═8a15e991-6a2e-4f2e-a5e7-f248a7b5b36f
# ╠═a9639e54-98a7-4293-a61c-ce4b93cbbdfc
# ╠═2f72900d-7714-4178-adb6-1289d642dc14
# ╠═839eb0ca-51cf-4889-ab19-e36661d1e7bc
# ╠═aa0e0474-7d40-4261-aca7-65ee0e5ecbf0
# ╠═020341d2-9558-4c40-a993-2fff8b1bdffe
# ╠═a95a294d-f91f-4ef6-b2df-ee6da0b69155
# ╠═bfc31410-0571-4066-b80d-1acfbc0ef7ee
# ╠═f61499cf-3f2f-489b-b826-135e1e116916
# ╠═314a2753-1248-4b10-94b9-7192c9b709cf
# ╠═651d02d2-e7bb-4316-b5c1-ac1001e3905b
# ╠═09f89004-4d5f-4094-8d1a-85dbdcdd560d
# ╠═1cace43c-2f6a-4a50-ba65-23b78eb55faa
# ╟─743f3cc4-98e9-4bc8-858e-54ce328aa628
