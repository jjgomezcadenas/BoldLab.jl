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
using ImageFiltering, ImageView

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

# ╔═╡ f163b606-952f-4d14-9e9c-e9d9ce8722fd
md"""
## Define Setup
"""

# ╔═╡ cbfd9e9c-2faa-4334-8324-dd7351e681b6
begin
	lsr = Laser(450nm, 1mW)
	I = uconvert(Hz*cm^-2, photon_density(450nm, 1mW, 15*μm, 15*μm))
	I2 = uconvert(Hz*cm^-2, photon_density(450nm, 3.115mW, 15*μm, 15*μm))

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
	of = CMOS("ORCA_Flash4", 6.5μm, 2049, 1.0, 0.06Hz, 4, oflash4_eff)
	md"""
	#### CMOS: 
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
end

# ╔═╡ c7d45989-5a77-4b05-805e-1ded0d9beb76
md"""
Typical noise for 100 photons and 0.5 s exposition = $(noise(of, 550nm, 100, 0.5s))

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
end

# ╔═╡ 8ad1d34e-7125-4cc9-934c-11522d6a74df
begin 
	pbcycles = 5e+4 # photobleaching cycles (nγ abosorbed before pB)
	dkcycles = 1e+5 # dark cycles (nγ abosorbed before DT)
	tdark    =2s    # time in dark states before decaying to ground
	
	naph3fB = FBMolecule(naph3f, "naph3-free", pbcycles, (dkcycles, tdark))
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
	sigma  = 15μm # sigma of gaussian beam in sample
	center = (15μm,15μm) # center of laser beam
	dimensions = [30μm, 30μm] #effectise "laser light spot"

	#sigma  = 1μm # sigma of gaussian beam in sample
	#center = (0μm,0μm) # center of laser beam
	#dimensions = [1μm, 1μm] #effectise "laser light spot"
	
	lexc = LaserExcitation(lsr, sigma, center)
end

# ╔═╡ be574290-a0c3-4cb1-990e-f1b8f66e1034
begin
	N = 1000
	sample=Sample_3D(N,dimensions)
end

# ╔═╡ 5c24b381-88ce-469e-8311-8136725261c4
md"""
### Generate Data
- Generate a sample of $(N) molecules with (x,y) coordinates distributed in a spot of  $(dimensions). Each molecule is oriented with orientation defined by angles (theta, phi)
"""

# ╔═╡ 18b12b7b-db7d-40e0-be9b-cbf5a6f77438
md"""
- Generate "true data" for this sample.
- The sample is iluminated with a laser excitation of intensity $(lexc.I), through an objective of numerical aperture $(obj.NA), and a filter defined between $(bf.λmin) and $(bf.λmax)-
- The data DF includes the emission rate R (for a Molecuar efficiency MOef) and the molecular trajectories before photobleaching. 
"""

# ╔═╡ c1cd36d1-e43a-476c-ac1c-833b11c576f9
begin
	data=generate_data(naph3fB, lexc, sample, obj.NA, bf.λmin, bf.λmax)
end

# ╔═╡ 50569fd9-c673-4dc2-8bbb-69b3c000fa3b
RunInfo

# ╔═╡ 020341d2-9558-4c40-a993-2fff8b1bdffe
md"""
#### Create traces
"""

# ╔═╡ bfc31410-0571-4066-b80d-1acfbc0ef7ee
dft = df_traces(data, RunInfo)

# ╔═╡ f61499cf-3f2f-489b-b826-135e1e116916
begin
	nf = 100
	imx = hr_image(nf, dft, RunInfo)
	heatmap(imx, colorbar=true, title="HR image frame=$(nf) ", aspect_ratio=1)
end

# ╔═╡ 09f89004-4d5f-4094-8d1a-85dbdcdd560d
begin
	dlrx = RunInfo["dl"] /2.0 # diffraction limit radius
    resx = RunInfo["res"]     # resolution in the HR image = pixelscale/pixeldiv

    # Pass a gaussian filter to simulate the effect of difractive limit. 
    sigmax = uconvert(nm, dlrx)/uconvert(nm, resx)
	imgfx = imfilter(imx, reflect(Kernel.gaussian(sigmax)))
	heatmap(imgfx, colorbar=true, title="HR + Gaussian Filter frame =$(nf) ", aspect_ratio=1)
end

# ╔═╡ a3be3612-7edd-4f3c-92a3-bf92d30ce8ae
begin
	#using NPZ
	#npzwrite("img.npy", imgf_trimmed)
end

# ╔═╡ ef259a4e-d420-47c8-b908-6677c0babfba
imfx = frame2D(nf, dft, RunInfo)

# ╔═╡ c6ffec58-094e-4e32-bd2f-259e14a62214
heatmap(imfx, colorbar=true, title="Frame $(nf): HR + GF + base-noise ", aspect_ratio=1)

# ╔═╡ 4fa40384-9520-4f49-90a6-d43317e577aa
md"""
### F
"""

# ╔═╡ 3f715111-af88-454f-bdc2-814515f9c79b
f2dn =frame2Dn(imfx, RunInfo, of)

# ╔═╡ b4f15620-d868-4827-963c-9334b52f4303
heatmap(f2dn, colorbar=true, title="Frame $(nf): HR + GF + all-noise ", aspect_ratio=1)

# ╔═╡ e8394f9d-b752-4d42-97e6-51488d5afe06
f3d = frame3D(dft, RunInfo)

# ╔═╡ 50d3c818-e3f9-47f1-b1cd-84961c7e5d11
size(f3d)

# ╔═╡ 1ee8e8fb-a757-458a-bbf8-a6c8559f8c01
begin
	fnn = 100
	heatmap(f3d[100,:,:], colorbar=true, title="f3d number $(fnn)", aspect_ratio=1)
end

# ╔═╡ 993688d6-d02a-49ca-8ffe-56393c61e9b7
f3dn = frame3Dn(dft, RunInfo, of)

# ╔═╡ 38cfe310-33cd-4172-b2fb-3538a9f51f10
heatmap(f3dn[150,:,:], colorbar=true, title="noisy f3dn number $(fnn)", aspect_ratio=1)

# ╔═╡ 2817d6df-7553-4d93-86f8-c11f668e6749
md"""
## Plot functions
"""

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

# ╔═╡ a2076949-49f3-4964-a2f8-6b3d6c6303fb
plot_molecule(naph3f, naph3c)

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

# ╔═╡ 8bcbc2c2-fb3e-458d-a643-d55d7e15173c
function plot_real_traces(data, t_meas, t_resol)
	P=[]
	for i in 1:9
		tx, ty = real_trace(data[i,:], t_meas, t_resol)
		#println("tx =$(tx)")
		#println("ty =$(ty)")
		push!(P, plot(tx, ty,
			  label = "NAPH3 steps",
		      xlabel = "t",
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
	RunInfo["gaussf"]=true
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
	RunInfo["NHR"] = Int.( dimensions .÷ res .-1)
	md"""
	- Time of measurement = $(t_meas) 
	- Time of exposition = $(t_exp)
	- Divisions = $(t_div)
	- Time (resolution) = $(t_resol)
	- Diffractive limit = $(dl)
	- pixel size = $(of.pixelsize)
	- bin pixel size = $(of.binPixelsize)
	- pixel scale = $(uconvert(nm, pixelscale))
	- pixel division = $(pixeldiv)
	- resolution =$(uconvert(nm, res))
	- constant background = $(RunInfo["bck"])
	"""

end

# ╔═╡ 743f3cc4-98e9-4bc8-858e-54ce328aa628
res

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
# ╠═8d8cc6ff-3006-4cd1-9308-15b4c449a998
# ╠═1d59335e-7c76-49c2-88f6-aef6e4312941
# ╠═a2076949-49f3-4964-a2f8-6b3d6c6303fb
# ╠═f163b606-952f-4d14-9e9c-e9d9ce8722fd
# ╠═cbfd9e9c-2faa-4334-8324-dd7351e681b6
# ╠═5ae75425-b152-4099-b2fd-d212a75cca51
# ╠═67c79039-c97e-4efb-a1d0-acb6b3667c49
# ╠═c7d45989-5a77-4b05-805e-1ded0d9beb76
# ╠═9543ff4c-7e48-4009-9a79-812fc3ffffcb
# ╠═9b37c1b4-8528-4b5e-8a30-ca7bfebda9ad
# ╠═8ad1d34e-7125-4cc9-934c-11522d6a74df
# ╠═55d8bfdf-8d6a-4d74-89d6-862fed0d2260
# ╠═7d2d9c04-7237-4177-8711-a67df460aec7
# ╠═5c24b381-88ce-469e-8311-8136725261c4
# ╠═be574290-a0c3-4cb1-990e-f1b8f66e1034
# ╠═18b12b7b-db7d-40e0-be9b-cbf5a6f77438
# ╠═c1cd36d1-e43a-476c-ac1c-833b11c576f9
# ╠═f683cfc7-9bf3-45e7-ab28-fb1532e829dc
# ╠═50569fd9-c673-4dc2-8bbb-69b3c000fa3b
# ╠═c54c8180-4f47-4c7c-b6a9-949c1685502a
# ╠═a95a294d-f91f-4ef6-b2df-ee6da0b69155
# ╠═020341d2-9558-4c40-a993-2fff8b1bdffe
# ╠═bfc31410-0571-4066-b80d-1acfbc0ef7ee
# ╟─743f3cc4-98e9-4bc8-858e-54ce328aa628
# ╠═f61499cf-3f2f-489b-b826-135e1e116916
# ╠═09f89004-4d5f-4094-8d1a-85dbdcdd560d
# ╠═a3be3612-7edd-4f3c-92a3-bf92d30ce8ae
# ╠═ef259a4e-d420-47c8-b908-6677c0babfba
# ╠═c6ffec58-094e-4e32-bd2f-259e14a62214
# ╠═4fa40384-9520-4f49-90a6-d43317e577aa
# ╠═3f715111-af88-454f-bdc2-814515f9c79b
# ╠═b4f15620-d868-4827-963c-9334b52f4303
# ╠═e8394f9d-b752-4d42-97e6-51488d5afe06
# ╠═50d3c818-e3f9-47f1-b1cd-84961c7e5d11
# ╠═1ee8e8fb-a757-458a-bbf8-a6c8559f8c01
# ╠═993688d6-d02a-49ca-8ffe-56393c61e9b7
# ╠═38cfe310-33cd-4172-b2fb-3538a9f51f10
# ╠═2817d6df-7553-4d93-86f8-c11f668e6749
# ╠═4badd96d-2461-4507-b771-1b4882cd239d
# ╠═9200aad7-9aa6-4dc6-930b-6a0e7552fb90
# ╠═8bcbc2c2-fb3e-458d-a643-d55d7e15173c
