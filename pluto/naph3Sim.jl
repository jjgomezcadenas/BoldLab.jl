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

# ╔═╡ e07ef09a-de10-4b0a-9448-904dc433f133
md"""
## Laser setup
- Define a Laser of 1 mW power and 450 nm wavelength
- Check that the results using the functions with Units and the Excitation functions (no Units) is consistent.
"""

# ╔═╡ cbfd9e9c-2faa-4334-8324-dd7351e681b6
lsr = Laser(450nm, 1mW)

# ╔═╡ b36df06b-58e9-4d66-a516-b06390350c85
photon_energy(lsr.λ)

# ╔═╡ 21dff539-4fbd-4b64-aaa0-c9879cb74055
let
	sizex=30. #μm
	sizey=30.  #μm
	dimensions=(sizex,sizey)
	I=5e20 # Intensity at the center of the beam ( nphot cm⁻² s⁻¹)
	wl=450 # Wavelength (nm)
	sigma=15 # Waist of the beam (microns)
	center=dimensions ./2 # Center of the beam ([x0, y0] in microns)
	exc=Excitation(I,wl,sigma,center)
	photon_energy(exc) # in eV
	total_power(exc) # in mW
end

# ╔═╡ b4ecf41f-2b5c-45ae-9d9d-3695f0314c05


# ╔═╡ 522ad34e-80e8-4387-8979-d66978fc6693


# ╔═╡ c3440909-a0c4-4e73-8bc2-953b2d26bbc3
md"""
- photon density for 1 mW
- photon density for 3.115 mW
"""

# ╔═╡ e67eb86d-ea64-4a5c-97c3-15fa67725b7a
uconvert(Hz*cm^-2, photon_density(450nm, 1mW, 15*μm, 15*μm))

# ╔═╡ e34083c1-526f-46bb-b75b-4717ac600fb7
uconvert(Hz*cm^-2, photon_density(450nm, 3.115mW, 15*μm, 15*μm))

# ╔═╡ ed02c818-6224-489d-ae1d-426f322b7ce5
typeof(photon_density(450nm, 3.115mW, 15*μm, 15*μm))

# ╔═╡ f163b606-952f-4d14-9e9c-e9d9ce8722fd
md"""
## Physics
"""

# ╔═╡ 3ae994ca-44ec-45bf-b23a-9a578b35c5e2
naph3f.xs(450.0)

# ╔═╡ 1a9a3d2e-077d-415d-ad42-3f21f82359a7
"""A fluorescent molecule with fotobleaching and Dark States"""
struct FBMolecule{T}
	fm::FMolecule
    name::String 
    QY::Float64
    pdf::T       # Interpolated emission PDF function
    xs::T        # Interpolated absorption cross section function
    M::typeof(1.0Hz)  # Average excitation/de-excitation cycles before photobleaching
    Darks::Tuple{typeof(1.0Hz), typeof(1.0s)}  # excitation to dark and time in dark.

    function FBMolecule(fm::FMolecule, name, M, Darks)
        # Define a helper function to modify the interpolation functions.
        gxs(f) = (λnm -> f(λnm/nm))
        # Explicitly specify the type parameter T as the type of gxs(fm.pdf)
        return new{typeof(gxs(fm.pdf))}(fm, name, fm.QY, gxs(fm.pdf), gxs(fm.xs), 
										M, Darks)
    end
end

# ╔═╡ ab03da2f-2cef-47a2-a1ea-5f3ae2d372a1
cross_section(fm::FBMolecule, λ::Unitful.Length) = fm.xs(λ) * fm.QY * barn

# ╔═╡ a8869a63-7cf4-4d68-86e6-963c0a3e74c9
naph3fB = FBMolecule(naph3f, "naph3-free", 5e+4Hz, (1e+5Hz, 2.0s))

# ╔═╡ bbca2d6f-a5dc-458e-9488-8ff01e0d92ec
cross_section(naph3fB, 450*nm)/barn

# ╔═╡ 880538a4-9012-44e1-b255-916ec302f891
"""
Represents the Laser excitation 
"""
struct LaserExcitation
    I::typeof(1.0Hz*cm^-2)  # Intensity at the center of the beam 
    λ::Unitful.Length       # Wavelength 
    sigma::Unitful.Length   # Waist of the beam (microns)
    center::Tuple{Unitful.Length,Unitful.Length} # Center of the beam ([x0, y0]

	function LaserExcitation(laser::Laser, sigma::Unitful.Length, 
							 center::Tuple{Unitful.Length,Unitful.Length} )
		
		I = photon_density(laser.λ, laser.P, sigma, sigma)
		new(I, laser.λ, sigma, center)
	end
end

# ╔═╡ 827dae28-8969-49e5-9fdb-d17338b79c40
"""Computes the intensity at a given (x, y) (cm⁻² s⁻¹) """
function intensity(ex::LaserExcitation, x::Unitful.Length, y::Unitful.Length)
    ex.I *
         exp(-((x - ex.center[1]) / (sqrt(2)*ex.sigma))^2) *
         exp(-((y - ex.center[2]) / (sqrt(2)*ex.sigma))^2)
end

# ╔═╡ 7d2d9c04-7237-4177-8711-a67df460aec7
begin
	sigma=15μm 
	center=(sigma,sigma)
	lexc = LaserExcitation(lsr, sigma, center )
end

# ╔═╡ cc17739b-4881-4d75-a740-147bb0bcfae5
function Sample_3D(N, dimensions)
    # dimensions is assumed to be a tuple or vector like (dim_x, dim_y)
    x = rand(N) .* dimensions[1]
    y = rand(N) .* dimensions[2]
    theta = acos.(1 .- 2 .* rand(N))
    phi = 2 * π .* rand(N)
    return DataFrame(x = x, y = y, theta = theta, phi = phi)
end

# ╔═╡ 220c1c90-370f-4162-94d4-d1a782e81e6b
"""Compute the radiative rate"""
function rad_rate(fm::FBMolecule, exc::LaserExcitation, row)
    B = cross_section(fm, exc.λ)
    I = intensity(exc, row[:x], row[:y])
    B * sin(row[:theta])^2 * sin(row[:phi])^2 * I
end

# ╔═╡ 1227065f-461d-4314-8c05-5359fb2082df
"""Compute the MO efficiency"""
function mo_eff(NA, row)
    delta = asin(NA)
    tetha = row[:theta]
    sincosm = sin(tetha - delta) * cos(tetha - delta)
    sincosp = sin(tetha + delta) * cos(tetha + delta)
    #sincos0 = sin(delta) * cos(delta)  # computed but not used further
    return (2 * delta - sincosp + sincosm) / π / 2
end

# ╔═╡ 347166ba-c7b0-493b-a710-1912c98067a4
""" Compute a time trajector based on exponential waiting times"""
function trajectory(molecule::FBMolecule, r::typeof(1.0Hz))
    t = 0.0
    hist = []
    dpb = r / molecule.M             # rate/PB-rate 
    ddark = r / molecule.Darks[1]    # rate/PB-rate
    
    tdark = rand(Exponential(ddark))
    tpb = rand(Exponential(dpb))
    
    while tdark < tpb
        t += tdark
        push!(hist, t*1e+3ms)
       
        tground = rand(Exponential(1.0s / molecule.Darks[2]))
        t += tground
        push!(hist, t*1e+3ms)
        tdark = rand(Exponential(ddark))
        tpb = rand(Exponential(dpb))
    end
    t += tpb
    push!(hist, t*ms)
    hist
end

# ╔═╡ f5156dc2-ebf4-457b-bbfc-2a635826f499
function spectral_fraction(fb::FBMolecule, λmin::typeof(1.0nm), λmax::typeof(1.0nm))
	quadgk(fb.fm.pdf, λmin/nm, λmax/nm)[1]
end

# ╔═╡ d91507d2-dbaf-4d5a-a946-8c7dd523dfba
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

# ╔═╡ be574290-a0c3-4cb1-990e-f1b8f66e1034
begin
	N= 10
	dimensions=(2*sigma, 2*sigma)
	sample=Sample_3D(N,dimensions)
end

# ╔═╡ 908da7e1-6d9c-4d84-8840-e028aafacac7
sample[1,:]

# ╔═╡ 1f36d879-602f-46cd-b3d0-f6d4c074f9d7
rad_rate(naph3fB, lexc, sample[1,:])

# ╔═╡ da8e8737-69a3-40e9-984c-5727a7c9f12f
mo_eff(1.0, sample[1,:])

# ╔═╡ c1cd36d1-e43a-476c-ac1c-833b11c576f9
begin
	NA=0.95
	λmin=450.0nm
	λmax=800.0nm
	data=generate_data(naph3fB, lexc, sample, NA, λmin, λmax)
end

# ╔═╡ c2389741-5323-4988-b208-12fe659de175
data.x

# ╔═╡ Cell order:
# ╠═34233420-116a-11f0-0e9c-75681c7ddb3f
# ╠═bf3dd6ae-7cf8-45b4-b407-f480cb9e7f30
# ╠═1b562d29-0419-4c3b-b699-271888b547a5
# ╠═5b04c22a-1a29-4d17-be3b-3a42fa1ecdf9
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
# ╠═e07ef09a-de10-4b0a-9448-904dc433f133
# ╠═cbfd9e9c-2faa-4334-8324-dd7351e681b6
# ╠═b36df06b-58e9-4d66-a516-b06390350c85
# ╠═21dff539-4fbd-4b64-aaa0-c9879cb74055
# ╠═b4ecf41f-2b5c-45ae-9d9d-3695f0314c05
# ╠═522ad34e-80e8-4387-8979-d66978fc6693
# ╠═c3440909-a0c4-4e73-8bc2-953b2d26bbc3
# ╠═e67eb86d-ea64-4a5c-97c3-15fa67725b7a
# ╠═e34083c1-526f-46bb-b75b-4717ac600fb7
# ╠═ed02c818-6224-489d-ae1d-426f322b7ce5
# ╠═f163b606-952f-4d14-9e9c-e9d9ce8722fd
# ╠═ab03da2f-2cef-47a2-a1ea-5f3ae2d372a1
# ╠═827dae28-8969-49e5-9fdb-d17338b79c40
# ╠═bbca2d6f-a5dc-458e-9488-8ff01e0d92ec
# ╠═3ae994ca-44ec-45bf-b23a-9a578b35c5e2
# ╠═1a9a3d2e-077d-415d-ad42-3f21f82359a7
# ╠═a8869a63-7cf4-4d68-86e6-963c0a3e74c9
# ╠═880538a4-9012-44e1-b255-916ec302f891
# ╠═7d2d9c04-7237-4177-8711-a67df460aec7
# ╠═cc17739b-4881-4d75-a740-147bb0bcfae5
# ╠═220c1c90-370f-4162-94d4-d1a782e81e6b
# ╠═1227065f-461d-4314-8c05-5359fb2082df
# ╠═347166ba-c7b0-493b-a710-1912c98067a4
# ╠═f5156dc2-ebf4-457b-bbfc-2a635826f499
# ╠═d91507d2-dbaf-4d5a-a946-8c7dd523dfba
# ╠═be574290-a0c3-4cb1-990e-f1b8f66e1034
# ╠═908da7e1-6d9c-4d84-8840-e028aafacac7
# ╠═1f36d879-602f-46cd-b3d0-f6d4c074f9d7
# ╠═da8e8737-69a3-40e9-984c-5727a7c9f12f
# ╠═c1cd36d1-e43a-476c-ac1c-833b11c576f9
# ╠═c2389741-5323-4988-b208-12fe659de175
