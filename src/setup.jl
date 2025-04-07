using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
using Interpolations
using DataFrames
using CSV

import Unitful:
    nm, μm, mm, cm, m, km, 
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W,
    A, N, mol, mmol, V, L, M


"""
A struct representing a Fluorescence Molecule
"""
struct FMolecule{T}
    QY::Float64  # Quantum Yield
    em_spectrum::Tuple{Vector{Float64}, Vector{Float64}}  # emission spectrum pdf
    abs_spectrum::Tuple{Vector{Float64}, Vector{Float64}}   # absorption spectrum
    pdf::T       # Interpolated emission PDF function
    xs::T        # Interpolated absorption cross section function
    
end


"""
Creates the FM structs from the DF. Assumes that the DF columns are:
for dfe (emission data frame):
c1: lambda, c2; free emission; c3 chelated emissions.
for dfa (absorbtion data frame):
c1: lambda, c2; free abs; c3 chelated abs.
"""
function fm_from_csv(datadir, emfile, absfile)
	
	function trapz(x::Vector{Float64}, y::Vector{Float64})
		s = 0.0
		for i in 1:length(x)-1
			s += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
		end
		return s
	end

	function itp(λ, ϕ)
		area = trapz(λ, ϕ)
		y     = ϕ ./area
		LinearInterpolation(λ, y, extrapolation_bc=Line())
	end

	eps_to_xsec = 3.82E+3 # eps in M*1 cm^-1, xsec in barn (10^-24 cm2)
		
	fn3e = joinpath(datadir, emfile)
	fn3a = joinpath(datadir, absfile)
	dfe = CSV.read(fn3e, DataFrame)
	dfa = CSV.read(fn3a, DataFrame)

	ndfe = names(dfe)
	le = Float64.(dfe[!,ndfe[1]])
	emf = Float64.(dfe[!,ndfe[2]])
	emc = Float64.(dfe[!,ndfe[3]])
	
	iemf = itp(le, emf)
	iemc = itp(le, emc)

	ndfa = names(dfa)
	la = Float64.(dfa[!,ndfa[1]])
	ema = Float64.(dfa[!,ndfa[2]]) 
	emca = Float64.(dfa[!,ndfa[3]]) 
	
	iema = LinearInterpolation(la, ema * eps_to_xsec, extrapolation_bc=Line())
	iemca = LinearInterpolation(la, emca * eps_to_xsec, extrapolation_bc=Line())
	
	FMolecule(1.0, (le,emf), (la,ema), iemf, iema), FMolecule(1.0, (le,emc), (la,emca), iemc, iemca)

end


"""
	struct Fov

Represent a field of view

# Fields
- `d::Unitful.Length`  : diameter of Fov
- `z::Unitful.Length`  : thickness
- `a::Unitful.Area`    : area (computed)
- `v::Unitful.Volume`  : volume (computed)

"""
struct Fov
    d::Unitful.Length
    z::Unitful.Length
	a::Unitful.Area
    v::Unitful.Volume

	function Fov(d,z)
		a = π * (d/2.)^2
		v = a * z
		new(d,z,a,v)
	end
end


"""
	struct Laser

Simple representation of a laser

# Fields
- `λ::typeof(1.0nm)`  : Laser wavelength
- `P::typeof(1.0mW)`  : Power

"""
struct Laser
	λ::typeof(1.0nm)
	P::typeof(1.0mW)
end


"""
	struct Objective

Simple representation of a microscope objective

# Fields
- `name::String` : identifies the objective
- `NA::Float64`  : Numerical aperture
- `M::Float64`   : Magnification

"""
struct Objective
    name::String
    NA::Float64
    M::Float64
end


"""
	GaussianLaser

Representation of a Gaussian laser

# Fields
- `laser::Laser`       : A laser type
- `obj::Objective`     : An objective type
- `w0::typeof(1.0nm)`  : Waist of laser at focusing point
- `zr::typeof(1.0nm)`  : z of laser at focusing point

"""
struct GaussianLaser
	laser::Laser
	obj::Objective
	w0::typeof(1.0nm)
	zr::typeof(1.0nm)

	function GaussianLaser(laser,obj)
		w0 = laser.λ/(π * obj.NA)
        zr  = laser.λ/(π * obj.NA^2)
		new(laser,obj, w0, zr)
	end

    function GaussianLaser(laser,w0, zr)
		obj = Objective("widefield", 1.0, 1.0)
		new(laser,obj, w0, zr)
	end
end


"""
Represents the Laser excitation 
"""
struct Excitation
    I::Float64                    # Intensity at the center of the beam ( nphot cm⁻² s⁻¹)
    wl::Float64                   # Wavelength (nm)
    sigma::Float64                # Waist of the beam (microns)
    center::Tuple{Float64,Float64}  # Center of the beam ([x0, y0] in microns)
end



"""Computes the intensity at a given (x, y) (cm⁻² s⁻¹) """
function intensity(ex::Excitation, x::Float64, y::Float64)
    return ex.I *
           exp(-((x - ex.center[1]) / (sqrt(2)*ex.sigma))^2) *
           exp(-((y - ex.center[2]) / (sqrt(2)*ex.sigma))^2)
end


"""Computes the photon energy in eV"""
function photon_energy(ex::Excitation)
    return 1239.8 / ex.wl
end


"""Computes the total power in mW"""
function total_power(ex::Excitation)
    # Note: ex.sigma is in microns, so we convert it (σ/1e4) to cm
    return ex.I * 2 * π * (ex.sigma / 1e4)^2 * photon_energy(ex) * 1.6e-19 * 1e+3
end

#FUNCTIONS

"""
	photon_energy(λ::Unitful.Length)

Given wavelength of photon return its energy.
# Fields

- `λ::Unitful.Length`  : Photon wavelength

"""
function photon_energy(λ::Unitful.Length)
	uconvert(eV, λ, Spectral())
end


"""
	delivered_energy(laser::Laser, t::Unitful.Time)

Delivered energy of a laser in a given time.
# Fields

- `laser::Laser`     : Laser
- `t::Unitful.Time`  : Time in which target is illuminated

"""
function delivered_energy(laser::Laser, t::Unitful.Time)
	laser.P * t
end


"""
	n_photons(laser::Laser)

Rate of photons (number of photons per unit time) produced by a laser
# Fields

- `laser::Laser`     : Laser
- `t::Unitful.Time`  : Time in which target is illuminated

"""
function n_photons(laser::Laser)
	uconvert(Hz, laser.P / photon_energy(laser.λ))
end


"""
	n_photons(λ::Unitful.Length, p::Unitful.Power)

Rate of photons (number of photons per unit time) corresponding to a wavelength
λ and a power P

# Fields

- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power

"""
function n_photons(λ::Unitful.Length, p::Unitful.Power)
	uconvert(Hz, p / photon_energy(λ))
end


"""
	n_photons_int(laser::Laser, t::Unitful.Time)

Integrated number of photons in a given time emitted by a laser

# Fields

- `laser::Laser`    : Laser
- `t::Unitful.Time` : time of measurement

"""
function n_photons_int(laser::Laser, t::Unitful.Time)
	uconvert(eV,delivered_energy(laser, t)) / photon_energy(laser.λ)
end


"""
	photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area)

number of photons per unit time per unit area

# Fields

- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power
- `a::Unitful.Area`   : Area

"""
function photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area)
	return n_photons(λ, p)/ a
end

function photon_density(λ::Unitful.Length, p::Unitful.Power, 
                        σx::Unitful.Length, σy::Unitful.Length)
	return n_photons(λ, p)/ (2*π*σx*σy)
end


"""
	photon_density(l::Laser, fov::Fov)

number of photons per unit time per unit area, in a Fov illuminated by a laser

# Fields

- `laser::Laser` : Laser
- `fov::Fov`     : Field of view

"""
function photon_density(laser::Laser, fov::Fov)
	return n_photons(laser) / fov.a
end


function cw(gl::GaussianLaser, z::Unitful.Length)
	return 1.0 + (z / gl.zr)^2
end


function w0wz2(gl::GaussianLaser, z::Unitful.Length)
	return 1.0 / cw(gl,z)
end


"""
	w(gl::GaussianLaser, z::Unitful.Length)

Waist of a laser at length z

# Fields

- `gl::GaussianLaser`  : A gaussian laser
- `z::Unitful.Length`  : z distance from focusing point
"""
function w(gl::GaussianLaser, z::Unitful.Length)
	return gl.w0 * sqrt(cw(gl, z))
end


"""
	gf(gl::GaussianLaser, z::Unitful.Length, r::Unitful.Length)

Gaussian distribution of a gaussian beam

# Fields

- `gl::GaussianLaser`  : A gaussian laser
- `z::Unitful.Length`  : z distance from focusing point
- `r::Unitful.Length`  : r distance from focusing point
"""
function gf(gl::GaussianLaser, z::Unitful.Length, r::Unitful.Length)
	wz = w(gl, z)
	return w0wz2(gl, z) * exp(-2.0 * (r/wz)^2)
end


"""
	gI(gl::GaussianLaser, z::Unitful.Length, r::Unitful.Length)

Intensity of a gaussian beam

# Fields

- `gl::GaussianLaser`  : A gaussian laser
- `z::Unitful.Length`  : z distance from focusing point
- `r::Unitful.Length`  : r distance from focusing point
"""
function gI(gl::GaussianLaser, z::Unitful.Length, r::Unitful.Length)
	wz = w(gl, z)
	return gl.I * (gl.w0/wz)^2 * gf(gl,z,r)
end


"""
	diffraction_limit(l::Laser, obj:: Objective)

Return the diameter of diffractive spot for a laser l
focused with an objective obj

# Fields

- `l::Laser`         : A laser
- `obj:: Objective`  : An objective
"""
function diffraction_limit(l::Laser, obj:: Objective)
    return 1.22 * l.λ/(2 * obj.NA)
end



"""
	diffraction_limit(gl::GaussianLaser)

Return the diameter of diffractive spot for a gaussian laser

# Fields

- `gl::GaussianLaser`  : A gaussian laser
"""
function diffraction_limit(gl::GaussianLaser)
    return 1.22 * gl.laser.λ/(2 * gl.obj.NA)
end


"""
	geometrical_acceptance(d::Float64, D::Float64)

Compute the fraction of photons that make it through an iris
of diameter D located at a distance d from the emission point.

# Fields

- `d::Float64`   : distance between emission point and iris
- `D::Float64`   : Diameter of iris

"""
function geometrical_acceptance(d::Float64, D::Float64)
	return 0.5(1. - d/sqrt(d^2 + (D/2.)^2))
end



"""
	transmission(objective::Objective)

Compute the transmission of an objective (depends only of NA).

# Fields

- `objective::Objective` : Objective

"""
function transmission(objective::Objective)
	A = objective.NA
	A <=1 ? (1 - sqrt(1 - A^2)) /2 : 0.5
end


"""
	transmission(A::Float64)

Compute the transmission as a function of NA.

# Fields

- `A::Float64` : Numerical acceptance (NA)

"""
function transmission(A::Float64)
	A <=1 ? (1 - sqrt(1 - A^2)) /2 : 0.5
end

#CCD is defined in terms of a function which returns the CCD response
#(e.g, function ccd returns a response function, which gives efficiency
# as a function of wavelength)
"""
	ccd(lmin::Float64=350.0, lmax::Float64=1000.0)

Return the efficiency of a CCD as a function of wavelength.

# Fields

- `lmin::Float64=350.0` : Minimum wavelength for which efficiency is defined
- `lmax::Float64=350.0` : Maximum wavelength for which efficiency is defined

"""
function ccd(lmin::Float64=350.0, lmax::Float64=1000.0)
	function eff(l::Float64)
		if l < lmin || l > lmax
			return 0.
		else
			wl = 350.:50.:1000.
			ϵ = [0.3, 0.4,0.65,0.78,0.82,0.82,0.8,0.72,0.62,0.5,0.37,
			  0.24,0.12,0.07]
			e = CubicSplineInterpolation(wl, ϵ)
			return e(l)
		end
	end
	return eff
end