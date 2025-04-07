using Interpolations
using JLD2  # For loading/saving data
using FileIO

mutable struct FMolecule
    name::String
    cs::Float64  # Cross section
    QY::Float64  # Quantum Yield
    em_spectrum::Tuple{Vector{Float64}, Vector{Float64}}  # emission spectrum
    abs_spectrum::Tuple{Vector{Float64}, Vector{Float64}} # absorption spectrum
    M::Float64
    Darks::Vector{Any}
    CS_spectrum::Tuple{Vector{Float64}, Vector{Float64}}
end

# Constructor
function Molecule(name, cs, QY, em_spectrum, abs_spectrum, M, Darks=[])
    CS_spectrum = (abs_spectrum[1], abs_spectrum[2] ./ maximum(abs_spectrum[2]) .* cs)
    return Molecule(name, cs, QY, em_spectrum, abs_spectrum, M, Darks, CS_spectrum)
end

# Methods for saving/loading
function save(m::Molecule, path::String)
    @save joinpath(path, "$(m.name).jld2") m
end

function load(path::String)
    data = load(path)
    return data["m"]
end

# Interpolated cross-section
function cross_section(m::Molecule, wl)
    interp = LinearInterpolation(m.CS_spectrum[1], m.CS_spectrum[2], extrapolation_bc=Line())
    return interp(wl)
end

# Absorption spectrum (normalized or not)
function abs_spectrum(m::Molecule, wl, norm=true)
    y = norm ? m.abs_spectrum[2] ./ maximum(m.abs_spectrum[2]) : m.abs_spectrum[2]
    interp = LinearInterpolation(m.abs_spectrum[1], y, extrapolation_bc=Line())
    return interp(wl)
end

# Emission spectrum (normalized or not)
function em_spectrum(m::Molecule, wl, norm=true)
    y = norm ? m.em_spectrum[2] ./ maximum(m.em_spectrum[2]) : m.em_spectrum[2]
    interp = LinearInterpolation(m.em_spectrum[1], y, extrapolation_bc=Line())
    return interp(wl)
end

# Spectral fraction
function spectral_fraction(m::Molecule, wl0, wl1)
    x_total = 100:0.1:1000
    interp_total = LinearInterpolation(m.em_spectrum[1], m.em_spectrum[2], extrapolation_bc=Line())
    y_total = interp_total.(x_total)

    x_window = wl0:0.1:wl1
    y_window = interp_total.(x_window)

    return sum(y_window) / sum(y_total)
end