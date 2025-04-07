module BoldLab
using Revise
#print_greeting() = print("Hello, world!")
include("setup.jl")
include("SampleGenerators.jl")
export FMolecule, fm_from_csv
export Excitation, intensity, photon_energy, total_power
export Fov, Laser, Objective, GaussianLaser, gf, gI, ccd
export photon_energy, delivered_energy, n_photons, n_photons_int, photon_density
export  diffraction_limit, geometrical_acceptance, transmission
export Sample_1D, Sample_2D, Sample_3D
#include("histos.jl")
#include("SimpleLogger.jl")
#include("AutoStepfinder.jl")
#using .SimpleLogger
#using .AutoStepfinder
#using .dffunctions
#export SimpleLogger
#export AutoStepfinder
#export dffunctions

end