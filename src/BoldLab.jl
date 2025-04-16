module BoldLab
using Revise
#print_greeting() = print("Hello, world!")
include("SimpleLogger.jl")
using .SimpleLogger
include("dffunctions.jl")

include("setup.jl")
include("SampleGenerators.jl")
include("GenerateData.jl")
export FMolecule, FBMolecule, fm_from_csv
export Laser, Objective, ccd, LaserExcitation, intensity, photon_energy, total_power
#export Fov,  Objective, GaussianLaser, gf, gI, ccd
export photon_energy, delivered_energy, n_photons, n_photons_int, photon_density
export diffraction_limit, geometrical_acceptance, transmission
export Sample_1D, Sample_2D, Sample_3D
export cross_section, generate_data, trajectory
#include("histos.jl")

#include("AutoStepfinder.jl")
#
#using .AutoStepfinder
#using .dffunctions
#export SimpleLogger
#export AutoStepfinder
#export dffunctions

end