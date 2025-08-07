module BoldLab
using Revise
include("SimpleLogger.jl")
using .SimpleLogger
include("dffunctions.jl")
include("util.jl")
include("setup.jl")
include("plot_functions.jl")
include("SampleGenerators.jl")
include("GenerateData.jl")
include("analysis_functions.jl")
export to_fstr, vect_to_fstr, vect_to_list
export find_max
export detect_local_maxima, fit_traces_exponential
export generate_filename
export plot_traces, plot_frames, plot_traces_stats, plot_stats, plot_traces_h2d
export FMolecule, FBMolecule, fm_from_csv, plot_molecule, plot_molecule_emission, compute_filter_coverage
export Laser, Objective, ccd, LaserExcitation, BandFilter, intensity, photon_energy, total_power
export CMOS, oflash4_eff, nphe, noise
#export Fov,  Objective, GaussianLaser, gf, gI, ccd
export photon_energy, delivered_energy, n_photons, n_photons_int, photon_density
export diffraction_limit, geometrical_acceptance, transmission
export Sample_1D, Sample_2D, Sample_3D
export cross_section, generate_data, trajectory
export traces, df_traces, real_trace, measured_trace, hr_image
export hr_image, frame2D, frame2Dn, frame3D, frame3Dn
#include("histos.jl")

#include("AutoStepfinder.jl")
#
#using .AutoStepfinder
#using .dffunctions
#export SimpleLogger
#export AutoStepfinder
#export dffunctions

end