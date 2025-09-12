module BoldLab
using Revise
include("SimpleLogger.jl")
using .SimpleLogger
include("dffunctions.jl")
include("util.jl")
include("setup.jl")
#include("SampleGenerators.jl")
#include("GenerateData.jl")
include("analysis_functions.jl")
include("LabStepAnalysis.jl")
using .LabStepAnalysis
include("JPELT.jl")
using .JPELT
include("StepAnalysis.jl")
using .StepAnalysis
include("StepPreprocessing.jl")
using .StepPreprocessing
include("fit_functions.jl")
include("plot_functions.jl")
export to_fstr, vect_to_fstr, vect_to_list
export find_max
export detect_local_maxima
export load_tif_stack_int16, get_tif_image, tif_to_matrix, regularize_img, regularize_stack!
export generate_filename
export plot_traces, plot_frames, plot_traces_stats, plot_stats, plot_traces_h2d, plot_exponential_fit, plot_pelt_fit, plot_asf_fit, plot_crops_results
export fit_pelt, PeltFit, fitpelt, pelt_fit
export fit_asf, AsfFit
export fit_traces_exponential, ExponentialFit, chi2_dof, elbow
export get_stats, pixel_trace, traces_stats, traces_h2d, traces_thr
export plot_steps, plot_trx, fit_data, fit_data3, fit_traces, getsteps, plotfit, plotfit2
export traces_above_thr, build_traces, renumber_nmol!
export subtract_background_from_stack, compute_background_from_stack, denoise_stack
export filter_regions_in_stack, detections_dataframe
export build_sparse_traces, unique_i_j_across_t
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