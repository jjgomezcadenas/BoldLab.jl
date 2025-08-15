#!/usr/bin/env julia
"""
NAPH Analysis Script
Analyzes fluorescence microscopy data for NAPH molecules
Based on fitTracesData2.jl Pluto notebook
"""

# ═══════════════════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════════════════

using Pkg
Pkg.activate(ENV["JBoldLab"])
push!(LOAD_PATH, ENV["JBoldLab"] * "/src")

using Revise
using BoldLab
using SimpleLogger
using JStepFinder
using StepAnalysis
using LabStepAnalysis
using StepPreprocessing
using JPELT
using histos

using CSV
using DataFrames
using StatsPlots 
using LsqFit
using Statistics
using StatsBase
using Images
using SparseArrays

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

# Set logging level
set_log_level(SimpleLogger.WARN)

# Filter definitions
ActiveFilters = Dict(
    "Filter3"  => [500.0, 520.0],
    "Filter4"  => [524.0, 544.0],
    "Filter5"  => [554.0, 568.0],
    "Filter6"  => [576.0, 596.0],
    "Filter7"  => [605.0, 625.0],
    "Filter8"  => [633.0, 647.0],
    "Filter12" => [500.0, 800.0]
)

dirs = ["Evap_test_Fakes", "Evap_test_NAPH"]

fake_samples =["FakeA1_0.5monolayer",
               "FakeA1_pre",
               "FakeA2_0.25monolayer",
               "FakeA2_pre",
               "FakeA3_0.1monolayer",
               "FakeA3_pre"]
signal_samples = [
    "SE0_pre_A1",    # NAPH3SIL 1e-6M pre evaporation, measured at day 0 
    "SE1_ba0.1m_A1", # NAPH3SIL 1e-6M after evaporation of 0.1 monolayer, measured on day+1
    "SE0_pre_B1",    # NAPH3SIL 1e-7M pre evaporation, measured at day 0
    "SE1_ba0.1m_B1", # NAPH3SIL 1e-7M after evaporation of 0.1 monolayer, measured on day+1
    "SE0_pre_B2",    # NAPH3SIL 1e-7M pre evaporation, measured at day 0, (replicate 2)
    "SE1_ba0.1m_B2", # NAPH3SIL 1e-7M after evaporation of 0.1 monolayer, measured on day+1 (replicate 2)
    "SE0_pre_B3",    # NAPH3SIL 1e-7M pre evaporation, measured at day 0, (replicate 3)
    "SE1_ba0.2m_B3", # NAPH3SIL 1e-7M after evaporation of 0.2 monolayer, measured on day+1 (replicate 3)
    "SE0_pre_C1",    # NAPH3SIL 1e-8M pre evaporation, measured at day 0
    "SE1_ba0.1m_C1", # NAPH3SIL 1e-8M after evaporation of 0.1 monolayer, measured on day+1
    "SE0_pre_C2",    # NAPH3SIL 1e-8M pre evaporation, measured at day 0, (replicate 2)
    "SE1_ba0.1m_C2", # NAPH3SIL 1e-8M after evaporation of 0.1monolayer, measured at day 1, (replicate 2)
    "SE0_pre_C3",    # NAPH3SIL 1e-8M pre evaporation, measured at day 0, (replicate 3)
    "SE2_ba0.2m_C3", # NAPH3SIL 1e-8M after evaporation of 0.2monolayer, measured at day 2, (replicate 3)
]

# Analysis parameters - modify these as needed
const ANALYSIS_CONFIG = Dict(
    # Path configuration
    "base_dir" => "/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign",
    "week" => "Evap_test_Fakes",     # Week1, Week2, Week3 where initial names, refers to main folder data
    "day" => "FakeA1_0.5monolayer",  # Second level folder data
    "scan" => "Set1",                # Scan/field folder
    "filter" => "Filter3",           # Filter to analyze
    
    # Background subtraction
    "σ_background" => 10.0,      # Sigma for background averaging
    "nlf_background" => 5,       # Number of frames for background
    
    # Denoising
    "σ_denoise" => 1.0,          # Sigma for Gaussian denoising
    
    # Peak detection
    "peak_threshold" => 50.0,    # Threshold for peak detection
    "peak_dx" => 0,              # Margin x for edge exclusion
    "peak_dy" => 0,              # Margin y for edge exclusion
    
    # Analysis frame
    "analysis_frame" => 1,       # Frame to analyze for peaks
    
    # Fitting options
    "fit_PELT" => true,          # Perform PELT fitting
    "fit_EXP" => true,           # Perform exponential fitting
    
    # Output
    "save_data" => true,        # Save plots to files
    "output_dir" => "/Users/jjgomezcadenas/Projects/BoldLab/naph_analysis"  # Output directory
)

# ═══════════════════════════════════════════════════════════════════════════════
# Helper Functions
# ═══════════════════════════════════════════════════════════════════════════════

function list_subfolders(dir)
    [basename(e) for e in readdir(dir; join=true) if isdir(e) && !startswith(basename(e), ".")]
end

function find_filter(path, filters)
    ff = split(path, "/")[end]
    i = findfirst("Filter", ff)
    suffix = i === nothing ? "" : ff[i.start:end]
    suffix, filters[suffix]
end

function png_from_path(path)
    parts = split(path, '/')
    folder1 = parts[end - 1]
    folder2 = parts[end]
    string("stepAnalysis", "/", folder1, "_", folder2, ".png")
end


function count_mol(df; cmol="nmol")
    length(unique(df[:, cmol]))
end


function build_data_path(config)
    """Build hierarchical data path from configuration"""
    root_dir = config["base_dir"]
    
    # Handle week folder structure (level 1)
    week = config["week"]
    path_week = occursin("Week", week) ? joinpath(root_dir, week, "Data") : joinpath(root_dir, week)
    
    # Build day and scan paths (levels 2-3)
    path_scan = joinpath(path_week, config["day"], config["scan"])
    
    # Check for filter subfolder (level 4)
    subdirs4, _, _ = scan_level(path_scan)
    path_tiff = config["filter"] in subdirs4 ? joinpath(path_scan, config["filter"]) : path_scan
    
    println("\nAnalyzing data from: $path_tiff")
    return path_tiff
end



# ═══════════════════════════════════════════════════════════════════════════════
# Main Analysis
# ═══════════════════════════════════════════════════════════════════════════════

function main(config=ANALYSIS_CONFIG)
    println("\n" * "="^80)
    println("NAPH ANALYSIS")
    println("="^80)
    
    # Create output directory if needed
    if config["save_data"]
        mkpath(config["output_dir"])
    end

    # Label this run
    lbl = string(config["week"], "_", config["day"], "_", config["scan"], "_", config["filter"])

    # Build hierarchical data path
    path_tiff = build_data_path(config)
    
    # Get filter information
    flt, fltv = find_filter(path_tiff, ActiveFilters)
    println("Filter: $flt")
    println("Filter band (nm): $(fltv)")
    
    # Load image stack
    println("\nLoading image stack...")
    imxt = load_tif_stack_int16(path_tiff, pedestal=0.0)
    println("Stack shape: $(size(imxt))")
    
    # Get global statistics
    totalI, meanI, stdI = get_stats(imxt; bck=0.0)
    println("Total intensity: $(sum(totalI))")
    println("Mean intensity: $(mean(meanI))")
    
    # Compute background
    println("\nComputing background...")
    background = compute_background_from_stack(imxt; 
                                               σ=config["σ_background"], 
                                               nlf=config["nlf_background"])
    
    # Subtract background
    println("Subtracting background...")
    imgbsub = subtract_background_from_stack(imxt, background)
    
    # Denoise stack
    println("Denoising stack...")
    imgd = denoise_stack(imgbsub; σ=config["σ_denoise"])

    # Display frames and wait for user input
    println("\nDisplaying processed frames...")
    plot_frames(imgd, nscale=10)
    println("\nPress ENTER to continue with peak detection...")
    readline()
    
    # Detect peaks in specified frame 
    nframe = config["analysis_frame"]
    println("\nDetecting peaks in frame $nframe...")
    peaks = detect_local_maxima(imgd[:, :, nframe]; 
                                threshold=config["peak_threshold"], 
                                dx=config["peak_dx"], 
                                dy=config["peak_dy"])
    
    println("Found $(size(peaks, 1)) molecule candidates")

    # Create and display intensity histogram
    pthr = config["peak_threshold"]
    _, pisp = step_hist(peaks.intensity;
              nbins=20,
              xlim=(0.0, 2*maximum(peaks.intensity)),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
              title="Intensity of selected peaks with cut = $(pthr)")
    
    # Display histogram and wait for user input
    println("\nDisplaying peak intensity histogram...")
    plt_intensity = plot(pisp, size=(600, 300))
    display(plt_intensity)
   
    println("\nPress ENTER to continue with peak visualization...")
    readline()

    CSV.write(joinpath(config["output_dir"], lbl * "_peaks.csv"), peaks)
    println("Saved Peaks results to CSV")

    # Display detected peaks on image
    println("\nDisplaying detected peaks on frame $nframe...")
    p1 = heatmap(imgd[:, :, nframe]; 
        color = cgrad(:grays, rev = true),
        aspect_ratio = 1, 
        title = "Detected Peaks (Frame $nframe)")
    
    scatter!(p1, 
        peaks.j,               # x-axis (columns)
        peaks.i,               # y-axis (rows)
        zcolor = peaks.intensity,
        colorbar = true,
        marker = (:circle, 3),
        c = :viridis,
        label = "Peaks")
    
    display(p1)
    println("\nPress ENTER to continue with trace building...")
    readline()
    
    # Build traces
    println("\nBuilding traces...")
    TRZS = build_sparse_traces(imgd, peaks)
    
    # PELT fitting
    if config["fit_PELT"]
        println("\nPerforming PELT analysis...")
        pdf, pI, pJ, pDX, pFX, pSC = find_fit_pelt(TRZS, peaks)
        println("PELT fitted $(count_mol(pdf)) molecules")
        

        # Save PELT results
        CSV.write(joinpath(config["output_dir"], lbl * "_pelt_fits.csv"), pdf)
        println("Saved PELT results to CSV")
    end
    
    # Exponential fitting
    if config["fit_EXP"]
        println("\nPerforming exponential fitting...")
        dfe, VI, VJ, VDX, VFX = fit_traces_exponential(TRZS, peaks)
        println("Exponential fitted $(size(dfe, 1)) traces")
                
        # Save exponential results
        CSV.write(joinpath(config["output_dir"], lbl * "_exp_fits.csv"), dfe)
        println("Saved exponential results to CSV")
        
    end
    
    # Save summary statistics
    summary = Dict(
        "data_path" => path_tiff,
        "filter" => flt,
        "filter_band_nm" => fltv,
        "total_intensity" => sum(totalI),
        "n_frames" => size(imxt, 3),
        "n_peaks_detected" => size(peaks, 1),
        "analysis_frame" => nframe,
        "peak_threshold" => config["peak_threshold"]
    )
    
    if config["fit_PELT"]
        summary["n_molecules_pelt"] = count_mol(pdf)
    end
    
    if config["fit_EXP"]
        summary["n_traces_exp"] = size(dfe, 1)
    end
    
    # Write summary
    ofile = lbl * "_summary.txt"
    
    open(joinpath(config["output_dir"], ofile), "w") do io
        println(io, "NAPH Analysis Summary")
        println(io, "=" ^ 50)
        for (key, value) in summary
            println(io, "$key: $value")
        end
    end
    
    println("\n" * "="^80)
    println("Analysis complete!")
    println("Results saved to: $(config["output_dir"])")
    println("="^80)
    
    return summary
end

# ═══════════════════════════════════════════════════════════════════════════════
# Script execution
# ═══════════════════════════════════════════════════════════════════════════════

if abspath(PROGRAM_FILE) == @__FILE__
    # Run with default configuration
    main()
    
    # Or modify configuration for specific analysis:
    # custom_config = copy(ANALYSIS_CONFIG)
    # custom_config["week"] = "Week2"
    # custom_config["day"] = "W2_SVE0_1"
    # custom_config["filter"] = "Filter4"
    # main(custom_config)
end