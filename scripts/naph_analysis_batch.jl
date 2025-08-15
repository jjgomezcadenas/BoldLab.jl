#!/usr/bin/env julia
"""
NAPH Batch Analysis Script
Processes multiple datasets without interactive plots
Runs analysis over directories, samples, and filters
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

# Set logging level
set_log_level(SimpleLogger.WARN)

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

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

# Directories and samples
dirs = ["Evap_test_Fakes", "Evap_test_NAPH"]

fake_samples = [
    "FakeA1_0.5monolayer",
    "FakeA1_pre",
    "FakeA2_0.25monolayer",
    "FakeA2_pre",
    "FakeA3_0.1monolayer",
    "FakeA3_pre"
]

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


# Filters to process
filters_to_process = ["Filter3", "Filter4", "Filter5", "Filter6", "Filter7", "Filter8"]

# Base configuration template
const BASE_CONFIG = Dict(
    # Path configuration (will be modified in loops)
    "base_dir" => "/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign",
    "scan" => "Set1",  # Fixed scan folder
    
    # Background subtraction
    "σ_background" => 10.0,
    "nlf_background" => 5,
    
    # Denoising
    "σ_denoise" => 1.0,
    
    # Peak detection
    "peak_threshold" => 50.0,
    "peak_dx" => 0,
    "peak_dy" => 0,
    
    # Analysis frame
    "analysis_frame" => 1,
    
    # Fitting options
    "fit_PELT" => true,
    "fit_EXP" => true
)

# ═══════════════════════════════════════════════════════════════════════════════
# Helper Functions from naph_analysis.jl
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
    
    return path_tiff
end

# ═══════════════════════════════════════════════════════════════════════════════
# Non-Interactive Analysis Function
# ═══════════════════════════════════════════════════════════════════════════════

function analyze_dataset(config, output_dir)
    """
    Run analysis without interactive plots
    Returns true if successful, false otherwise
    """
    try
        println("\n" * "="^60)
        println("Processing: $(config["week"]) / $(config["day"]) / $(config["filter"])")
        println("="^60)
        
        # Create output directory
        mkpath(output_dir)
        
        # Build data path
        path_tiff = build_data_path(config)
        
        # Check if path exists
        if !isdir(path_tiff)
            @warn "Path does not exist: $path_tiff"
            return false
        end
        
        println("Data path: $path_tiff")
        
        # Get filter information
        flt, fltv = find_filter(path_tiff, ActiveFilters)
        
        # Load image stack
        println("Loading image stack...")
        imxt = load_tif_stack_int16(path_tiff, pedestal=0.0)
        println("Stack shape: $(size(imxt))")
        
        # Get global statistics
        totalI, meanI, stdI = get_stats(imxt; bck=0.0)
        
        # Compute background
        println("Computing background...")
        background = compute_background_from_stack(imxt; 
                                                   σ=config["σ_background"], 
                                                   nlf=config["nlf_background"])
        
        # Subtract background
        println("Subtracting background...")
        imgbsub = subtract_background_from_stack(imxt, background)
        
        # Denoise stack
        println("Denoising stack...")
        imgd = denoise_stack(imgbsub; σ=config["σ_denoise"])
        
        # Detect peaks
        nframe = config["analysis_frame"]
        println("Detecting peaks in frame $nframe...")
        peaks = detect_local_maxima(imgd[:, :, nframe]; 
                                    threshold=config["peak_threshold"], 
                                    dx=config["peak_dx"], 
                                    dy=config["peak_dy"])
        
        println("Found $(size(peaks, 1)) molecule candidates")
        
        # Save peaks
        if size(peaks, 1) > 0
            CSV.write(joinpath(output_dir, "peaks.csv"), peaks)
            println("Saved peaks to CSV")
        else
            println("No peaks found, skipping peak CSV")
        end
        
        # Build traces
        println("Building traces...")
        TRZS = build_sparse_traces(imgd, peaks)
        
        # Initialize results variables
        pdf = DataFrame()
        dfe = DataFrame()
        
        # PELT fitting
        if config["fit_PELT"] && size(peaks, 1) > 0
            println("Performing PELT analysis...")
            try
                pdf, pI, pJ, pDX, pFX, pSC = find_fit_pelt(TRZS, peaks)
                println("PELT fitted $(count_mol(pdf)) molecules")
                CSV.write(joinpath(output_dir, "pelt_fits.csv"), pdf)
                println("Saved PELT results to CSV")
            catch e
                @warn "PELT fitting failed: $e"
            end
        end
        
        # Exponential fitting
        if config["fit_EXP"] && size(peaks, 1) > 0
            println("Performing exponential fitting...")
            try
                dfe, VI, VJ, VDX, VFX = fit_traces_exponential(TRZS, peaks)
                println("Exponential fitted $(size(dfe, 1)) traces")
                CSV.write(joinpath(output_dir, "exp_fits.csv"), dfe)
                println("Saved exponential results to CSV")
            catch e
                @warn "Exponential fitting failed: $e"
            end
        end
        
        # Save summary
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
        
        if !isempty(pdf)
            summary["n_molecules_pelt"] = count_mol(pdf)
        end
        
        if !isempty(dfe)
            summary["n_traces_exp"] = size(dfe, 1)
        end
        
        # Write summary
        open(joinpath(output_dir, "summary.txt"), "w") do io
            println(io, "NAPH Analysis Summary")
            println(io, "=" ^ 50)
            println(io, "Configuration:")
            println(io, "  Directory: $(config["week"])")
            println(io, "  Sample: $(config["day"])")
            println(io, "  Filter: $(config["filter"])")
            println(io, "=" ^ 50)
            for (key, value) in summary
                println(io, "$key: $value")
            end
        end
        
        println("Analysis complete for this configuration")
        return true
        
    catch e
        @error "Analysis failed: $e"
        println(stacktrace())
        return false
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Main Batch Processing
# ═══════════════════════════════════════════════════════════════════════════════

function main()
    println("\n" * "="^80)
    println("NAPH BATCH ANALYSIS")
    println("="^80)
    
    # Create base output directory
    base_output = "/Users/jjgomezcadenas/Projects/BoldLab/naph_analysis_batch"
    mkpath(base_output)
    
    # Track statistics
    total_configs = 0
    successful_runs = 0
    failed_runs = []
    
    # Process each directory
    for dir in dirs
        println("\n" * "="^80)
        println("Processing directory: $dir")
        println("="^80)
        
        # Create directory-specific output folder
        dir_output = joinpath(base_output, dir)
        mkpath(dir_output)
        
        # Select appropriate sample list
        samples = dir == "Evap_test_Fakes" ? fake_samples : signal_samples
        
        # Process each sample
        for sample in samples
            println("\n" * "-"^60)
            println("Processing sample: $sample")
            println("-"^60)
            
            # Create sample-specific output folder
            sample_output = joinpath(dir_output, sample)
            mkpath(sample_output)
            
            # Process each filter
            for filter in filters_to_process
                total_configs += 1
                
                # Create filter-specific output folder
                filter_output = joinpath(sample_output, filter)
                
                # Configure analysis
                config = copy(BASE_CONFIG)
                config["week"] = dir
                config["day"] = sample
                config["filter"] = filter
                
                # Run analysis
                success = analyze_dataset(config, filter_output)
                
                if success
                    successful_runs += 1
                    println("✓ Successfully processed: $dir/$sample/$filter")
                else
                    push!(failed_runs, "$dir/$sample/$filter")
                    println("✗ Failed to process: $dir/$sample/$filter")
                end
            end
        end
    end
    
    # Print summary
    println("\n" * "="^80)
    println("BATCH PROCESSING COMPLETE")
    println("="^80)
    println("Total configurations: $total_configs")
    println("Successful runs: $successful_runs")
    println("Failed runs: $(length(failed_runs))")
    
    if length(failed_runs) > 0
        println("\nFailed configurations:")
        for config in failed_runs
            println("  - $config")
        end
    end
    
    println("\nOutput directory: $base_output")
    println("="^80)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Script execution
# ═══════════════════════════════════════════════════════════════════════════════

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end