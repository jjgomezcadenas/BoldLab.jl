#!/usr/bin/env julia
"""
NAPH Batch Analysis Script (Mikel Samples version)
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
# The following are available from BoldLab but not used in this script:
# using StepAnalysis
# using LabStepAnalysis  
# using StepPreprocessing
# using JPELT

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

# Directories and samples

base_dir = "/Users/jjgomezcadenas/BoldLab/BOLD/Data1/ZFFL179_NAPH3_SIL_5e-7"
outdir = "/Users/jjgomezcadenas/Projects/BoldLab/naph_analysis_batch_MS/ZFFL179_NAPH3_SIL_5e-7"
  
phys_dirs = ["FREE", "EVAP"]

free_samples = [
    "A1",
    "A2",
    "A3",
    "A4",
]

evap_samples = [
    "A6",
    "A7",
    "A8",
    "A9",
    "A10",
]

# Base configuration 
const BASE_CONFIG = Dict(
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
    "fit_PELT" => false,
    "fit_ASF" => false
)

# ═══════════════════════════════════════════════════════════════════════════════
# Helper Functions from naph_analysis.jl
# ═══════════════════════════════════════════════════════════════════════════════

function list_subfolders(dir)
    [basename(e) for e in readdir(dir; join=true) if isdir(e) && !startswith(basename(e), ".")]
end


function count_mol(df; cmol="nmol")
    length(unique(df[:, cmol]))
end



# Add scan_level function that was missing
function scan_level(path)
    entries = readdir(path, join=false, sort=true)
    dirs = String[]
    files = String[]
    all_files = String[]
    
    for entry in entries
        full_path = joinpath(path, entry)
        if isdir(full_path)
            push!(dirs, entry)
        elseif isfile(full_path)
            push!(files, entry)
            push!(all_files, entry)
        end
    end
    
    return dirs, files, all_files
end


function process_sample(phys_folder, sample_folder, results)
    path_files = joinpath(base_dir, phys_folder, sample_folder)

    println("\n" * "="^80)
    println("Processing directory: $path_files")
    println("="^80)

    _, _, all_files = scan_level(path_files)
    file_numbers = sort(unique([parse(Int, match(r"Image(\d+)_", f).captures[1]) for f in all_files]))

    for field in file_numbers
        filtered_files = filter(f -> occursin("Image$(field)_", f), all_files)
        tif_files = [joinpath(path_files, ff) for ff in filtered_files]
        imxt = load_tif_stack_int16(tif_files, pedestal=0.0)
        igreen = imxt[1:399, :, :]
        ired = imxt[400:800, :, :]

        process_images(igreen, ired, field, results)
    end
end


function process_images(igreen, ired, field, results)
    totalI, meanI, nframes, npeaks = process_image(igreen, results)

    results[string("total_intensity_green_", field)] = totalI
    results[string("mean_intensity_green_", field)] = meanI
    results[string("n_frames_green_", field)] =nframes
    results[string("n_peaks_detected_green_", field)] =npeaks
     
     totalI, meanI, nframes, npeaks  = process_image(ired, results)

    results[string("total_intensity_red_", field)] = totalI
    results[string("mean_intensity_red_", field)] = meanI
    results[string("n_frames_red_", field)] =nframes
    results[string("n_peaks_detected_red_", field)] =npeaks
end


function process_image(imxt, results)
    config = BASE_CONFIG
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
    
    return sum(totalI), sum(meanI),  size(imxt, 3), size(peaks, 1)
       
end

function write_sample(io, results)
    println(io, "  sample = $(results["sample"])")
    println(io, "  total_intensity_red = $(results["total_intensity_red"])")
    println(io, "  mean_intensity_red = $(results["mean_intensity_red"])")
    println(io, "  n_frames_red = $(results["n_frames_red"])")
    println(io, "  n_peaks_detected_red = $(results["n_peaks_detected_red"])")
    println(io, "  total_intensity_green = $(results["total_intensity_green"])")
    println(io, "  mean_intensity_green = $(results["mean_intensity_green"])")
    println(io, "  n_frames_green = $(results["n_frames_green"])")
    println(io, "  n_peaks_detected_green = $(results["n_peaks_detected_green"])")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Main Batch Processing
# ═══════════════════════════════════════════════════════════════════════════════

function main()

    config = BASE_CONFIG

    println("\n" * "="^80)
    println("NAPH BATCH (MS) ANALYSIS")
    println("="^80)
    
    # Create base output directory
    mkpath(outdir)

    # output file
    fout = joinpath(outdir, "results_13_september_2025.csv")

    # Summary
    results = Dict()

   # Write summary
    open(fout, "w") do io
        println(io, "NAPH (MS) Analysis Summary")
        println(io, "=" ^ 50)
        println(io, "Configuration:")
        println(io, "  σ_background =  $(config["σ_background"])")
        println(io, "  nlf_background: $(config["nlf_background"])")
        println(io, "  peak_threshold: $(config["peak_threshold"])")
        println(io, "  σ_denoise: $(config["σ_denoise"])")
        println(io, "  analysis_frame: $(config["analysis_frame"])")
        println(io, "=" ^ 50)

    # loop over directories

        for phys_dir in phys_dirs
            if phys_dir == "FREE"
                results["phys_dir"] = "FREE"
                println(io, " phys_dir = FREE")
                for sample in free_samples   
                    results["sample"] = sample
                    process_sample(phys_dir, sample, results)
                    write_sample(io, results)    
                end
            else
                results["phys_dir"] = "EVAP"
                println(io, " phys_dir = EVAP")
                for sample in evap_samples
                    results["sample"] = sample
                    process_sample(phys_dir, sample, results)
                    write_sample(io, results)
                end
            end
        end
    end
    
    println("Analysis completed")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Script execution
# ═══════════════════════════════════════════════════════════════════════════════

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end