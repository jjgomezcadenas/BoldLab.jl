#!/usr/bin/env julia
"""
NAPH Batch Analysis Script (Mikel Samples version)
Processes multiple datasets and outputs results as DataFrames
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

using CSV
using DataFrames
using Statistics
using StatsBase
using Images

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
# Helper Functions
# ═══════════════════════════════════════════════════════════════════════════════

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
    totalI, meanI, nframes, npeaks = process_image(igreen)

    results[string("total_intensity_green_", field)] = totalI
    results[string("mean_intensity_green_", field)] = meanI
    results[string("n_frames_green_", field)] = nframes
    results[string("n_peaks_detected_green_", field)] = npeaks
     
    totalI, meanI, nframes, npeaks = process_image(ired)

    results[string("total_intensity_red_", field)] = totalI
    results[string("mean_intensity_red_", field)] = meanI
    results[string("n_frames_red_", field)] = nframes
    results[string("n_peaks_detected_red_", field)] = npeaks
end

function process_image(imxt)
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
    
    return sum(totalI), sum(meanI), size(imxt, 3), size(peaks, 1)
end

# ═══════════════════════════════════════════════════════════════════════════════
# DataFrame Construction Functions
# ═══════════════════════════════════════════════════════════════════════════════

function extract_results_to_dataframe(all_results)
    # Create empty vectors for DataFrame columns
    physdir_vec = String[]
    sample_vec = String[]
    color_vec = String[]
    field_vec = Int[]
    total_intensity_vec = Float64[]
    mean_intensity_vec = Float64[]
    n_frames_vec = Int[]
    n_peaks_vec = Int[]
    
    # Process each sample's results
    for result in all_results
        phys_dir = result["phys_dir"]
        sample = result["sample"]
        
        # Extract field numbers from keys
        field_numbers = Set{Int}()
        for key in keys(result)
            # Match patterns like "total_intensity_green_880"
            m = match(r"_(green|red)_(\d+)$", key)
            if m !== nothing
                push!(field_numbers, parse(Int, m.captures[2]))
            end
        end
        
        # Process each field and color combination
        for field in sort(collect(field_numbers))
            for color in ["green", "red"]
                # Extract values for this field and color
                total_key = "total_intensity_$(color)_$(field)"
                mean_key = "mean_intensity_$(color)_$(field)"
                frames_key = "n_frames_$(color)_$(field)"
                peaks_key = "n_peaks_detected_$(color)_$(field)"
                
                if haskey(result, total_key)
                    push!(physdir_vec, phys_dir)
                    push!(sample_vec, sample)
                    push!(color_vec, color)
                    push!(field_vec, field)
                    push!(total_intensity_vec, result[total_key])
                    push!(mean_intensity_vec, result[mean_key])
                    push!(n_frames_vec, result[frames_key])
                    push!(n_peaks_vec, result[peaks_key])
                end
            end
        end
    end
    
    # Create DataFrame
    df = DataFrame(
        physdir = physdir_vec,
        sample = sample_vec,
        color = color_vec,
        field = field_vec,
        total_intensity = total_intensity_vec,
        mean_intensity = mean_intensity_vec,
        n_frames = n_frames_vec,
        n_peaks = n_peaks_vec
    )
    
    return df
end

function save_configuration(config, outdir)
    config_file = joinpath(outdir, "config_13_september_2025.csv")
    
    # Convert config dict to DataFrame
    config_df = DataFrame(
        parameter = collect(keys(config)),
        value = collect(values(config))
    )
    
    # Sort by parameter name for consistency
    sort!(config_df, :parameter)
    
    # Save to CSV
    CSV.write(config_file, config_df)
    println("Configuration saved to: $config_file")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Main Batch Processing
# ═══════════════════════════════════════════════════════════════════════════════

function main()
    config = BASE_CONFIG

    println("\n" * "="^80)
    println("NAPH BATCH (MS) ANALYSIS - DataFrame Version")
    println("="^80)
    
    # Create base output directory
    mkpath(outdir)
    
    # Store all results
    all_results = []
    
    # Process all samples
    for phys_dir in phys_dirs
        samples = phys_dir == "FREE" ? free_samples : evap_samples
        
        for sample in samples
            println("\n--- Processing $phys_dir / $sample ---")
            
            # Create results dictionary for this sample
            results = Dict()
            results["phys_dir"] = phys_dir
            results["sample"] = sample
            
            # Process the sample
            process_sample(phys_dir, sample, results)
            
            # Store results
            push!(all_results, results)
        end
    end
    
    # Convert results to DataFrame
    df = extract_results_to_dataframe(all_results)
    
    # Save configuration
    save_configuration(config, outdir)
    
    # Save DataFrame
    df_file = joinpath(outdir, "df_13_september_2025.csv")
    CSV.write(df_file, df)
    println("\nResults saved to: $df_file")
    
    # Print summary statistics
    println("\n" * "="^80)
    println("SUMMARY STATISTICS")
    println("="^80)
    
    # Group by physdir and color for summary
    gdf = groupby(df, [:physdir, :color])
    summary_df = combine(gdf, 
        :total_intensity => mean => :mean_total_intensity,
        :mean_intensity => mean => :mean_mean_intensity,
        :n_peaks => mean => :mean_n_peaks,
        nrow => :n_measurements
    )
    
    println(summary_df)
    
    println("\nTotal measurements: $(nrow(df))")
    println("Unique samples: $(length(unique(df.sample)))")
    println("Unique fields: $(length(unique(df.field)))")
    
    println("\nAnalysis completed successfully!")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Script execution
# ═══════════════════════════════════════════════════════════════════════════════

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end