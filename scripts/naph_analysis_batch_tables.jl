#!/usr/bin/env julia
"""
NAPH Batch Analysis Tables Script
Reads summary.txt files from batch analysis and creates organized tables
Produces separate table blocks for Fakes and Signal samples
"""

using DataFrames
using CSV
using Printf

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

# Base directories
const BATCH_DIR = "/Users/jjgomezcadenas/Projects/BoldLab/naph_analysis_batch"
const OUTPUT_DIR = "/Users/jjgomezcadenas/Projects/BoldLab/naph_analysis_summary"

# Sample definitions
const FAKE_SAMPLES = [
    "FakeA1_0.5monolayer",
    "FakeA1_pre",
    "FakeA2_0.25monolayer",
    "FakeA2_pre",
    "FakeA3_0.1monolayer",
    "FakeA3_pre"
]

const SIGNAL_SAMPLES = [
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
const FILTERS = ["Filter3", "Filter4", "Filter5", "Filter6", "Filter7", "Filter8"]

# ═══════════════════════════════════════════════════════════════════════════════
# Helper Functions
# ═══════════════════════════════════════════════════════════════════════════════

function parse_summary_file(filepath)
    """
    Parse a summary.txt file and extract key metrics
    Returns a Dict with the metrics or nothing if file doesn't exist
    """
    if !isfile(filepath)
        return nothing
    end
    
    metrics = Dict{String, Any}()
    
    try
        open(filepath, "r") do io
            for line in eachline(io)
                if contains(line, "total_intensity:")
                    value_str = strip(split(line, ":")[2])
                    metrics["total_intensity"] = parse(Float64, value_str)
                elseif contains(line, "n_peaks_detected:")
                    value_str = strip(split(line, ":")[2])
                    metrics["n_peaks_detected"] = parse(Int, value_str)
                elseif contains(line, "n_molecules_pelt:")
                    value_str = strip(split(line, ":")[2])
                    metrics["n_molecules_pelt"] = parse(Int, value_str)
                end
            end
        end
        
        # Set default values if not found
        if !haskey(metrics, "n_molecules_pelt")
            metrics["n_molecules_pelt"] = 0
        end
        
        return metrics
    catch e
        @warn "Error parsing file $filepath: $e"
        return nothing
    end
end

function create_table_for_filter(dir_type, samples, filter)
    """
    Create a DataFrame table for a specific filter and sample set
    """
    # Initialize DataFrame
    df = DataFrame(
        Sample = String[],
        Total_Intensity = Float64[],
        N_Peaks = Int[],
        N_Molecules_PELT = Int[]
    )
    
    # Process each sample
    for sample in samples
        # Construct path to summary file
        summary_path = joinpath(BATCH_DIR, dir_type, sample, filter, "summary.txt")
        
        # Parse summary file
        metrics = parse_summary_file(summary_path)
        
        if metrics !== nothing
            push!(df, (
                sample,
                metrics["total_intensity"],
                metrics["n_peaks_detected"],
                metrics["n_molecules_pelt"]
            ))
        else
            # Add row with missing data
            push!(df, (sample, NaN, -1, -1))
        end
    end
    
    return df
end

function format_intensity(val)
    """Format intensity values in scientific notation"""
    if isnan(val)
        return "N/A"
    else
        return @sprintf("%.2e", val)
    end
end

function format_count(val)
    """Format count values"""
    if val < 0
        return "N/A"
    else
        return string(val)
    end
end

function save_tables_to_file(io, title, dir_type, samples)
    """
    Save formatted tables for a sample set to an IO stream
    """
    println(io, "\n" * "="^80)
    println(io, title)
    println(io, "="^80)
    
    for filter in FILTERS
        println(io, "\n" * "-"^60)
        println(io, filter)
        println(io, "-"^60)
        
        # Create table for this filter
        df = create_table_for_filter(dir_type, samples, filter)
        
        # Format and print table
        println(io, @sprintf("%-25s %15s %12s %15s", 
                            "Sample", "Total_Intensity", "N_Peaks", "N_Molecules_PELT"))
        println(io, "-"^70)
        
        for row in eachrow(df)
            println(io, @sprintf("%-25s %15s %12s %15s",
                                row.Sample,
                                format_intensity(row.Total_Intensity),
                                format_count(row.N_Peaks),
                                format_count(row.N_Molecules_PELT)))
        end
        
        # Also save as CSV for easy loading
        csv_path = joinpath(OUTPUT_DIR, "$(dir_type)_$(filter).csv")
        CSV.write(csv_path, df)
    end
end

function create_summary_tables()
    """
    Create comprehensive summary tables combining all data
    """
    println("\nCreating comprehensive summary tables...")
    
    # Create summary DataFrames for all filters
    for filter in FILTERS
        # Fakes summary
        fakes_df = create_table_for_filter("Evap_test_Fakes", FAKE_SAMPLES, filter)
        fakes_df[!, :Type] .= "Fake"
        
        # Signal summary
        signal_df = create_table_for_filter("Evap_test_NAPH", SIGNAL_SAMPLES, filter)
        signal_df[!, :Type] .= "Signal"
        
        # Combine and save
        combined_df = vcat(fakes_df, signal_df)
        csv_path = joinpath(OUTPUT_DIR, "combined_$(filter).csv")
        CSV.write(csv_path, combined_df)
        
        println("  Saved combined table for $filter")
    end
    
    # Create master summary with all filters
    master_df = DataFrame()
    
    for filter in FILTERS
        for (dir_type, samples, sample_type) in [
            ("Evap_test_Fakes", FAKE_SAMPLES, "Fake"),
            ("Evap_test_NAPH", SIGNAL_SAMPLES, "Signal")
        ]
            for sample in samples
                summary_path = joinpath(BATCH_DIR, dir_type, sample, filter, "summary.txt")
                metrics = parse_summary_file(summary_path)
                
                if metrics !== nothing
                    row_data = DataFrame(
                        Type = [sample_type],
                        Sample = [sample],
                        Filter = [filter],
                        Total_Intensity = [metrics["total_intensity"]],
                        N_Peaks = [metrics["n_peaks_detected"]],
                        N_Molecules_PELT = [metrics["n_molecules_pelt"]]
                    )
                    master_df = vcat(master_df, row_data, cols=:union)
                end
            end
        end
    end
    
    # Save master summary
    if !isempty(master_df)
        csv_path = joinpath(OUTPUT_DIR, "master_summary.csv")
        CSV.write(csv_path, master_df)
        println("  Saved master summary table")
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Main Function
# ═══════════════════════════════════════════════════════════════════════════════

function main()
    println("\n" * "="^80)
    println("NAPH BATCH ANALYSIS TABLES")
    println("="^80)
    
    # Create output directory
    mkpath(OUTPUT_DIR)
    println("\nOutput directory: $OUTPUT_DIR")
    
    # Open main summary file
    summary_file = joinpath(OUTPUT_DIR, "analysis_summary.txt")
    
    open(summary_file, "w") do io
        println(io, "NAPH Batch Analysis Summary Tables")
        println(io, "Generated: $(Dates.now())")
        println(io, "="^80)
        
        # Process Fakes
        println("\nProcessing Fake samples...")
        save_tables_to_file(io, "FAKE SAMPLES", "Evap_test_Fakes", FAKE_SAMPLES)
        
        # Process Signal samples
        println("\nProcessing Signal samples...")
        save_tables_to_file(io, "SIGNAL SAMPLES", "Evap_test_NAPH", SIGNAL_SAMPLES)
        
        # Summary statistics
        println(io, "\n" * "="^80)
        println(io, "SUMMARY STATISTICS")
        println(io, "="^80)
        
        # Calculate and print overall statistics
        for filter in FILTERS
            println(io, "\n$filter Statistics:")
            
            # Fakes
            fakes_df = create_table_for_filter("Evap_test_Fakes", FAKE_SAMPLES, filter)
            valid_fakes = Base.filter(row -> !isnan(row.Total_Intensity), fakes_df)
            if nrow(valid_fakes) > 0
                avg_peaks_fake = mean(valid_fakes.N_Peaks)
                avg_molecules_fake = mean(valid_fakes.N_Molecules_PELT)
                println(io, "  Fakes - Avg Peaks: $(round(avg_peaks_fake, digits=1)), Avg PELT Molecules: $(round(avg_molecules_fake, digits=1))")
            end
            
            # Signals
            signal_df = create_table_for_filter("Evap_test_NAPH", SIGNAL_SAMPLES, filter)
            valid_signals = Base.filter(row -> !isnan(row.Total_Intensity), signal_df)
            if nrow(valid_signals) > 0
                avg_peaks_signal = mean(valid_signals.N_Peaks)
                avg_molecules_signal = mean(valid_signals.N_Molecules_PELT)
                println(io, "  Signal - Avg Peaks: $(round(avg_peaks_signal, digits=1)), Avg PELT Molecules: $(round(avg_molecules_signal, digits=1))")
            end
        end
    end
    
    println("\nSummary file saved: $summary_file")
    
    # Create additional CSV tables for easy analysis
    create_summary_tables()
    
    # Print summary to console
    println("\n" * "="^80)
    println("PROCESSING COMPLETE")
    println("="^80)
    println("Output files created in: $OUTPUT_DIR")
    println("  - analysis_summary.txt: Formatted text tables")
    println("  - *_Filter*.csv: Individual CSV tables per filter")
    println("  - combined_Filter*.csv: Combined tables with type labels")
    println("  - master_summary.csv: Master table with all data")
    println("="^80)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Script execution
# ═══════════════════════════════════════════════════════════════════════════════

# Add required packages
using Pkg
using Statistics
using Dates

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end