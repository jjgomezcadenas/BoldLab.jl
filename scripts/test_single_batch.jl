#!/usr/bin/env julia
"""
Test script to run a single configuration from batch analysis
"""

include("naph_analysis_batch.jl")

# Test with a single configuration
config = copy(BASE_CONFIG)
config["week"] = "Evap_test_Fakes"
config["day"] = "FakeA3_0.1monolayer"
config["filter"] = "Filter3"

output_dir = "/Users/jjgomezcadenas/Projects/BoldLab/test_output"
mkpath(output_dir)

println("Testing single configuration:")
println("  Week: $(config["week"])")
println("  Day: $(config["day"])")
println("  Filter: $(config["filter"])")

success = analyze_dataset(config, output_dir)

if success
    println("✓ Test successful!")
else
    println("✗ Test failed!")
end