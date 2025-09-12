using BoldLab
using Test
using DataFrames
using SparseArrays
using Random
using Statistics

@testset "Fit Functions" begin
    
    # Helper function to generate exponential decay data with noise
    function generate_exponential_data(A0::Float64, tau::Float64, C0::Float64, length::Int; noise_level::Float64=0.02, seed::Int=42)
        Random.seed!(seed)
        t = collect(1:length)
        # Generate true exponential decay: A₀*exp(-t/τ) + C₀
        true_values = A0 .* exp.(-t ./ tau) .+ C0
        # Add Gaussian noise (2% of signal range by default)
        signal_range = maximum(true_values) - minimum(true_values)
        noise = randn(length) .* (noise_level * signal_range)
        noisy_data = true_values .+ noise
        return noisy_data, true_values
    end
    
    # Helper function to generate step trace data with noise
    function generate_step_trace(step_heights::Vector{Float64}, step_lengths::Vector{Int}; 
                                 noise_level::Float64=0.02, seed::Int=42)
        """
        Generate a trace with multiple steps of specified heights and lengths.
        
        Parameters:
        - step_heights: Vector of step heights (values for each plateau)
        - step_lengths: Vector of step durations (number of points in each step)
        - noise_level: Fraction of signal range to use as noise standard deviation
        - seed: Random seed for reproducibility
        
        Returns:
        - noisy_data: The trace with added noise
        - true_values: The clean step trace
        - changepoints: Indices where steps occur
        """
        @assert length(step_heights) == length(step_lengths) "Heights and lengths must have same size"
        
        Random.seed!(seed)
        
        # Build the clean trace
        true_values = Float64[]
        changepoints = Int[1]  # First changepoint is at the start
        
        for (height, len) in zip(step_heights, step_lengths)
            append!(true_values, fill(height, len))
            if length(true_values) < sum(step_lengths)  # Not the last step
                push!(changepoints, length(true_values) + 1)
            end
        end
        
        # Add noise
        if noise_level > 0
            signal_range = maximum(true_values) - minimum(true_values)
            if signal_range > 0
                noise = randn(length(true_values)) .* (noise_level * signal_range)
            else
                # If all steps have same height, use absolute noise
                noise = randn(length(true_values)) .* noise_level
            end
            noisy_data = true_values .+ noise
        else
            noisy_data = copy(true_values)
        end
        
        return noisy_data, true_values, changepoints
    end
    
    # Helper function to generate random step parameters
    function generate_random_steps(; n_steps::Int=3, min_length::Int=10, max_length::Int=30,
                                   min_height::Float64=10.0, max_height::Float64=100.0, 
                                   seed::Int=42)
        """
        Generate random step parameters for testing.
        
        Parameters:
        - n_steps: Number of steps to generate
        - min_length, max_length: Range for step lengths
        - min_height, max_height: Range for step heights
        - seed: Random seed
        
        Returns:
        - step_heights: Vector of random step heights
        - step_lengths: Vector of random step lengths
        """
        Random.seed!(seed)
        
        step_heights = min_height .+ (max_height - min_height) .* rand(n_steps)
        step_lengths = rand(min_length:max_length, n_steps)
        
        return step_heights, step_lengths
    end
    
    @testset "Exponential data generation" begin
        # Test the data generation function itself
        data, true_vals = generate_exponential_data(100.0, 10.0, 5.0, 50)
        @test length(data) == 50
        @test length(true_vals) == 50
        # Check that data has some noise but follows general trend
        @test cor(data, true_vals) > 0.95  # High correlation despite noise
    end
    
    @testset "Step trace generation" begin
        # Test basic step trace generation
        step_heights = [10.0, 30.0, 20.0, 40.0]
        step_lengths = [15, 20, 10, 15]
        
        noisy_trace, clean_trace, changepoints = generate_step_trace(step_heights, step_lengths)
        
        @test length(noisy_trace) == sum(step_lengths)
        @test length(clean_trace) == sum(step_lengths)
        @test length(changepoints) == length(step_heights)
        
        # Check that clean trace has correct structure
        idx = 1
        for (height, len) in zip(step_heights, step_lengths)
            segment = clean_trace[idx:idx+len-1]
            @test all(segment .== height)
            idx += len
        end
        
        # Test that noise is applied correctly
        @test cor(noisy_trace, clean_trace) > 0.95  # Should be highly correlated
        @test noisy_trace != clean_trace  # But not identical
        
        # Test random step generation
        heights, lengths = generate_random_steps(n_steps=5, min_length=10, max_length=25)
        @test length(heights) == 5
        @test length(lengths) == 5
        @test all(10 .<= lengths .<= 25)
        @test all(10.0 .<= heights .<= 100.0)
        
        # Test with no noise
        noisy_trace2, clean_trace2, _ = generate_step_trace([50.0, 100.0], [20, 30]; noise_level=0.0)
        @test noisy_trace2 == clean_trace2
    end
    
    @testset "fit_traces_exponential - single trace" begin
        
        @testset "Parameter recovery - various trace lengths" begin
            # Test with different trace lengths (all >= 10)
            trace_lengths = [10, 15, 25, 50, 100]
            
            for trace_length in trace_lengths
                # Known parameters for testing
                A0_true, tau_true, C0_true = 100.0, 15.0, 10.0
                
                # Generate noisy exponential data
                trace, _ = generate_exponential_data(A0_true, tau_true, C0_true, trace_length)
                
                # Fit the trace
                result = fit_traces_exponential(trace)
                
                @test result !== nothing  # Fit should succeed
                @test isa(result, BoldLab.ExponentialFit)
                
                # Check parameter recovery within 15% tolerance (more relaxed for short traces)
                @test abs(result.A0 - A0_true) / A0_true < 0.15
                @test abs(result.tau - tau_true) / tau_true < 0.15  
                @test abs(result.C0 - C0_true) / abs(C0_true) < 0.95  # C0 can be harder to recover
                
                # Check that chi-squared is reasonable (should be close to 1 for good fit)
                @test result.χ2_red < 5.0  # Allow some flexibility due to noise
                @test result.χ2_red > 0.1  # But not too small (would indicate overfitting)
                
                # Check that fitted data has same length as input
                @test length(result.fitX) == trace_length
                @test length(result.dataX) == trace_length
            end
        end
        
        @testset "Different parameter sets" begin
            # Test with various parameter combinations
            parameter_sets = [
                (50.0, 5.0, 0.0),    # Short tau, no offset
                (200.0, 25.0, 50.0), # Long tau, large offset  
                (10.0, 8.0, 2.0),    # Small amplitude
                (500.0, 40.0, 100.0) # Large values
            ]
            
            for (A0_true, tau_true, C0_true) in parameter_sets
                trace, _ = generate_exponential_data(A0_true, tau_true, C0_true, 30)
                result = fit_traces_exponential(trace)
                
                @test result !== nothing
                
                # Parameter recovery within 20% (can be harder with limited data)
                @test abs(result.A0 - A0_true) / A0_true < 0.20
                @test abs(result.tau - tau_true) / tau_true < 0.20
                if C0_true != 0
                    @test abs(result.C0 - C0_true) / abs(C0_true) < 0.20
                else
                    @test abs(result.C0) < 1.0  # Should be small if true value is 0
                end
                
                # Reasonable chi-squared (allow higher values for difficult fits)
                @test result.χ2_red < 30.0
            end
        end
        
        @testset "Custom initial parameters" begin
            A0_true, tau_true, C0_true = 80.0, 12.0, 15.0
            trace, _ = generate_exponential_data(A0_true, tau_true, C0_true, 25)
            
            # Test with custom initial parameters (close to true values)
            initial_params = [85.0, 10.0, 12.0]
            result = fit_traces_exponential(trace; initial_params=initial_params)
            
            @test result !== nothing
            @test abs(result.A0 - A0_true) / A0_true < 0.15
            @test abs(result.tau - tau_true) / tau_true < 0.15
            @test abs(result.C0 - C0_true) / abs(C0_true) < 0.15
        end
        
        @testset "Error handling" begin
            # Test with empty trace
            empty_trace = Float64[]
            result = fit_traces_exponential(empty_trace)
            @test result === nothing
            
            # Test with very short trace (should still work but might be less accurate)
            short_trace, _ = generate_exponential_data(50.0, 10.0, 5.0, 5)
            result = fit_traces_exponential(short_trace)
            # Should either succeed or fail gracefully
            @test (result === nothing) || isa(result, BoldLab.ExponentialFit)
        end
    end
    
    @testset "fit_traces_exponential - DataFrame/SparseMatrix" begin
        
        @testset "Multiple traces fitting" begin
            # Create a sparse matrix with multiple traces
            n_traces = 5
            trace_length = 30
            
            # Generate test data for multiple traces with different parameters
            test_params = [
                (100.0, 15.0, 10.0),
                (80.0, 12.0, 5.0), 
                (150.0, 20.0, 15.0),
                (60.0, 8.0, 3.0),
                (120.0, 18.0, 8.0)
            ]
            
            # Create sparse matrix and DataFrame
            TRZS = SparseMatrixCSC{Vector{Float64}, Int}(n_traces, n_traces, 
                                                         collect(1:n_traces+1), 
                                                         collect(1:n_traces), 
                                                         Vector{Vector{Float64}}(undef, n_traces))
            
            df = DataFrame(i=Int[], j=Int[])
            
            # Fill with test traces
            for k in 1:n_traces
                A0, tau, C0 = test_params[k]
                trace, _ = generate_exponential_data(A0, tau, C0, trace_length; seed=k*10)
                TRZS[k, k] = trace
                push!(df, (i=k, j=k))
            end
            
            # Fit all traces
            fit_results, indices_i, indices_j, original_data, fitted_data = fit_traces_exponential(TRZS, df)
            
            # Check that we got results
            @test nrow(fit_results) > 0
            @test length(indices_i) == length(indices_j)
            @test length(original_data) == length(fitted_data)
            
            # Check that parameters are recovered within tolerance for fitted traces
            for row in eachrow(fit_results)
                trace_idx = findfirst(x -> x == row.i, indices_i)
                if trace_idx !== nothing
                    A0_true, tau_true, C0_true = test_params[trace_idx]
                    
                    # Parameter recovery within 20% tolerance (slightly more lenient for batch fitting)
                    @test abs(row.A0 - A0_true) / A0_true < 0.20
                    @test abs(row.tau - tau_true) / tau_true < 0.20
                    @test abs(row.C0 - C0_true) / abs(C0_true) < 0.60  # C0 recovery can be harder
                    
                    # Reasonable chi-squared
                    @test row.chi2 < 10.0
                end
            end
            
            # Check data arrays consistency
            @test length(original_data) == nrow(fit_results)
            @test length(fitted_data) == nrow(fit_results)
            
            for k in 1:length(original_data)
                @test length(original_data[k]) == trace_length
                @test length(fitted_data[k]) == trace_length
            end
        end
        
        @testset "DataFrame structure validation" begin
            # Test with minimal dataset
            TRZS = SparseMatrixCSC{Vector{Float64}, Int}(2, 2, [1, 2, 3], [1, 2], 
                                                         Vector{Vector{Float64}}(undef, 2))
            
            trace1, _ = generate_exponential_data(50.0, 10.0, 5.0, 20)
            trace2, _ = generate_exponential_data(75.0, 15.0, 8.0, 20)
            
            TRZS[1, 1] = trace1
            TRZS[2, 2] = trace2
            
            df = DataFrame(i=[1, 2], j=[1, 2])
            
            fit_results, indices_i, indices_j, original_data, fitted_data = fit_traces_exponential(TRZS, df)
            
            # Check DataFrame column names and types
            expected_columns = ["i", "j", "A0", "tau", "C0", "chi2"]
            @test names(fit_results) == expected_columns
            @test eltype(fit_results.i) == Int
            @test eltype(fit_results.j) == Int
            @test eltype(fit_results.A0) == Float64
            @test eltype(fit_results.tau) == Float64  
            @test eltype(fit_results.C0) == Float64
            @test eltype(fit_results.chi2) == Float64
        end
        
        @testset "Failed fits handling" begin
            # Create traces that might fail to fit (e.g., constant values, very noisy)
            TRZS = SparseMatrixCSC{Vector{Float64}, Int}(3, 3, [1, 2, 3, 4], [1, 2, 3], 
                                                         Vector{Vector{Float64}}(undef, 3))
            
            # Good trace
            good_trace, _ = generate_exponential_data(100.0, 15.0, 10.0, 25)
            TRZS[1, 1] = good_trace
            
            # Constant trace (should fail or give poor fit)  
            constant_trace = fill(50.0, 25)
            TRZS[2, 2] = constant_trace
            
            # Very noisy trace (might fail)
            noisy_trace, _ = generate_exponential_data(20.0, 5.0, 2.0, 25; noise_level=0.5)
            TRZS[3, 3] = noisy_trace
            
            df = DataFrame(i=[1, 2, 3], j=[1, 2, 3])
            
            fit_results, indices_i, indices_j, original_data, fitted_data = fit_traces_exponential(TRZS, df)
            
            # Should have at least the good trace fitted
            @test nrow(fit_results) >= 1
            
            # Check that failed fits are properly excluded from results
            @test length(indices_i) == nrow(fit_results)
            @test length(original_data) == nrow(fit_results)
            @test length(fitted_data) == nrow(fit_results)
        end
    end
    
    @testset "fit_pelt - PELT step detection" begin
        
        @testset "Simple step detection" begin
            # Create a simple trace with clear steps
            step_heights = [20.0, 50.0, 30.0, 60.0]
            step_lengths = [20, 25, 15, 20]
            
            trace, clean_trace, true_changepoints = generate_step_trace(step_heights, step_lengths; 
                                                                        noise_level=0.01, seed=100)
            
            # Run PELT fit
            cro, noise_i, noise_f = pelt_fit(trace; i=1, f=5, auto_noise_region=true)
            
            @test cro !== nothing
            @test haskey(cro, "changepoints")
            @test haskey(cro, "penalty")
            @test haskey(cro, "constrained")
            @test haskey(cro, "number")
            
            # Test with first changepoint set (ic=1)
            if length(cro["changepoints"]) >= 1
                result = fit_pelt(trace, cro, 1)
                
                if result !== nothing
                    @test isa(result, BoldLab.PeltFit)
                    @test result.nstep > 0  # Should detect at least one step
                    @test length(result.stepHeight) == result.nstep
                    @test length(result.stepTime) == result.nstep
                    @test length(result.stepLength) == result.nstep
                    @test result.χ2_red < 10.0  # Should have reasonable fit
                    
                    # Check that fitted trace has same length as input
                    @test length(result.fitX) == length(trace)
                    @test length(result.dataX) == length(trace)
                end
            end
        end
        
        @testset "Multiple complexity levels" begin
            # Test with different numbers of steps
            test_cases = [
                ([30.0, 60.0], [30, 30]),           # 2 steps
                ([10.0, 40.0, 20.0], [20, 20, 20]), # 3 equal-length steps
                ([100.0, 50.0, 75.0, 25.0, 60.0], [15, 20, 15, 10, 20]) # 5 varied steps
            ]
            
            for (heights, lengths) in test_cases
                trace, _, _ = generate_step_trace(heights, lengths; noise_level=0.02)
                
                cro, noise_i, noise_f = pelt_fit(trace; auto_noise_region=true)
                
                if cro !== nothing && length(cro["changepoints"]) >= 1
                    result = fit_pelt(trace, cro, 1)
                    
                    if result !== nothing
                        # Should detect approximately the right number of steps
                        # Allow some flexibility due to noise
                        @test abs(result.nstep - (length(heights) - 1)) <= 2
                        @test result.χ2_red < 20.0
                    end
                end
            end
        end
        
        @testset "Edge cases" begin
            # Test with constant trace (no steps)
            constant_trace = fill(50.0, 50) .+ randn(50) * 0.5
            cro, noise_i, noise_f = pelt_fit(constant_trace; auto_noise_region=true)
            
            if cro !== nothing && length(cro["changepoints"]) >= 1
                result = fit_pelt(constant_trace, cro, 1)
                # Should either fail or detect very few steps
                @test (result === nothing) || (result.nstep <= 1)
            end
            
            # Test with very noisy trace
            heights = [20.0, 80.0, 40.0]
            lengths = [20, 20, 20]
            noisy_trace, _, _ = generate_step_trace(heights, lengths; noise_level=0.5)
            
            cro, noise_i, noise_f = pelt_fit(noisy_trace; auto_noise_region=true)
            if cro !== nothing && length(cro["changepoints"]) >= 1
                result = fit_pelt(noisy_trace, cro, 1)
                # Should handle gracefully
                @test (result === nothing) || isa(result, BoldLab.PeltFit)
            end
        end
        
        @testset "fitpelt optimization" begin
            # Test the fitpelt function that optimizes over multiple ic values
            step_heights = [25.0, 75.0, 50.0]
            step_lengths = [25, 30, 25]
            trace, _, _ = generate_step_trace(step_heights, step_lengths; noise_level=0.015)
            
            # The fitpelt function may return nothing if PELT fitting fails
            result = fitpelt(trace; max_ic=3, auto_noise_region=true)
            
            if result !== nothing
                icx, pfit_results = result
                @test icx >= 1  # Should find at least one changepoint set
                @test length(pfit_results) <= min(icx, 3)  # Should not exceed max_ic
                
                # Check that all results are valid PeltFit objects
                for pfit in pfit_results
                    @test isa(pfit, BoldLab.PeltFit)
                    @test pfit.χ2_red >= 0
                    @test pfit.nstep >= 0
                end
                
                # Results should be ordered by ic (implicitly by construction)
                if length(pfit_results) > 1
                    # Later fits might have more changepoints
                    @test pfit_results[end].nstep >= pfit_results[1].nstep
                end
            end
        end
    end
    
    @testset "fit_asf - Automatic Step Finder" begin
        
        @testset "Basic step finding" begin
            # Create a trace with clear steps
            step_heights = [15.0, 45.0, 25.0, 55.0]
            step_lengths = [20, 30, 20, 25]
            
            trace, clean_trace, _ = generate_step_trace(step_heights, step_lengths; 
                                                        noise_level=0.01, seed=200)
            
            # Run ASF fit
            result = fit_asf(trace; niter=10, tresH=0.15)
            
            if result !== nothing
                @test isa(result, BoldLab.AsfFit)
                @test result.best_shot > 0  # Should have positive quality measure
                # iter can be negative in ASF (indicates early stopping)
                @test result.iter != 0  # Should have performed some processing
                @test length(result.stepHeights) > 0  # Should detect steps
                @test length(result.stepTimes) == length(result.stepHeights)
                @test length(result.steps) == length(result.stepHeights)
                @test result.χ2_red < 25.0  # Should have reasonable fit (slightly relaxed)
                
                # Check data arrays
                @test length(result.fitX) == length(trace)
                @test length(result.dataX) == length(trace)
                @test length(result.S_curve) > 0
            end
        end
        
        @testset "Different step patterns" begin
            # Test various step configurations
            test_patterns = [
                ([20.0, 60.0], [40, 40]),              # Simple two-level
                ([10.0, 50.0, 30.0], [25, 25, 25]),    # Three equal steps
                ([5.0, 25.0, 15.0, 35.0, 20.0], [15, 18, 20, 17, 15])  # Complex pattern
            ]
            
            for (heights, lengths) in test_patterns
                trace, _, _ = generate_step_trace(heights, lengths; noise_level=0.02, seed=300)
                
                result = fit_asf(trace; niter=15, tresH=0.10)
                
                if result !== nothing
                    @test isa(result, BoldLab.AsfFit)
                    # Should detect reasonable number of steps
                    @test 0 < length(result.stepHeights) <= length(heights)
                    @test result.χ2_red < 30.0
                    @test result.best_shot > 0
                end
            end
        end
        
        @testset "Parameter variations" begin
            # Test with different ASF parameters
            step_heights = [30.0, 70.0, 40.0]
            step_lengths = [30, 35, 30]
            trace, _, _ = generate_step_trace(step_heights, step_lengths; noise_level=0.015)
            
            # Test with different iteration counts
            for niter in [5, 10, 20]
                result = fit_asf(trace; niter=niter, tresH=0.15)
                if result !== nothing
                    @test result.iter <= niter
                    @test length(result.stepHeights) > 0
                end
            end
            
            # Test with different thresholds
            for threshold in [0.05, 0.15, 0.25]
                result = fit_asf(trace; niter=10, tresH=threshold)
                if result !== nothing
                    @test isa(result, BoldLab.AsfFit)
                    @test result.χ2_red >= 0
                end
            end
        end
        
        @testset "Challenging traces" begin
            # Test with noisy trace
            heights = [25.0, 75.0]
            lengths = [40, 40]
            noisy_trace, _, _ = generate_step_trace(heights, lengths; noise_level=0.1, seed=400)
            
            result = fit_asf(noisy_trace; niter=20, tresH=0.20)
            # Should handle noisy data gracefully
            @test (result === nothing) || (isa(result, BoldLab.AsfFit) && result.χ2_red < 100)
            
            # Test with many small steps
            many_steps_heights = collect(range(10.0, 90.0, length=8))
            many_steps_lengths = fill(10, 8)
            complex_trace, _, _ = generate_step_trace(many_steps_heights, many_steps_lengths; 
                                                      noise_level=0.02)
            
            result = fit_asf(complex_trace; niter=15, tresH=0.10)
            if result !== nothing
                @test length(result.stepHeights) > 0
                @test result.best_shot > 0
            end
            
            # Test with gradual changes (ramp-like)
            ramp_heights = collect(range(20.0, 80.0, length=6))
            ramp_lengths = fill(15, 6)
            ramp_trace, _, _ = generate_step_trace(ramp_heights, ramp_lengths; noise_level=0.01)
            
            result = fit_asf(ramp_trace; niter=10, tresH=0.15)
            if result !== nothing
                @test isa(result, BoldLab.AsfFit)
                # May detect multiple steps in the ramp
                @test length(result.stepHeights) >= 1
            end
        end
    end
    
    @testset "elbow - Elbow/Knee Point Detection" begin
        
        @testset "Basic elbow detection" begin
            # Create a simple L-shaped curve (classic elbow)
            x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            y = [10.0, 9.0, 7.8, 6.2, 4.8, 3.8, 3.2, 2.8, 2.6, 2.5, 2.4]
            
            # Test kneedle method
            elbow_idx = elbow(x, y; method=:kneedle, return_index=true)
            elbow_point = elbow(x, y; method=:kneedle, return_index=false)
            
            @test isa(elbow_idx, Int)
            @test 1 ≤ elbow_idx ≤ length(x)
            @test isa(elbow_point, Tuple)
            @test length(elbow_point) == 2
            @test elbow_point[1] == x[elbow_idx]
            @test elbow_point[2] == y[elbow_idx]
            
            # Test curvature method
            elbow_idx_curv = elbow(x, y; method=:curvature, return_index=true)
            @test isa(elbow_idx_curv, Int)
            @test 1 ≤ elbow_idx_curv ≤ length(x)
        end
        
        @testset "Exponential decay curve" begin
            # Create exponential decay (common in fitting contexts)
            x = collect(1:20)
            y = 100.0 .* exp.(-x ./ 5.0) .+ 5.0
            
            elbow_idx = elbow(x, y; method=:kneedle, return_index=true)
            @test 2 ≤ elbow_idx ≤ length(x)-1  # Should not be at endpoints
            
            # Test both methods give reasonable results
            elbow_kn = elbow(x, y; method=:kneedle, return_index=false)
            elbow_cv = elbow(x, y; method=:curvature, return_index=false)
            
            @test elbow_kn[1] in x
            @test elbow_cv[1] in x
        end
        
        @testset "Data cleaning and preprocessing" begin
            # Test with NaN values
            x = [1.0, 2.0, NaN, 4.0, 5.0, 6.0]
            y = [10.0, 8.0, 6.0, NaN, 2.0, 1.0]
            
            # Should handle NaNs gracefully
            elbow_point = elbow(x, y; method=:kneedle, return_index=false)
            @test !any(isnan.([elbow_point[1], elbow_point[2]]))
            
            # Test with duplicate x values
            x_dup = [1.0, 2.0, 2.0, 3.0, 4.0, 5.0]
            y_dup = [10.0, 8.0, 7.0, 5.0, 3.0, 1.0]
            
            elbow_idx = elbow(x_dup, y_dup; method=:kneedle, return_index=true)
            @test isa(elbow_idx, Int)
            
            # Test with unsorted x values
            x_unsorted = [3.0, 1.0, 5.0, 2.0, 4.0]
            y_unsorted = [5.0, 10.0, 1.0, 8.0, 3.0]
            
            elbow_idx = elbow(x_unsorted, y_unsorted; method=:kneedle, return_index=true)
            @test isa(elbow_idx, Int)
        end
        
        @testset "Smoothing parameters" begin
            # Create noisy curve
            x = collect(1:20)
            y_clean = 50.0 .* exp.(-x ./ 8.0) .+ 2.0
            y_noisy = y_clean .+ 0.5 .* randn(length(x))
            
            # Test different smoothing windows (curvature method only)
            for smooth_window in [3, 5, 7]
                elbow_idx = elbow(x, y_noisy; method=:curvature, smooth_window=smooth_window, return_index=true)
                @test isa(elbow_idx, Int)
            end
            
            # Test multiple smoothing passes
            elbow_idx = elbow(x, y_noisy; method=:curvature, smooth_passes=3, return_index=true)
            @test isa(elbow_idx, Int)
        end
        
        @testset "Error handling" begin
            # Test insufficient points
            x_short = [1.0, 2.0]
            y_short = [10.0, 5.0]
            @test_throws AssertionError elbow(x_short, y_short)
            
            # Test mismatched lengths
            x = [1.0, 2.0, 3.0]
            y = [10.0, 5.0]
            @test_throws AssertionError elbow(x, y)
            
            # Test unknown method
            x = [1.0, 2.0, 3.0, 4.0]
            y = [10.0, 8.0, 5.0, 2.0]
            @test_throws ErrorException elbow(x, y; method=:unknown)
            
            # Test even smooth_window
            @test_throws AssertionError elbow(x, y; method=:curvature, smooth_window=4)
        end
        
        @testset "Edge cases" begin
            # Test with constant y values (no elbow)
            x = collect(1:10)
            y = fill(5.0, 10)
            
            # Should still return an index
            elbow_idx = elbow(x, y; method=:kneedle, return_index=true)
            @test isa(elbow_idx, Int)
            
            # Test with linear relationship
            y_linear = x
            elbow_idx = elbow(x, y_linear; method=:kneedle, return_index=true)
            @test isa(elbow_idx, Int)
            
            # Test minimal case (3 points)
            x_min = [1.0, 2.0, 3.0]
            y_min = [10.0, 5.0, 1.0]
            elbow_point = elbow(x_min, y_min; method=:kneedle, return_index=false)
            @test elbow_point[1] == 2.0  # Middle point should be the elbow
        end
    end
    
    @testset "movavg - Moving Average Helper Function" begin
        
        @testset "Basic moving average" begin
            # Test with simple sequence
            v = [1.0, 2.0, 3.0, 4.0, 5.0]
            
            # Window size 3
            result = BoldLab.movavg(v, 3)
            @test length(result) == length(v)
            @test result[3] ≈ 3.0  # Middle should be exact average
            
            # Check edge handling (should replicate edges)
            @test result[1] ≈ (1.0 + 2.0) / 2  # First element uses [1,2]
            @test result[end] ≈ (4.0 + 5.0) / 2  # Last element uses [4,5]
        end
        
        @testset "Different window sizes" begin
            v = collect(1.0:10.0)
            
            # Test odd window sizes
            for w in [1, 3, 5, 7]
                result = BoldLab.movavg(v, w)
                @test length(result) == length(v)
                @test eltype(result) == Float64
            end
            
            # Window size 1 should return original (as Float64)
            result1 = BoldLab.movavg(v, 1)
            @test result1 ≈ Float64.(v)
            
            # Large window size
            result_large = BoldLab.movavg(v, 15)  # Larger than vector
            # Should still work (all points get close to mean due to edge replication)
            @test length(result_large) == length(v)
            # With large window, middle values should be close to the mean
            @test abs(result_large[5] - mean(v)) < 1.0
        end
        
        @testset "Smoothing noisy data" begin
            # Create noisy signal
            t = collect(1:100)
            clean_signal = sin.(t ./ 10)
            noisy_signal = clean_signal .+ 0.1 .* randn(100)
            
            # Apply smoothing
            smoothed = BoldLab.movavg(noisy_signal, 5)
            
            # Smoothed should be closer to clean signal
            error_original = sum(abs2, noisy_signal - clean_signal)
            error_smoothed = sum(abs2, smoothed - clean_signal)
            @test error_smoothed < error_original
        end
        
        @testset "Type handling" begin
            # Test with different input types
            v_int = [1, 2, 3, 4, 5]
            result_int = BoldLab.movavg(v_int, 3)
            @test eltype(result_int) == Float64
            
            # Test with Float32
            v_f32 = Float32[1, 2, 3, 4, 5]
            result_f32 = BoldLab.movavg(v_f32, 3)
            @test eltype(result_f32) == Float64
        end
        
        @testset "Edge cases" begin
            # Single element
            v_single = [5.0]
            result_single = BoldLab.movavg(v_single, 3)
            @test result_single == [5.0]
            
            # Two elements
            v_two = [1.0, 3.0]
            result_two = BoldLab.movavg(v_two, 3)
            @test result_two ≈ [2.0, 2.0]  # Both should be average
            
            # Empty vector
            v_empty = Float64[]
            result_empty = BoldLab.movavg(v_empty, 3)
            @test isempty(result_empty)
        end
    end
    
end