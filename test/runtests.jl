using BoldLab
using Test
using DataFrames
using SparseArrays
using Random


@testset "Analysis Functions" begin
    
    @testset "detect_local_maxima" begin
        
        # Test case 1: Simple 3x3 matrix with one peak
        @testset "Single peak detection" begin
            frame = [1.0 2.0 1.0; 2.0 5.0 2.0; 1.0 2.0 1.0]
            peaks = detect_local_maxima(frame; threshold=3.0)
            
            @test size(peaks, 1) == 1  # Should find 1 peak
            @test peaks[1, :i] == 2   # Peak at row 2
            @test peaks[1, :j] == 2   # Peak at column 2
            @test peaks[1, :intensity] == 5.0  # Peak intensity
        end
        
        # Test case 2: Multiple peaks
        @testset "Multiple peaks detection" begin
            frame = [5.0 1.0 6.0; 1.0 1.0 1.0; 7.0 1.0 4.0]
            peaks = detect_local_maxima(frame; threshold=3.0)
            
            @test size(peaks, 1) == 4  # Should find 4 peaks (all corners are local maxima)
            peak_intensities = sort(peaks.intensity)
            @test peak_intensities == [4.0, 5.0, 6.0, 7.0]  # All four corner values
        end
        
        # Test case 3: No peaks above threshold
        @testset "No peaks above threshold" begin
            frame = [1.0 2.0 1.0; 2.0 3.0 2.0; 1.0 2.0 1.0]
            peaks = detect_local_maxima(frame; threshold=10.0)
            
            @test size(peaks, 1) == 0  # Should find no peaks
        end
        
        # Test case 4: Edge exclusion
        @testset "Edge exclusion" begin
            frame = [5.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 6.0]
            peaks = detect_local_maxima(frame; threshold=3.0, dx=0, dy=0)
            
            @test size(peaks, 1) == 2  # Should find both corner peaks without exclusion
            
            # With edge exclusion
            peaks_excluded = detect_local_maxima(frame; threshold=3.0, dx=1, dy=1)
            @test size(peaks_excluded, 1) == 0  # Should find no peaks with exclusion
        end
        
        # Test case 5: Custom window size
        @testset "Custom window size" begin
            frame = [1.0 1.0 1.0 1.0 1.0;
                    1.0 2.0 3.0 2.0 1.0;
                    1.0 3.0 9.0 3.0 1.0;
                    1.0 2.0 3.0 2.0 1.0;
                    1.0 1.0 1.0 1.0 1.0]
            
            # Test with default 3x3 window
            peaks_3x3 = detect_local_maxima(frame; threshold=5.0, window_size=(3,3))
            @test size(peaks_3x3, 1) == 1
            @test peaks_3x3[1, :intensity] == 9.0
            
            # Test with 5x5 window
            peaks_5x5 = detect_local_maxima(frame; threshold=5.0, window_size=(5,5))
            @test size(peaks_5x5, 1) == 1
            @test peaks_5x5[1, :intensity] == 9.0
        end
        
        # Test case 6: Input validation
        @testset "Input validation" begin
            frame = [1.0 2.0 1.0; 2.0 5.0 2.0; 1.0 2.0 1.0]
            
            # Should throw error for even window size
            @test_throws ArgumentError detect_local_maxima(frame; window_size=(2,2))
            @test_throws ArgumentError detect_local_maxima(frame; window_size=(3,4))
            
            # Should work fine with odd window sizes
            @test_nowarn detect_local_maxima(frame; window_size=(3,3))
            @test_nowarn detect_local_maxima(frame; window_size=(5,5))
        end
        
        # Test case 7: Return type and columns
        @testset "Return type validation" begin
            frame = [1.0 2.0 1.0; 2.0 5.0 2.0; 1.0 2.0 1.0]
            peaks = detect_local_maxima(frame; threshold=0.0)
            
            @test isa(peaks, DataFrame)
            @test names(peaks) == ["i", "j", "intensity"]
            @test eltype(peaks.i) == Int
            @test eltype(peaks.j) == Int
            @test eltype(peaks.intensity) == Float64
        end
        
        # Test case 8: Larger realistic test
        @testset "Larger matrix test" begin
            # Create 10x10 matrix with known peaks
            frame = ones(10, 10)
            frame[3, 3] = 10.0  # Peak 1
            frame[7, 8] = 8.0   # Peak 2
            frame[5, 5] = 12.0  # Peak 3 (highest)
            
            peaks = detect_local_maxima(frame; threshold=5.0)
            
            @test size(peaks, 1) == 3  # Should find 3 peaks
            
            # Check that highest peak is found
            max_peak_idx = argmax(peaks.intensity)
            @test peaks[max_peak_idx, :i] == 5
            @test peaks[max_peak_idx, :j] == 5
            @test peaks[max_peak_idx, :intensity] == 12.0
        end
    end
    
    @testset "fit_traces_exponential" begin
        
        # Create test data - sparse matrix with exponential decay traces
        @testset "Basic exponential fitting" begin
            # Generate synthetic exponential decay data
            t = 1:20
            true_A0, true_tau, true_C0 = 100.0, 5.0, 10.0
            y_true = true_A0 .* exp.(-t ./ true_tau) .+ true_C0
            
            # Add small amount of noise
            y_noisy = y_true .+ 0.1 .* randn(length(t))
            
            # Create sparse matrix and DataFrame
            TRZS = spzeros(Vector{Float64}, 2, 2)
            TRZS[1, 1] = y_noisy
            TRZS[2, 2] = y_noisy .+ 5.0  # Second trace with different offset
            
            peaks_df = DataFrame(i = [1, 2], j = [1, 2])
            
            # Fit the traces
            results, i_vals, j_vals, orig_data, fit_data = fit_traces_exponential(TRZS, peaks_df)
            
            # Check basic structure
            @test size(results, 1) == 2  # Should fit both traces
            @test length(i_vals) == 2
            @test length(j_vals) == 2
            @test length(orig_data) == 2
            @test length(fit_data) == 2
            
            # Check first trace parameters are reasonable
            fit1 = results[1, :]
            @test fit1.A0 > 80.0  # Should be close to true_A0 = 100
            @test fit1.tau > 3.0 && fit1.tau < 8.0  # Should be close to true_tau = 5
            @test fit1.C0 > 5.0   # Should be close to true_C0 = 10
            @test fit1.converged == true
            
            # Check that fitted data has same length as original
            @test length(fit_data[1]) == length(orig_data[1])
        end
        
        # Test with custom initial parameters
        @testset "Custom initial parameters" begin
            t = 1:10
            y = [50.0, 40.0, 32.0, 25.0, 20.0, 16.0, 13.0, 11.0, 9.0, 8.0]
            
            TRZS = spzeros(Vector{Float64}, 1, 1)
            TRZS[1, 1] = y
            peaks_df = DataFrame(i = [1], j = [1])
            
            # Fit with custom initial parameters
            results, _, _, _, _ = fit_traces_exponential(
                TRZS, peaks_df; initial_params=(45.0, 3.0, 5.0)
            )
            
            @test size(results, 1) == 1
            @test results[1, :converged] == true
            @test results[1, :A0] > 30.0  # Should find reasonable amplitude
        end
        
        # Test error handling
        @testset "Error handling" begin
            # Empty trace should be skipped
            TRZS = spzeros(Vector{Float64}, 2, 2) 
            TRZS[1, 1] = Float64[]  # Empty trace
            TRZS[2, 2] = [10.0, 8.0, 6.0, 4.0, 2.0]  # Valid trace
            
            peaks_df = DataFrame(i = [1, 2], j = [1, 2])
            
            results, i_vals, j_vals, orig_data, fit_data = fit_traces_exponential(TRZS, peaks_df)
            
            # Should only fit the valid trace
            @test size(results, 1) == 1
            @test length(orig_data) == 1
            @test results[1, :i] == 2  # Should be the second trace
            @test results[1, :j] == 2
        end
        
        # Test convergence tracking
        @testset "Convergence tracking" begin
            # Create a trace that should converge
            t = 1:15
            y = 20.0 .* exp.(-t ./ 4.0) .+ 3.0
            
            TRZS = spzeros(Vector{Float64}, 1, 1)
            TRZS[1, 1] = y
            peaks_df = DataFrame(i = [1], j = [1])
            
            results, _, _, _, _ = fit_traces_exponential(TRZS, peaks_df)
            
            @test size(results, 1) == 1
            @test results[1, :converged] == true
            @test results[1, :chi2] >= 0  # Chi-squared should be non-negative
            @test results[1, :chi2_red] >= 0  # Reduced chi-squared should be non-negative
        end
        
        # Test DataFrame structure
        @testset "DataFrame structure validation" begin
            t = 1:8
            y = [15.0, 12.0, 9.0, 7.0, 5.5, 4.5, 3.5, 3.0]
            
            TRZS = spzeros(Vector{Float64}, 1, 1)
            TRZS[1, 1] = y
            peaks_df = DataFrame(i = [1], j = [1])
            
            results, i_vals, j_vals, orig_data, fit_data = fit_traces_exponential(TRZS, peaks_df)
            
            # Check DataFrame columns
            expected_cols = [:i, :j, :A0, :tau, :C0, :chi2, :chi2_red, :n, :converged]
            @test names(results) == string.(expected_cols)
            
            # Check data types
            @test eltype(results.i) == Int
            @test eltype(results.j) == Int
            @test eltype(results.A0) == Float64
            @test eltype(results.tau) == Float64
            @test eltype(results.C0) == Float64
            @test eltype(results.chi2) == Float64
            @test eltype(results.chi2_red) == Float64
            @test eltype(results.n) == Int
            @test eltype(results.converged) == Bool
            
            # Check that n matches the trace length
            @test results[1, :n] == length(y)
        end
        
        # Test multiple traces
        @testset "Multiple traces fitting" begin
            t = 1:12
            
            # Create three different exponential traces
            y1 = 50.0 .* exp.(-t ./ 3.0) .+ 5.0
            y2 = 30.0 .* exp.(-t ./ 6.0) .+ 2.0 
            y3 = 80.0 .* exp.(-t ./ 2.0) .+ 8.0
            
            TRZS = spzeros(Vector{Float64}, 2, 2)
            TRZS[1, 1] = y1
            TRZS[1, 2] = y2
            TRZS[2, 1] = y3
            
            peaks_df = DataFrame(i = [1, 1, 2], j = [1, 2, 1])
            
            results, i_vals, j_vals, orig_data, fit_data = fit_traces_exponential(TRZS, peaks_df)
            
            @test size(results, 1) == 3  # Should fit all three traces
            @test all(results.converged)  # All should converge
            
            # Check that we have different parameter values (within reason)
            A0_values = results.A0
            @test length(unique(A0_values)) >= 2  # Should have some variation
            
            # Check coordinate tracking
            @test Set(zip(i_vals, j_vals)) == Set([(1,1), (1,2), (2,1)])
        end
        
        # Test max_iterations parameter
        @testset "Max iterations parameter" begin
            t = 1:5
            y = [10.0, 8.0, 6.0, 4.0, 2.0]
            
            TRZS = spzeros(Vector{Float64}, 1, 1)
            TRZS[1, 1] = y
            peaks_df = DataFrame(i = [1], j = [1])
            
            # Test with very low max_iterations (may or may not converge)
            results_low, _, _, _, _ = fit_traces_exponential(
                TRZS, peaks_df; max_iterations=5
            )
            
            # Test with high max_iterations (should converge)
            results_high, _, _, _, _ = fit_traces_exponential(
                TRZS, peaks_df; max_iterations=5000
            )
            
            # The high iteration case should be more likely to converge
            @test size(results_high, 1) >= size(results_low, 1)
        end
    end
end
