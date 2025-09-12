using BoldLab
using Test
using DataFrames
using SparseArrays
using Random
using Images
using FileIO
using TiffImages
using Statistics

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
    

    @testset "load_tif_stack_int16 functions" begin
        
        # Create a temporary directory for test files
        test_dir = mktempdir()
        
        # Helper function to create test TIFF files
        function create_test_tif(filepath, width=10, height=8, intensity=0.5)
            # Create a Gray{N0f16} image with specified intensity
            img = Gray{N0f16}.(fill(intensity, height, width))
            save(filepath, img)
        end
        
        @testset "Setup test files" begin
            # Create test TIFF files with different numbering
            test_files = [
                "Test_00001.tif",
                "Test_00003.tif", 
                "Test_00002.tif",
                "Image776_00001.tif",
                "Image776_00002.tif",
                "Image780_00001.tif",
                "NotTif_00001.txt",  # Non-TIFF file
                "NoUnderscore.tif"   # No underscore file
            ]
            
            intensities = [0.1, 0.3, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8]
            
            for (i, filename) in enumerate(test_files)
                if endswith(filename, ".tif")
                    filepath = joinpath(test_dir, filename)
                    create_test_tif(filepath, 10, 8, intensities[i])
                end
            end
            
            @test isfile(joinpath(test_dir, "Test_00001.tif"))
            @test isfile(joinpath(test_dir, "Image776_00001.tif"))
        end
        
        @testset "Directory-based loading" begin
            # Test loading from directory
            stack = load_tif_stack_int16(test_dir; pedestal=0)
            
            @test isa(stack, Array{Float32, 3})
            @test size(stack) == (8, 10, 6)  # height × width × frames (6 valid TIF files with underscores)
            
            # Check that files are sorted by number
            # Test_00001, Test_00002, Test_00003, Image776_00001, Image776_00002, Image780_00001
            # Should be ordered by the number after underscore: 1, 2, 3, 1, 2, 1
            # But within same prefix, should be sorted: Test_001, Test_002, Test_003, then Image776_001, Image776_002, then Image780_001
            
            # Test dimensions are correct
            @test ndims(stack) == 3
            @test size(stack, 1) == 8   # height
            @test size(stack, 2) == 10  # width
            
            # Test pedestal subtraction
            stack_ped = load_tif_stack_int16(test_dir; pedestal=1000)
            @test all(stack_ped .< stack)  # With pedestal, values should be lower
        end
        
        @testset "Filename list with directory" begin
            # Test with filtered filename list
            filenames = ["Test_00001.tif", "Test_00003.tif", "Image776_00001.tif"]
            stack = load_tif_stack_int16(filenames, test_dir; pedestal=0)
            
            @test isa(stack, Array{Float32, 3})
            @test size(stack) == (8, 10, 3)  # Should load only the 3 specified files
            
            # Test error handling for empty list
            empty_list = String[]
            @test_throws ErrorException load_tif_stack_int16(empty_list, test_dir)
            
            # Test with non-TIFF files (should be filtered out)
            mixed_files = ["Test_00001.tif", "NotTif_00001.txt", "Test_00002.tif"]
            stack_mixed = load_tif_stack_int16(mixed_files, test_dir; pedestal=0)
            @test size(stack_mixed, 3) == 2  # Should only load the 2 TIFF files
        end
        
        @testset "Full path list" begin
            # Create full paths
            full_paths = [
                joinpath(test_dir, "Test_00001.tif"),
                joinpath(test_dir, "Test_00002.tif"),
                joinpath(test_dir, "Image776_00001.tif")
            ]
            
            stack = load_tif_stack_int16(full_paths; pedestal=0)
            
            @test isa(stack, Array{Float32, 3})
            @test size(stack) == (8, 10, 3)
            
            # Test error handling for non-existent files
            bad_paths = [joinpath(test_dir, "NonExistent_00001.tif")]
            @test_throws ArgumentError load_tif_stack_int16(bad_paths)
            
            # Test with .tiff extension (should work)
            tiff_path = joinpath(test_dir, "Test_00001.tiff")
            cp(joinpath(test_dir, "Test_00001.tif"), tiff_path)
            tiff_paths = [tiff_path]
            stack_tiff = load_tif_stack_int16(tiff_paths; pedestal=0)
            @test size(stack_tiff, 3) == 1
        end
        
        @testset "Sorting and numbering" begin
            # Create files with specific numbering to test sorting
            sort_test_dir = mktempdir()
            
            # Create files in non-sequential order
            files_order = [
                ("Data_00005.tif", 0.5),
                ("Data_00001.tif", 0.1), 
                ("Data_00010.tif", 0.9),
                ("Data_00002.tif", 0.2)
            ]
            
            for (filename, intensity) in files_order
                create_test_tif(joinpath(sort_test_dir, filename), 4, 4, intensity)
            end
            
            stack = load_tif_stack_int16(sort_test_dir; pedestal=0)
            
            @test size(stack, 3) == 4  # Should load all 4 files
            
            # Check that files are sorted by number (1, 2, 5, 10)
            # The intensity should increase accordingly: 0.1, 0.2, 0.5, 0.9
            max16 = Float32(typemax(UInt16))
            
            # Calculate expected values (intensity * max16)
            expected_vals = [0.1 * max16, 0.2 * max16, 0.5 * max16, 0.9 * max16]
            
            for i in 1:4
                avg_intensity = mean(stack[:, :, i])
                @test abs(avg_intensity - expected_vals[i]) < 1000  # Allow some tolerance
            end
        end
        
        @testset "Error handling" begin
            # Test empty directory
            empty_dir = mktempdir()
            @test_throws ErrorException load_tif_stack_int16(empty_dir)
            
            # Test directory with no underscore files
            no_underscore_dir = mktempdir()
            create_test_tif(joinpath(no_underscore_dir, "NoUnderscore.tif"))
            @test_throws ErrorException load_tif_stack_int16(no_underscore_dir)
            
            # Test non-existent directory
            @test_throws Base.IOError load_tif_stack_int16("/non/existent/directory")
        end
        
        @testset "Image processing pipeline" begin
            # Test that the processing pipeline works correctly
            # Create a test image with known values
            processing_dir = mktempdir()
            
            # Create image with specific Gray values
            img = Gray{N0f16}.(fill(0.25, 6, 8))  # 25% intensity
            save(joinpath(processing_dir, "Process_00001.tif"), img)
            
            stack = load_tif_stack_int16(processing_dir; pedestal=1000)
            
            # Calculate expected value: 0.25 * typemax(UInt16) - 1000
            expected_val = Float32(0.25 * typemax(UInt16) - 1000)
            
            # Check that all pixels have approximately the expected value
            @test all(abs.(stack[:, :, 1] .- expected_val) .< 10)  # Small tolerance for floating point
            
            @test eltype(stack) == Float32
            @test size(stack) == (6, 8, 1)
        end
        
        @testset "Multiple dispatch validation" begin
            # Verify that all three dispatch methods exist and work
            
            # Method 1: Directory only
            stack1 = load_tif_stack_int16(test_dir)
            @test isa(stack1, Array{Float32, 3})
            
            # Method 2: Filenames + directory  
            files = ["Test_00001.tif", "Test_00002.tif"]
            stack2 = load_tif_stack_int16(files, test_dir)
            @test isa(stack2, Array{Float32, 3})
            @test size(stack2, 3) == 2
            
            # Method 3: Full paths
            full_paths = [joinpath(test_dir, f) for f in files]
            stack3 = load_tif_stack_int16(full_paths)
            @test isa(stack3, Array{Float32, 3})
            @test size(stack3, 3) == 2
            
            # Verify that methods 2 and 3 produce the same result
            @test stack2 ≈ stack3
        end
        
        # Cleanup
        @testset "Cleanup test files" begin
            rm(test_dir, recursive=true)
            @test !isdir(test_dir)
        end
    end

    
end