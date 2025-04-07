### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 9292025e-0d7d-11f0-365e-f1724fc39b3c
using Pkg; Pkg.activate(ENV["JBoldLab"])

# ╔═╡ f8573e77-38d8-467c-a123-9a83f73e8970
push!(LOAD_PATH, ENV["JBoldLab"] * "/src")

# ╔═╡ 19476521-78bd-44d6-a302-6e971f9fc095
begin
	using Revise
	using BoldLab
	using SimpleLogger
	using AutoStepfinder
end

# ╔═╡ 5b3b5f83-058a-4582-9482-9b2fc5892692
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Plots 
	using Measures
	using Printf
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using DelimitedFiles
	using Images, FileIO
end

# ╔═╡ 583a9aee-08eb-4f5a-95ef-d0087eb98cbc
names(SimpleLogger)

# ╔═╡ 39b011c6-f511-42dd-befc-eaf3fd17ea1a
names(AutoStepfinder)

# ╔═╡ 8e8cc911-5a3c-4ec5-9498-cd60a09fd276
set_log_level(SimpleLogger.DEBUG)

# ╔═╡ b7508873-d85d-4502-8845-b779feb4666c
debug("test")

# ╔═╡ 88b1aca9-3390-40c2-a9c6-cdce13528b07
begin
	rdir = "/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign/24_Monday"
	sample = "Dirty_coverslip"
	field = "Field1"
end

# ╔═╡ 89114d7b-298d-41cb-bc48-bdf32865a735
function getfiles(rdir, sample, field; pixelsize=512)

	function convert_to_gray(img)
	# Convert to grayscale if needed
    	if ndims(img) == 2
        	img_array = img
    	elseif ndims(img) == 3
        	# Use first channel if image has multiple channels
        	img_array = channelview(img)[1, :, :]
    	elseif ndims(img) == 1
       		error("Image appears to be a 1D array. Unexpected format in file: $file")
    	else
        	error("Unhandled image dimension: $(ndims(img)) in file: $file")
    	end
		img_array
	end
	
	folder_path = joinpath(rdir, sample, field)
	
	# Get a list of all .tif files, sorted by name (assumes filenames are ordered correctly)
	image_files = filter(f -> occursin(r"Image1_\d+\.tif", f), 
						 readdir(folder_path))


	# Extract the number from each filename for sorting
	sorted_files = sort(image_files, 
						by = fname -> parse(Int, 
											match(r"Image1_(\d+)\.tif", 
												  fname).captures[1]))


	# Prepend folder path to filenames
	sorted_paths = joinpath.(folder_path, sorted_files)

	num_images = length(sorted_paths)
	image_stack = Array{Float32}(undef, pixelsize, pixelsize, num_images)
	
	# Load each image into the 3D array
	for (i, file) in enumerate(sorted_paths)
	    img = load(file)
		debug("Loading image $i: $(basename(file))")
		debug("Raw image size: $(size(img))")
		debug("Type:  $(typeof(img))")
		img_array = convert_to_gray(img)

		# Convert and assign
    	image_stack[:, :, i] = Float32.(img_array) 
	end

	println("Loaded image stack with size: ", size(image_stack))  # (512, 512, 400)
	image_stack
end

# ╔═╡ 7009519c-c468-4547-9d26-32e54e90404f
"""
Computes the gloabl statistics of the stack. Sum, mean and stds. Also creates
a vector of means and stds.
imst is the stack  
"""
function get_stats(imst::Array{Float32, 3})
	
	n_frames = size(imst, 3)
	
	# Preallocate vectors
	total_intensity = zeros(Float32, n_frames)
	mean_intensity = zeros(Float32, n_frames)
	std_intensity = zeros(Float32, n_frames)
	
	# Loop over each image in the stack
	for i in 1:n_frames
	    frame = imst[:, :, i]
	    total_intensity[i] = sum(frame)
	    mean_intensity[i] = mean(frame)
	    std_intensity[i] = std(frame)
	end
	total_intensity, mean_intensity, std_intensity 
end
	


# ╔═╡ fe61875d-8e06-48b0-b67b-5ed01b6dcc6a
function plot_stats(total_intensity, mean_intensity, std_intensity)
	# Plot mean and std over time
	n_frames = length(total_intensity)
	# Create a 3-row, 1-column plot layout
	plot_layout = @layout [a; b; c]

	p1 = plot(1:n_frames, total_intensity,
          label="Total Intensity", xlabel="Frame", ylabel="Total",
          title="Total Intensity Over Time", lw=1)

	p2 = plot(1:n_frames, mean_intensity,
          label="Mean Intensity", xlabel="Frame", ylabel="Mean",
          title="Mean Intensity Over Time", lw=1)

	p3 = plot(1:n_frames, std_intensity,
          label="Std Dev", xlabel="Frame", ylabel="Std",
          title="Standard Deviation Over Time", lw=1)

	# Combine them in a single figure
	plot(p1, p2, p3, layout=plot_layout, size=(700, 900))
end

# ╔═╡ b01acdf0-9a1e-4021-bf72-baf9a80be1e0
"""
Get the time trace for a given pixel 
"""
function pixel_trace(imstack::Array{Float32, 3}, i::Int, j::Int)
    # Get the number of frames (3rd dimension)
    n_frames = size(imstack, 3)
    
    # Extract intensity at (i,j) across all frames
    trace = [imstack[i, j, t] for t in 1:n_frames]
    
    return trace
end


# ╔═╡ b6847d1f-b953-4ed0-b63b-ffc86d053685
imst = getfiles(rdir, sample, field; pixelsize=330)

# ╔═╡ 4e220360-8eeb-4743-8ea1-4a38eaaadad1
total_intensity, mean_intensity, std_intensity = get_stats(imst)

# ╔═╡ d96023e9-9540-4b83-9a33-6898eb73764f
plot_stats(total_intensity, mean_intensity, std_intensity)

# ╔═╡ 04ca3b4b-06ea-48c8-9bd0-564f8cc286cb
begin
	trace = pixel_trace(imst, 63, 32)
	FitX, Fits, S_curves, best_shots, steptable = auto_stepfinder(trace; N_iter=30)
end


# ╔═╡ 8055159d-cf4b-4ace-a6c8-32865d43ddfc
typeof(trace)

# ╔═╡ af7cc72a-cf5b-4430-8890-031457af5d3b
begin
	p1 = plot(trace, xlabel="Frame", ylabel="Intensity",
		 title="Pixel Intensity Trace", label="trace")
	p1 = plot!(p1, FitX, label = "fit")
end

# ╔═╡ 708b25fc-fa8d-4b11-bd67-87a4f66f7768
begin
	step_positions = steptable[:, 1]
	step_heights = steptable[:, 4]
	bar(step_positions, step_heights, label="Steps", xlabel="Index", ylabel="Step Size")
end

# ╔═╡ ff866e71-7f47-41f8-bdde-25c1e3e1978c
function get_traces(imst; imin=1, imax=300, jmin=1, jmax=330)
	for i in imin:imax
		for j in jmin:jmax
			trace = pixel_trace(imst, i, j)
			#debug(" mean trace = $(mean(trace)), std trace =$(std(trace))")
			FitX, Fits, S_curves, best_shots, steptable = auto_stepfinder(trace; N_iter=30)
			if sum(FitX) > 0
				println("found step for pixel (i,j) = ($(i), $(j))")
				step_positions = steptable[:, 1]
				step_heights = steptable[:, 4]
				println("step_positions = $(step_positions)")
				println("step_heights = $(step_heights)")
			end
		end
	end
end

# ╔═╡ cfea7a8b-8e7d-452a-bcf8-d53264b77317
get_traces(imst; imin=1, imax=100, jmin=1, jmax=100)

# ╔═╡ Cell order:
# ╠═9292025e-0d7d-11f0-365e-f1724fc39b3c
# ╠═5b3b5f83-058a-4582-9482-9b2fc5892692
# ╠═f8573e77-38d8-467c-a123-9a83f73e8970
# ╠═19476521-78bd-44d6-a302-6e971f9fc095
# ╠═583a9aee-08eb-4f5a-95ef-d0087eb98cbc
# ╠═39b011c6-f511-42dd-befc-eaf3fd17ea1a
# ╠═8e8cc911-5a3c-4ec5-9498-cd60a09fd276
# ╠═b7508873-d85d-4502-8845-b779feb4666c
# ╠═88b1aca9-3390-40c2-a9c6-cdce13528b07
# ╠═89114d7b-298d-41cb-bc48-bdf32865a735
# ╠═7009519c-c468-4547-9d26-32e54e90404f
# ╠═fe61875d-8e06-48b0-b67b-5ed01b6dcc6a
# ╠═b01acdf0-9a1e-4021-bf72-baf9a80be1e0
# ╠═b6847d1f-b953-4ed0-b63b-ffc86d053685
# ╠═4e220360-8eeb-4743-8ea1-4a38eaaadad1
# ╠═d96023e9-9540-4b83-9a33-6898eb73764f
# ╠═04ca3b4b-06ea-48c8-9bd0-564f8cc286cb
# ╠═8055159d-cf4b-4ace-a6c8-32865d43ddfc
# ╠═af7cc72a-cf5b-4430-8890-031457af5d3b
# ╠═708b25fc-fa8d-4b11-bd67-87a4f66f7768
# ╠═ff866e71-7f47-41f8-bdde-25c1e3e1978c
# ╠═cfea7a8b-8e7d-452a-bcf8-d53264b77317
