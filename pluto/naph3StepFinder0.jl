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
	using NPZ
end

# ╔═╡ 5b3b5f83-058a-4582-9482-9b2fc5892692
begin
	using PlutoUI
	using CSV
	using DataFrames
	#using Plots 
	using StatsPlots 
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
import Measures

# ╔═╡ abc940c8-68ac-466f-af96-1315dc89317a
5 * Measures.mm

# ╔═╡ 7009519c-c468-4547-9d26-32e54e90404f
"""
Computes the gloabl statistics of the stack. Sum, mean and stds. Also creates
a vector of means and stds.
imst is the stack  
"""
function get_stats(imst::Array{T, 3}; bck=200.0) where {T<:AbstractFloat}
	
	n_frames = size(imst, 3)
	
	# Preallocate vectors
	total_intensity = zeros(T, n_frames)
	mean_intensity = zeros(T, n_frames)
	std_intensity = zeros(T, n_frames)

	
	# Loop over each image in the stack
	for i in 1:n_frames
	    frame = imst[:, :, i] .- bck
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
	#plot(p1, p2, p3, layout=plot_layout, size=(700, 900))
	plot(p1, p2, p3;
     layout=(3, 1),
     size=(700, 900),
     margin=5 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=10,
     guidefontsize=9,
     tickfontsize=8,
     legendfontsize=8)
end

# ╔═╡ b01acdf0-9a1e-4021-bf72-baf9a80be1e0
"""
Get the time trace for a given pixel 
"""
function pixel_trace(imstack::Array{T, 3}, i::Int, j::Int; bck=200.0) where {T<:AbstractFloat}
    # Get the number of frames (3rd dimension)
    n_frames = size(imstack, 3)
    
    # Extract intensity at (i,j) across all frames
    trace = [imstack[i, j, t] .- bck for t in 1:n_frames]
    
    return trace
end


# ╔═╡ f41f04da-5666-4c28-97c4-16eba8e6e7cf


# ╔═╡ ff866e71-7f47-41f8-bdde-25c1e3e1978c
function get_traces(imst::Array{T, 3}; imin=1, imax=300, jmin=1, jmax=330) where {T<:AbstractFloat}
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

# ╔═╡ 04ca3b4b-06ea-48c8-9bd0-564f8cc286cb
function traces_stats(imst::Array{T, 3}) where {T<:AbstractFloat}
	imx = size(imst)[1]
	jmx = size(imst)[2]
	xsum = zeros(imx,jmx)
	xmean = zeros(imx,jmx)
	xstd = zeros(imx,jmx)
	for i in  1:imx
		for j in  1:jmx
			trace = pixel_trace(imst, i, j)
			xsum[i,j] = sum(trace)
			xmean[i,j] = mean(trace)
			xstd[i,j] =  std(trace)
		end
	end
	xsum, xmean, xstd
	#FitX, Fits, S_curves, best_shots, steptable = auto_stepfinder(trace; N_iter=30)
end


# ╔═╡ 3f4f0f81-ec9b-48a2-839b-3e9cc569b53f
function traces_thr(imst::Array{T, 3}; xmthr=0.0, xsthr=10.0) where {T<:AbstractFloat}
	imx = size(imst)[1]
	jmx = size(imst)[2]

	M = Matrix{Vector{T}}(undef, imx, jmx) 
	for i in  1:imx
		for j in  1:jmx
			trace = pixel_trace(imst, i, j)
			xmean = mean(trace)
			xstd =  std(trace)
			if xmean > xmthr && xstd > xsthr
				M[i,j] = trace
			else
				M[i,j] = [0.0]
			end
		end
	end
	[v for v in M if v != [0.0]]
end

# ╔═╡ 96b62a2c-15d0-4ed4-a34d-dc7428d57af9
function plot_traces_stats(xsum, xmean, xstd; meanmx=25.0, summx=1e+3, stdmx=50.0)
	xmin = 0.0
	xmax = meanmx
	data = vec(xmean)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p1 = histogram(filtered, bins=collect(edges), title="Mean of traces", xlabel="trace index", ylabel="Frequency")

	xmin = 0.0
	xmax = summx
	data = vec(xsum)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p2 = histogram(filtered, bins=collect(edges), title="Sum of traces", xlabel="trace index", ylabel="Frequency")

	xmin = 0.0
	xmax = stdmx
	data = vec(xstd)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p3 = histogram(filtered, bins=collect(edges), title="std of traces", xlabel="trace index", ylabel="Frequency")

	plot(p1, p2, p3;
     layout=(3, 1),
     size=(700, 900),
     margin=5 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=10,
     guidefontsize=9,
     tickfontsize=8,
     legendfontsize=8)

	
end

# ╔═╡ 2d10dd39-5192-469b-bd5d-4ec3dc38f337
function traces_h2d(xsum, xmean, xstd; bins=50)
	function h2d(x1, x2)
		# Ensure matrices have the same shape
		@assert size(x1) == size(x2)
	
		# Flatten both to vectors
		x = vec(x1)
		y = vec(x2)

		# Compute the 2D histogram with 50×50 bins
		h = fit(Histogram, (x, y), nbins=(bins, bins))

	end
	h2d(xsum, xstd), h2d(xmean, xstd)
end

# ╔═╡ 54a0af5d-988e-45a4-819f-a9c2a87a47e6
function plot_traces_h2d(xsum, xmean, xstd; bins=50)
	function h2d(x1, x2, x1l, x2l, xt)
		# Ensure matrices have the same shape
		@assert size(x1) == size(x2)
	
		# Flatten both to vectors
		x = vec(x1)
		y = vec(x2)
	
		# 2D histogram
		histogram2d(x, y;
		       bins=bins,
		       xlabel=x1l,
		       ylabel=x2l,
		       title=xt)
	end
	p1 = h2d(xsum, xstd, "xsum", "xstd", "Sum vs Std")
	p2 = h2d(xmean, xstd, "mean", "xstd", "Mean vs Std")
	plot(p1, p2;
     layout=(1, 2),
     margin=1 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=8,
     guidefontsize=9,
     tickfontsize=6,
     legendfontsize=8)
end

# ╔═╡ 84678b4f-847b-4883-ab0d-3086468e30b8
md"""
## Analysis 
"""

# ╔═╡ b6847d1f-b953-4ed0-b63b-ffc86d053685
begin
	MC = true 
	if MC
		imst = npzread("nap3d.npy")
	else
		imst = getfiles(rdir, sample, field; pixelsize=330)
	end
	total_intensity, mean_intensity, std_intensity = get_stats(imst)
	plot_stats(total_intensity, mean_intensity, std_intensity)
end


# ╔═╡ 58a8b07c-cf1b-498f-925b-f40973f1a9d5


# ╔═╡ 7f737b3f-7980-4fb9-9a97-c2a2bdf2ffc9
xsum, xmean, xstd = traces_stats(imst)

# ╔═╡ 32505f27-b184-4df1-af46-9619c4824232
xstd

# ╔═╡ a8c7e174-8913-47ff-a4d4-4b17954b3c02
plot_traces_stats(xsum, xmean, xstd; meanmx=250.0, summx=5e+4, stdmx=300.0)

# ╔═╡ cfbf4c54-78ee-46e6-ace2-0877e04543a5
plot_traces_h2d(xsum, xmean, xstd, bins=50)

# ╔═╡ 6b54a654-7850-4a27-8ec4-2dcbd3566d82
hsum2d, hmean2d = traces_h2d(xsum, xmean, xstd; bins=50)

# ╔═╡ a5880b0f-1996-4760-9e1a-9d88ec2b66fe
hsum2d.weights

# ╔═╡ ab92b991-1c40-48b6-8524-7f55c19dec27
hsum2d.edges[2]

# ╔═╡ b88350ba-400c-454c-b310-d6a249f03660
hsum2d.edges[1]

# ╔═╡ 4d1615ff-570d-4ae2-86c7-5be653ae6bd0
trz = traces_thr(imst; xmthr=0.0, xsthr=10.0)

# ╔═╡ 964d94ad-22b3-4468-ac11-46d31593b902
function plt_ft_trace(trz, nt; fit=false, niter=30, thr=1.0)
	
	pt1 =plot(trz[nt], xlabel="trace number $(nt)", ylabel="Intensity",
		 title="Pixel Intensity Trace", label="trace")

	if fit
		FitX, Fits, S_curves, best_shots, steptable = auto_stepfinder(trz[nt]; N_iter=niter, tresH=thr)
		pt1 = plot!(pt1, FitX, label = "fit")
	end
	FitX, steptable, pt1
end

# ╔═╡ 062f2b66-7080-4e89-bc39-90ab8363a480
begin
	FitX, steptable, pt1 = plt_ft_trace(trz, 4; fit=true, niter=30, thr=0.3)
	pt1
end

# ╔═╡ 46e2907d-41e4-40bf-86f8-07d4a784f966


# ╔═╡ 708b25fc-fa8d-4b11-bd67-87a4f66f7768
begin
	step_positions = steptable[:, 1]
	step_heights = steptable[:, 4]
	bar(step_positions, step_heights, label="Steps", xlabel="Index", ylabel="Step Size")
end

# ╔═╡ ffe345e1-94fd-4f70-bcf7-bfeab6e732b6
FitX

# ╔═╡ 7a3dae87-bfa3-4e95-81b1-e119ddf15108
histogram(step_positions, bins=10, title="step positions", xlabel="step position", ylabel="Frequency")

# ╔═╡ 146ffd87-8a28-4589-9bc4-840216ec4d92
histogram(step_heights, bins=10, title="step heights", xlabel="step height", ylabel="Frequency")

# ╔═╡ 73e198cc-b1a6-4594-9354-5dc3a9d96648
function get_traces(trz::Vector{Vector{T}}; niter=30, thr=0.3) where {T<:AbstractFloat}
    step_positions = Float64[]
    step_heights = Float64[]

    for nt in 1:length(trz)
        try
            FitX, Fits, S_curves, best_shots, steptable = auto_stepfinder(trz[nt]; N_iter=niter, tresH=thr)
            append!(step_positions, steptable[:, 1])
            append!(step_heights, steptable[:, 4])
        catch e
            if e isa ArgumentError
                @warn "Skipped trace $nt due to ArgumentError: $(e.msg)"
            else
                rethrow(e)  # let other exceptions bubble up
            end
        end
    end

    return step_positions, step_heights
end

# ╔═╡ ec8c7b95-a59a-4804-ae46-01519102f82d
stpx, stph = get_traces(trz; niter=30, thr=0.3)

# ╔═╡ 68f8312e-a4df-4f1d-983d-9834586ec6e9
typeof(trz)

# ╔═╡ cfea7a8b-8e7d-452a-bcf8-d53264b77317
begin
	xmin=-5e+2
	xmax=5e+2
	stpf = stph[(stph .>= xmin) .& (stph .<= xmax)]
	histogram(stpf, bins=20, title="step heights", xlabel="step height", ylabel="Frequency")
end

# ╔═╡ 16f31587-317b-4481-9354-7ac2494b20be
minimum(stph)

# ╔═╡ b2b19c65-b45d-4404-8ce4-847c813596b3
histogram(stpx, bins=20, title="step position", xlabel="step position", ylabel="Frequency")

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
# ╠═abc940c8-68ac-466f-af96-1315dc89317a
# ╠═7009519c-c468-4547-9d26-32e54e90404f
# ╠═fe61875d-8e06-48b0-b67b-5ed01b6dcc6a
# ╠═b01acdf0-9a1e-4021-bf72-baf9a80be1e0
# ╠═f41f04da-5666-4c28-97c4-16eba8e6e7cf
# ╠═ff866e71-7f47-41f8-bdde-25c1e3e1978c
# ╠═04ca3b4b-06ea-48c8-9bd0-564f8cc286cb
# ╠═3f4f0f81-ec9b-48a2-839b-3e9cc569b53f
# ╠═96b62a2c-15d0-4ed4-a34d-dc7428d57af9
# ╠═2d10dd39-5192-469b-bd5d-4ec3dc38f337
# ╠═54a0af5d-988e-45a4-819f-a9c2a87a47e6
# ╠═84678b4f-847b-4883-ab0d-3086468e30b8
# ╠═b6847d1f-b953-4ed0-b63b-ffc86d053685
# ╠═58a8b07c-cf1b-498f-925b-f40973f1a9d5
# ╟─7f737b3f-7980-4fb9-9a97-c2a2bdf2ffc9
# ╠═32505f27-b184-4df1-af46-9619c4824232
# ╠═a8c7e174-8913-47ff-a4d4-4b17954b3c02
# ╠═cfbf4c54-78ee-46e6-ace2-0877e04543a5
# ╠═6b54a654-7850-4a27-8ec4-2dcbd3566d82
# ╠═a5880b0f-1996-4760-9e1a-9d88ec2b66fe
# ╠═ab92b991-1c40-48b6-8524-7f55c19dec27
# ╠═b88350ba-400c-454c-b310-d6a249f03660
# ╠═4d1615ff-570d-4ae2-86c7-5be653ae6bd0
# ╠═964d94ad-22b3-4468-ac11-46d31593b902
# ╠═062f2b66-7080-4e89-bc39-90ab8363a480
# ╠═46e2907d-41e4-40bf-86f8-07d4a784f966
# ╠═708b25fc-fa8d-4b11-bd67-87a4f66f7768
# ╠═ffe345e1-94fd-4f70-bcf7-bfeab6e732b6
# ╠═7a3dae87-bfa3-4e95-81b1-e119ddf15108
# ╠═146ffd87-8a28-4589-9bc4-840216ec4d92
# ╠═73e198cc-b1a6-4594-9354-5dc3a9d96648
# ╠═ec8c7b95-a59a-4804-ae46-01519102f82d
# ╠═68f8312e-a4df-4f1d-983d-9834586ec6e9
# ╠═cfea7a8b-8e7d-452a-bcf8-d53264b77317
# ╠═16f31587-317b-4481-9354-7ac2494b20be
# ╠═b2b19c65-b45d-4404-8ce4-847c813596b3
