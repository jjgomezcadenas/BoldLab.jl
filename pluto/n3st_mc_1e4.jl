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
	using JStepFinder
	using StepAnalysis
	using histos
	import Measures
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
names(StepAnalysis)

# ╔═╡ b5f399bf-5713-4f26-afb0-2d5771dbbc6f
names(JStepFinder)

# ╔═╡ 11e43f25-3fa8-4f2d-bba5-1773a4989178
names(histos)

# ╔═╡ 8e8cc911-5a3c-4ec5-9498-cd60a09fd276
set_log_level(SimpleLogger.DEBUG)

# ╔═╡ b7508873-d85d-4502-8845-b779feb4666c
debug("test")

# ╔═╡ f41f04da-5666-4c28-97c4-16eba8e6e7cf
function select_file(mc, noise, stats)
	xmc = "n3d_mc"
	if mc
		if noise 
			xf = string(xmc, "_", stats, "_noise.npy")
		else
			xf = string(xmc, "_", stats, "_no_noise.npy")
		end
		println("reading file $(xf)")
		npzread(xf)
	else
		imst = getfiles(rdir, sample, field; pixelsize=330)
	end
	
	
end

# ╔═╡ 84678b4f-847b-4883-ab0d-3086468e30b8
md"""
## Analysis 
"""

# ╔═╡ 77fa8f9e-2b6f-4743-a4ad-247c9e08ef54
readdir()

# ╔═╡ b6847d1f-b953-4ed0-b63b-ffc86d053685
begin
	MC = true 
	Noise = true
	Stats = "1e4"
	imst = select_file(MC, Noise, Stats)
	total_intensity, mean_intensity, std_intensity = get_stats(imst)
	plot_stats(total_intensity, mean_intensity, std_intensity)
end


# ╔═╡ cfbf4c54-78ee-46e6-ace2-0877e04543a5
begin
	xsum, xmean, xstd = traces_stats(imst)
	#plot_traces_stats(xsum, xmean, xstd; meanmx=250.0, summx=5e+4, stdmx=300.0)
	plot_traces_h2d(xsum, xmean, xstd, bins=50)
end

# ╔═╡ 6b54a654-7850-4a27-8ec4-2dcbd3566d82
#hsum2d, hmean2d = traces_h2d(xsum, xmean, xstd; bins=50)

# ╔═╡ 74deedfc-a344-4373-9fc9-db22e83d48ac
md"""
#### Get traces from the imst matrix
"""

# ╔═╡ 4d1615ff-570d-4ae2-86c7-5be653ae6bd0
begin
	trz = traces_thr(imst; xmthr=0.0, xsthr=0.0)
	md"""
	Vector of traces has lengh $(length(trz))
	"""
end

# ╔═╡ 3a8e7300-d209-4e50-ace8-c4b841d71f42
md"""
#### Plot traces
"""

# ╔═╡ 08dd1429-4c2d-4def-8282-e7abe469f318
md"""
### Fit data
"""

# ╔═╡ e28e14f4-8221-4224-b0a6-1aa5ec598227
function plot_trx(DX::Vector{Vector{Float64}}, 
				  FX::Vector{Vector{Float64}};
                  nx=4, ny=4, figsize=(1500, 500))
	
    nplots = min(length(DX), nx * ny)
    plots_array = Vector{Any}(undef, nplots)

    for i in 1:nplots
        plt = plot(1:length(DX[i]), DX[i],
				   label="Noisy Signal", color=:gray, lw=1,
                   xlabel="time steps", ylabel="Intensity", title="Trace", 
				   legend=:topright, grid=true)
		
        plot!(plt, 1:length(FX[i]), FX[i], label="Step Fit thr=0", 
			  color=:blue, lw=2)
        plots_array[i] = plt
    end

    # Fill unused plots if any
    for i in nplots+1:nx*ny
        plots_array[i] = plot(title="(empty)")
    end

    plot(plots_array...; layout=(nx, ny), size=figsize)
end

# ╔═╡ 241d979d-42f4-41bc-be0e-f3f8fbe3141f
function fit_traces(trz::Vector{Vector{T}};
					n::Int=1, niter::Int=10) where {T<:AbstractFloat}
    dataX = trz[n]
    thx = range(0.05, 0.5; length=30)

    FitX0, _, _, S_curve, best_shot0 = stepfindcore(dataX; tresH=0.0, 
													N_iter=niter, demo=false)

	
    thr = 0.0
    for (i, t) in enumerate(thx)
        FitX, _, _, S_curve, best_shot = stepfindcore(dataX; tresH=t,
															  N_iter=niter, demo=false)
        if best_shot == 0
            thr = i == 1 ? t : thx[i-1]
            break
        end
    end

    if thr > 0
        FitX, _, _, S_curve, best_shot = stepfindcore(dataX; 
													  tresH=thr, N_iter=niter,
											        demo=false)
    end

    return dataX, FitX0, FitX, S_curve, thr, best_shot 
end

# ╔═╡ 24c6456d-922a-4ea2-aa62-5250f8561598
function getsteps(arr::Vector{Float64}; atol=1e-8)
    # Find where the difference between consecutive elements is not close to zero
    change_mask = .!isapprox.(diff(arr), 0.0, atol=atol)
    change_indices = findall(identity, change_mask) .+ 1

    # Get the constant values of each segment
    segment_values = arr[change_indices]
    constant_values = vcat(arr[1], segment_values)

    # Compute segment start and end indices
    segment_starts = vcat(1, change_indices)
    segment_ends = vcat(change_indices, length(arr) + 1)
    segment_lengths = segment_ends .- segment_starts

    return constant_values, segment_starts, segment_lengths
end

# ╔═╡ a2fc0f98-6c5e-406d-bae4-7d8fdc8b51a1
function fit_data(trz::Vector{Vector{Float64}}; 
				  niter=10, thr=0.0, rmax=1000, fplot=50)
	
    DX = Vector{Vector{Float64}}()
    FX = Vector{Vector{Float64}}()
    StepHeight = Float64[]
    StepTime = Int[]
    StepLength = Int[]
    TL = Int[]

    nmax = min(rmax, length(trz))

    for n in 1:nmax
        dataX = trz[n]
        FitX, _, _, S_curve, best_shot0 = stepfindcore(dataX; 
													   tresH=thr, N_iter=niter, 
													   demo=false)

        if mod(n, fplot) == 0
            push!(DX, dataX)
            push!(FX, FitX)
        end

        sth, stt, stl = getsteps(FitX)
        push!(TL, length(sth))

        if length(sth) > 1
            append!(StepHeight, sth)
            append!(StepTime, stt)
            append!(StepLength, stl)
        end
    end

    return DX, FX, TL, StepHeight, StepTime, StepLength
end

# ╔═╡ f6678c7f-3560-4ca6-9331-9b702be30b17
DX, FX, TL, StepHeight, StepTime, StepLength =  fit_data(trz, niter=10, thr=0.1, rmax=10000, fplot=500)

# ╔═╡ 1c93276e-1de3-4784-8d3f-f520762af68c
plot_trx(DX,FX, nx=5, ny=4, figsize=(1500, 1500))

# ╔═╡ d37ff23d-15f6-4c96-b0ff-c82f6119c03b
plot_steps(TL, StepHeight, StepTime, StepLength,
		   nbins = 50, shlims = (-750.0, 100.0),
		   figsize=(1200, 1000))

# ╔═╡ ac712a9e-eb6b-4d86-afa6-7233ff512b67
begin
	htl, htlp = hist1d(TL, "number of steps", 10, 0, 10, datap=false)
	htlp
end

# ╔═╡ ac554889-acab-4483-bee3-c9cbe6fee9bf
md"""
- Fraction of 1 step (SM) in sample $(htl.weights[2]/sum(htl.weights))
"""

# ╔═╡ 507db1a2-3ddd-4461-a652-057773348648

function plotfit2(dataX::Vector{Float64}, 
				  FitX1::Vector{Float64}, 
				  FitX2::Vector{Float64}, S_curve::Vector{Float64};
                  figsize=(1000, 800))  # Size in pixels, similar to figsize in matplotlib
    
    plt1 = plot(1:length(dataX), dataX, 
				label="Noisy Signal", color=:gray, lw=1,
                xlabel="time steps", ylabel="Intensity", title="Trace", 
				legend=:topright, grid=true)
	
    plot!(plt1, 1:length(FitX1), FitX1, 
		  label="Step Fit thr=0", color=:blue, lw=2)
    plot!(plt1, 1:length(FitX2), FitX2, 
		  label="Step Fit", color=:red, lw=2)

    plt2 = plot(1:length(S_curve), S_curve, 
				marker=:circle, label="S-curve",
                xlabel="Iteration", ylabel="S-value (MSE ratio)", 
				title="Goodness of Fit (S-curve)", 
				grid=true, legend=:topright)

    plot(plt1, plt2, layout=(2, 1), size=figsize)
end

# ╔═╡ a7542de4-5e78-4c55-91b7-d68c4e2f7985
begin
	nt = 1
	dataX, FitX0, FitX1, S_curve, thx, bshot = fit_traces(trz, n=nt, niter=10)
	sth, stt, stl = getsteps(FitX1)
	
	println("threshold = $(thx), best shot =$(bshot)")
	println("Step heights: ", sth)
	println("Step times: ", stt)
	println("Segment lengths: ", stl)
	plotfit2(dataX, FitX0, FitX1, S_curve)
end

# ╔═╡ Cell order:
# ╠═9292025e-0d7d-11f0-365e-f1724fc39b3c
# ╠═5b3b5f83-058a-4582-9482-9b2fc5892692
# ╠═f8573e77-38d8-467c-a123-9a83f73e8970
# ╠═19476521-78bd-44d6-a302-6e971f9fc095
# ╠═583a9aee-08eb-4f5a-95ef-d0087eb98cbc
# ╠═39b011c6-f511-42dd-befc-eaf3fd17ea1a
# ╠═b5f399bf-5713-4f26-afb0-2d5771dbbc6f
# ╠═11e43f25-3fa8-4f2d-bba5-1773a4989178
# ╠═8e8cc911-5a3c-4ec5-9498-cd60a09fd276
# ╠═b7508873-d85d-4502-8845-b779feb4666c
# ╠═f41f04da-5666-4c28-97c4-16eba8e6e7cf
# ╠═84678b4f-847b-4883-ab0d-3086468e30b8
# ╠═77fa8f9e-2b6f-4743-a4ad-247c9e08ef54
# ╠═b6847d1f-b953-4ed0-b63b-ffc86d053685
# ╠═cfbf4c54-78ee-46e6-ace2-0877e04543a5
# ╠═6b54a654-7850-4a27-8ec4-2dcbd3566d82
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═4d1615ff-570d-4ae2-86c7-5be653ae6bd0
# ╠═3a8e7300-d209-4e50-ace8-c4b841d71f42
# ╠═a7542de4-5e78-4c55-91b7-d68c4e2f7985
# ╠═08dd1429-4c2d-4def-8282-e7abe469f318
# ╠═f6678c7f-3560-4ca6-9331-9b702be30b17
# ╠═1c93276e-1de3-4784-8d3f-f520762af68c
# ╠═d37ff23d-15f6-4c96-b0ff-c82f6119c03b
# ╠═ac712a9e-eb6b-4d86-afa6-7233ff512b67
# ╠═ac554889-acab-4483-bee3-c9cbe6fee9bf
# ╠═e28e14f4-8221-4224-b0a6-1aa5ec598227
# ╠═a2fc0f98-6c5e-406d-bae4-7d8fdc8b51a1
# ╠═241d979d-42f4-41bc-be0e-f3f8fbe3141f
# ╠═24c6456d-922a-4ea2-aa62-5250f8561598
# ╠═507db1a2-3ddd-4461-a652-057773348648
