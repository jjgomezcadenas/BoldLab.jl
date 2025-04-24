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
	using StepAnalysis
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

# ╔═╡ 8e8cc911-5a3c-4ec5-9498-cd60a09fd276
set_log_level(SimpleLogger.DEBUG)

# ╔═╡ b7508873-d85d-4502-8845-b779feb4666c
debug("test")

# ╔═╡ f41f04da-5666-4c28-97c4-16eba8e6e7cf


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
# ╠═f41f04da-5666-4c28-97c4-16eba8e6e7cf
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
# ╠═062f2b66-7080-4e89-bc39-90ab8363a480
# ╠═46e2907d-41e4-40bf-86f8-07d4a784f966
# ╠═708b25fc-fa8d-4b11-bd67-87a4f66f7768
# ╠═ffe345e1-94fd-4f70-bcf7-bfeab6e732b6
# ╠═7a3dae87-bfa3-4e95-81b1-e119ddf15108
# ╠═146ffd87-8a28-4589-9bc4-840216ec4d92
# ╠═ec8c7b95-a59a-4804-ae46-01519102f82d
# ╠═68f8312e-a4df-4f1d-983d-9834586ec6e9
# ╠═cfea7a8b-8e7d-452a-bcf8-d53264b77317
# ╠═16f31587-317b-4481-9354-7ac2494b20be
# ╠═b2b19c65-b45d-4404-8ce4-847c813596b3
