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
begin
	MC = true 
	Noise = true
	Stats = "1e3"
	
	md"""
	## Analysis 
	
	- MC = $MC 
	- Noise = $Noise
	- Stats = $(Stats) molecules in 30 x 30 micron. 
	"""
end

# ╔═╡ 77fa8f9e-2b6f-4743-a4ad-247c9e08ef54
readdir()

# ╔═╡ b6847d1f-b953-4ed0-b63b-ffc86d053685
begin
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

# ╔═╡ a7542de4-5e78-4c55-91b7-d68c4e2f7985
begin
	nt = 1
	dataX, FitX0, FitX1, S_curve, thx, bshot = fit_traces(trz, n=nt, niter=10)
	sth, stt, stl = getsteps(FitX0)
	
	println("threshold = $(thx), best shot =$(bshot)")
	println("Step heights: ", sth)
	println("Step times: ", stt)
	println("Segment lengths: ", stl)
	plotfit2(dataX, FitX0, FitX1, S_curve)
end

# ╔═╡ 08dd1429-4c2d-4def-8282-e7abe469f318
md"""
### Fit data
"""

# ╔═╡ f6678c7f-3560-4ca6-9331-9b702be30b17
#DX, FX, TL, StepHeight, StepTime, StepLength =  fit_data(trz, niter=10, thr=0.05, rmax=120000, fplot=400)

# ╔═╡ 12af547f-37c6-4954-b512-08b369b00966
#length(DX)

# ╔═╡ 1c93276e-1de3-4784-8d3f-f520762af68c
#plot_trx(DX,FX, nx=5, ny=4, figsize=(1500, 1500))

# ╔═╡ d37ff23d-15f6-4c96-b0ff-c82f6119c03b
#plot_steps(TL, StepHeight, StepTime, StepLength,
#		   nbins = 50, shlims = (-750.0, 100.0),
#		   figsize=(1200, 1000))

# ╔═╡ 4429672c-d953-41f3-b43e-beb0e285e680
#TL2 .-1

# ╔═╡ 8e8fd656-17df-44d0-9633-474155d7dc7f
md"""
- remove blinking (negative steps)
"""

# ╔═╡ b40f0395-9252-4381-9ac9-32e7d9084a12
md"""
- Only stories with at least 1 step
"""

# ╔═╡ 614c7468-34a4-4d39-8820-4d6aef99bd83
md"""
- Compute min step height for each story. This is done to define a baseline
"""

# ╔═╡ 267bd80d-165a-4445-90d8-ab70fd716569


# ╔═╡ 309c89fa-a70d-43fa-adc2-6eacd10369f8
begin
	#sdf4 = filter(row -> row.nstep ==2, sdf3)
end

# ╔═╡ 1f908e9b-855a-4ecc-a8d7-8ef92b502c54
md"""
## Functions
"""

# ╔═╡ 9893562b-800d-4384-aa91-d45e10de1cb6
function renumber_nmol!(df::DataFrame)
    unique_ids = unique(df.nmol)
    id_map = Dict(id => i for (i, id) in enumerate(sort(unique_ids)))
    df.nmol .= [id_map[id] for id in df.nmol]
    return df
end

# ╔═╡ 418f05eb-9b20-418f-9a35-783eac94e501
function remove_bad_molecules(df::DataFrame)
    grouped = groupby(df, :nmol)
    filtered_groups = []

    for g in grouped
        if nrow(g) >= 2 && g.stepHeight[1] ≥ g.stepHeight[2]
            push!(filtered_groups, g)
        end
    end

    sdf = vcat(filtered_groups...)
	renumber_nmol!(sdf)
end

# ╔═╡ 776c8d57-15c9-45e1-8ddb-c70afe4240a2

"""
    plot_trace(df::DataFrame, mol_id::Int)

Plot the trace of a given molecule from a dataframe.

# Arguments
- `df`: DataFrame with columns `nmol`, `nstep`, `stepHeight`, `stepTime`, `stepLength`.
- `mol_id`: Integer representing the renumbered molecule index to plot.

# Description
Assumes that `nmol` values have been renumbered to start from 1. Plots a step function with heights from `stepHeight`,
transitions at `stepTime`, and final segment length taken from `stepLength`.
"""
function plot_trace(df::DataFrame, mol_id::Int)
    # Filter rows for the specified molecule
    mol_df = filter(:nmol => ==(mol_id), df)
    sort!(mol_df, :stepTime)

    xs = Float64[]
    ys = Float64[]

    for i in 1:nrow(mol_df)
        h = mol_df.stepHeight[i]
        t_start = mol_df.stepTime[i]
        if i < nrow(mol_df)
            t_end = mol_df.stepTime[i+1]
        else
            t_end = t_start + mol_df.stepLength[i]
        end
        append!(xs, [t_start, t_end])
        append!(ys, [h, h])
    end

    plot(xs, ys, label="Trace $mol_id", lw=2, xlabel="Time", ylabel="Step Height", title="")
end

# ╔═╡ b46f0c62-45f4-44ca-946c-bc9447306c45
function plot_traces(df::DataFrame, mol_ids::Vector{Int}; size=(1200, 800))
	P = []
	for mol in mol_ids
		push!(P, plot_trace(df, mol))
	end
	plot(P..., size=size)
end

# ╔═╡ 4f558c30-f5ba-4b15-a72a-da16939c0bf5
function plot_sdf(sdf)
	pnst = step_hist(sdf.nstep .-1;
	                   nbins= 10,
	                   xlim = (0.0,10.0),
	                   xlabel = "number of steps",
	                   ylabel= "Frequency")
		
	psth = step_hist(sdf.stepHeight;
	                   nbins= 50,
	                   xlim = (-100.0, 800.0),
					   logy=true,
	                   xlabel = "step heigth",
	                   ylabel= "Frequency")
	pstt = step_hist(sdf.stepTime;
	                   nbins= 50,
	                   xlim = (0.0, 200.0),
					   logy=true,
	                   xlabel = "step time",
	                   ylabel= "Frequency")
	pstl = step_hist(sdf.stepLength;
	                   nbins= 50,
	                   xlim = (0.0, 200.0),
					   logy=true,
	                   xlabel = "step length",
	                   ylabel= "Frequency")

	plot(pnst, psth, pstt, pstl, layout=(2, 2), size=(1000, 800))
	
end

# ╔═╡ 8b5448b2-73b2-4c9b-9f64-398a7a7a7720
"""
    fit_data2(trz::Vector{Vector{Float64}}; niter=10, thr=0.0, rmax=1000, fplot=50)

Fits step functions to a vector of time traces using the `stepfindcore` algorithm 
and collects statistics of each detected step in a DataFrame.

# Arguments
- `trz`: A vector of time traces (`Vector{Vector{Float64}}`), where each trace corresponds to a pixel or measurement location.
- `niter`: Number of iterations used in the step-fitting algorithm (default: 10).
- `thr`: Threshold value for step detection (default: 0.0).
- `rmax`: Maximum number of traces to process (default: 1000).
- `fplot`: Frequency at which traces are selected for plotting (default: every 50th trace).

# Returns
- `DX`: Vector of raw traces selected for plotting.
- `FX`: Vector of fitted traces corresponding to `DX`.
- `TL`: Vector containing the number of steps found in each trace.
- `df`: A `DataFrame` with one row per detected step and the following columns:
    - `nstep`: Number of steps in the corresponding trace.
    - `stepHeight`: Height of the individual step.
    - `stepTime`: Time index at which the step occurs.
    - `stepLength`: Length of the segment following the step.
"""
function fit_data2(trz::Vector{Vector{Float64}}; 
                   niter=10, thr=0.0, rmax=1000, fplot=50)

    DX = Vector{Vector{Float64}}()
    FX = Vector{Vector{Float64}}()
	TL = Int[]

    df = DataFrame(nmol=Int[], nstep=Int[],
				   stepHeight=Float64[], stepTime=Int[], stepLength=Int[])

    nmax = min(rmax, length(trz))

    for n in 1:nmax
        dataX = trz[n]
        FitX, _, _, S_curve, best_shot0 = stepfindcore(dataX; 
                                                       tresH=thr, N_iter=niter, 
                                                       demo=false)
        sth, stt, stl = getsteps(FitX)
        nsteps = length(sth)

		
        if mod(n, fplot) == 0
            push!(DX, dataX)
            push!(FX, FitX)
			push!(TL, nsteps)
        end


        for i in 1:nsteps
            push!(df, (n, nsteps, sth[i], stt[i], stl[i]))
        end
        
    end
    return DX, FX, TL, df
end



# ╔═╡ 0ef14507-a762-445a-a506-ba05bf469456
begin
	DX2, FX2, TL2, sdf =  fit_data2(trz, niter=10, thr=0.05, rmax=120000, fplot=400)
	plot_trx(DX2,FX2, nx=5, ny=4, figsize=(1500, 1500))
end

# ╔═╡ afd8b98c-aee7-44bb-96ac-3e7acd7b2f00
let
	htl, htlp = hist1d(sdf.nstep .-1, "number of steps", 10, 0, 10, datap=false)
	htl
	xf0 = htl.weights[1]/sum(htl.weights)
	xf1 = sum(htl.weights[2:4])/sum(htl.weights)
	md"""
	- Fraction of zero steps =$(xf0)
	- Fraction of one/three steps =$(xf1)
	- number of mol trajectories $(length(unique(sdf.nmol)))
	"""
end

# ╔═╡ 4fc30548-5366-4663-9ec1-d5b6529d4e3d
sdf

# ╔═╡ 5fe515e5-23af-4e04-a1da-944c6f464062
plot_sdf(sdf)

# ╔═╡ 94562765-2067-4111-b343-206f404cc599
sdf1 = filter(row -> row.stepHeight > 0, sdf)

# ╔═╡ 7ea511fb-c91a-4b2b-9914-3a4fb29540eb
#sdf1 = filter(row -> 0 < row.nstep .-1 < 10, sdf)
sdf2 = filter(row -> row.nstep .-1 > 0, sdf1)

# ╔═╡ 5b4f4962-9118-4c85-b005-08c9f46324a2
begin
	sdf3 = filter(row -> row.stepHeight >5, sdf2)
	sdf3 = transform(groupby(sdf3, :nmol), nrow => :nstep)
end

# ╔═╡ 4014329c-e371-4ef8-a77e-738a668bf790
let
	htl = hist1d(sdf3.nstep .-1, 10, 0, 10)
	htl
	xf0 = htl.weights[1]/sum(htl.weights)
	xf1 = sum(htl.weights[2:4])/sum(htl.weights)
	md"""
	- Fraction of zero steps =$(xf0)
	- Fraction of one/three steps =$(xf1)
	"""
	plot(step_hist(sdf3.nstep .-1;
              nbins=10,
              xlim=(0.0, 10.0),
              logy=false,
              xlabel="# of steps",
              ylabel=" #entries "), size=(400,200))
end

# ╔═╡ c7f19b2c-921f-4126-8b09-4c3ac2186b81
renumber_nmol!(sdf3)

# ╔═╡ 7c41be24-ae9a-4a30-95da-f7f851d7a859
plot_traces(sdf3, collect(25:49), size=(1200, 1200))

# ╔═╡ 39aac501-6d96-4ab2-a5bc-ca8b2376a8a1
sdff = remove_bad_molecules(sdf3)

# ╔═╡ 07bf3357-72a7-4acd-acd0-3af77d1c35f9
plot_traces(sdff, collect(101:125), size=(1200, 1200))

# ╔═╡ 45e35a41-3c15-4ab6-b266-e1dca5bde06a
plot_sdf(sdff)

# ╔═╡ e3097b6e-5377-4eaf-b648-9fdda139c088
length(unique(sdff.nmol))

# ╔═╡ 52d7b937-1ffd-40cb-9dc5-6c25fcce2de0
plot_sdf(sdf1)

# ╔═╡ e6da6f93-3123-4d52-ad47-8812847fbef8
msthdf = combine(groupby(sdf1, :nmol), :stepHeight => minimum => :min_stepHeight)

# ╔═╡ cc02c3cd-cc1f-4084-ac53-1f274571843b
begin
	hmsh = hist1d(msthdf.min_stepHeight, 100, 0.0, 200.0)
	hmsh
	plot(step_hist(msthdf.min_stepHeight;
              nbins=100,
              xlim=(0.0, 200.0),
              logy=false,
              xlabel="step length",
              ylabel=" #events "), size=(400,200))
end

# ╔═╡ 9aa87760-ac29-49a9-93a7-8993330e683b
let
	hh = get_histo1d(hmsh)
	xf1 = sum(hh.weights[1])/sum(hh.weights)
	
	md"""
	- Fraction of step height in first three bins =$(xf1)
	- cutoff at $(hh.edges[2])
	"""
	
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
# ╠═12af547f-37c6-4954-b512-08b369b00966
# ╠═1c93276e-1de3-4784-8d3f-f520762af68c
# ╠═d37ff23d-15f6-4c96-b0ff-c82f6119c03b
# ╠═0ef14507-a762-445a-a506-ba05bf469456
# ╠═4429672c-d953-41f3-b43e-beb0e285e680
# ╠═afd8b98c-aee7-44bb-96ac-3e7acd7b2f00
# ╠═4fc30548-5366-4663-9ec1-d5b6529d4e3d
# ╠═5fe515e5-23af-4e04-a1da-944c6f464062
# ╠═8e8fd656-17df-44d0-9633-474155d7dc7f
# ╠═94562765-2067-4111-b343-206f404cc599
# ╠═b40f0395-9252-4381-9ac9-32e7d9084a12
# ╠═7ea511fb-c91a-4b2b-9914-3a4fb29540eb
# ╠═52d7b937-1ffd-40cb-9dc5-6c25fcce2de0
# ╠═614c7468-34a4-4d39-8820-4d6aef99bd83
# ╠═e6da6f93-3123-4d52-ad47-8812847fbef8
# ╠═cc02c3cd-cc1f-4084-ac53-1f274571843b
# ╠═267bd80d-165a-4445-90d8-ab70fd716569
# ╠═9aa87760-ac29-49a9-93a7-8993330e683b
# ╠═5b4f4962-9118-4c85-b005-08c9f46324a2
# ╠═4014329c-e371-4ef8-a77e-738a668bf790
# ╠═309c89fa-a70d-43fa-adc2-6eacd10369f8
# ╠═c7f19b2c-921f-4126-8b09-4c3ac2186b81
# ╠═7c41be24-ae9a-4a30-95da-f7f851d7a859
# ╠═39aac501-6d96-4ab2-a5bc-ca8b2376a8a1
# ╠═07bf3357-72a7-4acd-acd0-3af77d1c35f9
# ╠═45e35a41-3c15-4ab6-b266-e1dca5bde06a
# ╠═e3097b6e-5377-4eaf-b648-9fdda139c088
# ╠═1f908e9b-855a-4ecc-a8d7-8ef92b502c54
# ╠═418f05eb-9b20-418f-9a35-783eac94e501
# ╠═9893562b-800d-4384-aa91-d45e10de1cb6
# ╠═b46f0c62-45f4-44ca-946c-bc9447306c45
# ╠═776c8d57-15c9-45e1-8ddb-c70afe4240a2
# ╠═4f558c30-f5ba-4b15-a72a-da16939c0bf5
# ╠═8b5448b2-73b2-4c9b-9f64-398a7a7a7720
