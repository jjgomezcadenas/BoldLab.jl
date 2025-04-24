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
	using Unitful
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
	using SparseArrays
end

# ╔═╡ cba7fc75-8363-4ffa-b5a6-6e7d34363813
import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W

# ╔═╡ 583a9aee-08eb-4f5a-95ef-d0087eb98cbc
names(SimpleLogger)

# ╔═╡ 39b011c6-f511-42dd-befc-eaf3fd17ea1a
names(StepAnalysis)

# ╔═╡ b5f399bf-5713-4f26-afb0-2d5771dbbc6f
names(JStepFinder)

# ╔═╡ 11e43f25-3fa8-4f2d-bba5-1773a4989178
names(histos)

# ╔═╡ b3b16805-b5b8-4782-a49c-15029b3a749d
names(BoldLab)

# ╔═╡ 8e8cc911-5a3c-4ec5-9498-cd60a09fd276
set_log_level(SimpleLogger.DEBUG)

# ╔═╡ b7508873-d85d-4502-8845-b779feb4666c
debug("test")

# ╔═╡ 84678b4f-847b-4883-ab0d-3086468e30b8
begin
	MC = true 
	NM = 1e+1
	
	with_noise=true
	laser_power = 5mW
	pbcycles = 1e+5 # photobleaching cycles (nγ abosorbed before pB)
	dkcycles = 1e+6 # dark cycles (nγ abosorbed before DT)
	tdark    =2s    # time in dark states before decaying to ground
	file = generate_filename(NM, laser_power, pbcycles, dkcycles, tdark; with_noise=with_noise)

	
	md"""
	## Analysis 
	
	- MC = $MC 
	- Noise = $with_noise
	- Stats = $(NM) molecules in 30 x 30 micron. 
	- Laser Power = $(laser_power)
	- PB cycles = $(pbcycles)
	- DK cycles = $(dkcycles)
	- file = $file
	"""
end

# ╔═╡ c8dfaaf4-1430-419d-b048-93dccc92885b
if MC
	imst = npzread(file)
	FF = []

    for i in 1:9
        fn = (i-1) * 20 + i
        push!(FF, heatmap(imst[:, :, fn],
            colorbar=false,  # Optional: removes extra space
            title="Frame $fn",
            titlefontsize=7,
            tickfontsize=6,
            guidefontsize=6,
            titlelocation=:left,
            aspect_ratio=:equal))
    end

    plot(FF...;
        layout=(3, 3),
        size=(900, 900),
        margin=1.0*Measures.mm,
        top_margin=1.0*Measures.mm,
        bottom_margin=1.0*Measures.mm,
        left_margin=1.0*Measures.mm,
        right_margin=1.0*Measures.mm,
        plot_titlefontsize=7,
        legendfontsize=6)
end

# ╔═╡ 77fa8f9e-2b6f-4743-a4ad-247c9e08ef54
readdir()

# ╔═╡ bfc6bdfa-9e79-4be2-b2f9-c264ea27cfb4
begin
	total_intensity, mean_intensity, std_intensity = get_stats(imst)
	plot_stats(total_intensity, mean_intensity, std_intensity)
end

# ╔═╡ cfbf4c54-78ee-46e6-ace2-0877e04543a5
begin
	#xsum, xmean, xstd = traces_stats(imst)
	#plot_traces_stats(xsum, xmean, xstd; meanmx=250.0, summx=5e+4, stdmx=300.0)
	#plot_traces_h2d(xsum, xmean, xstd, bins=50)
end

# ╔═╡ 6b54a654-7850-4a27-8ec4-2dcbd3566d82
#hsum2d, hmean2d = traces_h2d(xsum, xmean, xstd; bins=50)

# ╔═╡ 74deedfc-a344-4373-9fc9-db22e83d48ac
md"""
#### Get traces from the imst matrix
"""

# ╔═╡ 709ca0ac-aaba-440c-a3c7-39e8dd0fee2e
begin
	trsm = traces_above_thr(imst; nsigma=2.0)
	is, js, trnz = findnz(trsm)
	md"""
	Vector of traces above threshold has lenght of $(length(trnz))
	"""
end

# ╔═╡ de42a1a5-aa7e-49a6-98b8-d27b3ee16849
begin
	trI =[sum(tr) for tr in trnz]
	plot(step_hist(trI;
              nbins=10,
              xlim=(0.0, 12000.0),
              logy=false,
              xlabel="Intensity",
              ylabel=" #events "), size=(400,200))
end

# ╔═╡ 3a8e7300-d209-4e50-ace8-c4b841d71f42
md"""
#### Plot traces
"""

# ╔═╡ a7542de4-5e78-4c55-91b7-d68c4e2f7985
begin
	nt = 1
	println("pixel i = $(is[nt]), j=$(js[nt]), I = $(sum(trnz[nt]))")
	dataX, FitX0, FitX1, S_curve, thx, bshot = fit_traces(trnz, n=nt, niter=10)
	sth, stt, stl = getsteps(FitX0)

	thxr = floor(thx, digits=2)
	println("threshold = $(thxr), best shot =$(bshot)")
	println("Step heights: ", sth)
	println("Step times: ", stt)
	println("Segment lengths: ", stl)
	plotfit2(dataX, FitX0, FitX1, S_curve, thxr)
end

# ╔═╡ 08dd1429-4c2d-4def-8282-e7abe469f318
md"""
### Fit data
"""

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

# ╔═╡ 1be93848-5757-48e6-8d5b-638cb11c4a61


# ╔═╡ 7f3d36a1-ac31-4cbb-b994-cbaf946d25fd
"""
    plot_traces_by_nmolx(df::DataFrame, molx_id::Int; lyt=(3, 3), size=(1200, 800))

Plot traces for all molecules with a given `nmolx` group in the provided DataFrame.

# Arguments
- `df`: DataFrame containing columns `nmolx` and `nmol`.
- `molx_id`: The `nmolx` group to select.
- `lyt`: Tuple specifying the subplot layout.
- `size`: Tuple specifying the plot size.
"""
function plot_traces_by_nmolx(df::DataFrame, molx_id::Int; lyt=(3, 3), size=(1200, 800))
    # Filter rows for molecules with the same nmolx
    mol_df = filter(:nmolx => ==(molx_id), df)
    # Get unique nmol values from this subset
    nmol_ids = unique(mol_df.nmol)
    # Plot the corresponding traces
    plot_traces(df, nmol_ids; lyt=lyt, size=size)
end

# ╔═╡ a3c2e539-eda7-4784-b90c-43fcf6982347
md"""
- Look for molecules with the same decay time
"""

# ╔═╡ 9b37fd0a-99eb-492a-a271-f584b043ef89
function unify_molecules_by_decay(df::DataFrame; time_column=:stepTime, id_column=:nmol, threshold=1.0)
    # Filter decay times greater than threshold
    df_valid = df[df[!, time_column] .> threshold, :]

    # Group by decay time and collect unique nmol per group
    grouped = groupby(df_valid, time_column)
    
    decay_map = Dict{Int, Int}()
    label = 1
    
    for g in grouped
        nmols = unique(g[!, id_column])
        for nm in nmols
            decay_map[nm] = get(decay_map, nm, label)
        end
        label += 1
    end

    # Create the new column nmolx
    df.nmolx = [get(decay_map, nm, nm) for nm in df[!, id_column]]
    
    return df
end

# ╔═╡ 83df15d6-5493-48ef-a8f4-29b8fc375da3
count_mol(df; cmol="nmol") = length(unique(df[:, cmol]))

# ╔═╡ 1f908e9b-855a-4ecc-a8d7-8ef92b502c54
md"""
## Functions
"""

# ╔═╡ baad60ca-3e4f-480c-a11a-c5b1d83e1eb6
function find_max(imst; tf=1)
	# Extract the frame at t 
	frame1 = imst[:, :, tf]
	
	# Find the maximum value and its linear index
	max_val = maximum(frame1)
	idx = argmax(frame1)
	
	# Convert linear index to (i, j)
	(i, j) = Tuple(CartesianIndices(frame1)[idx])
	i,j,max_val
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

# ╔═╡ 9893562b-800d-4384-aa91-d45e10de1cb6


# ╔═╡ b46f0c62-45f4-44ca-946c-bc9447306c45


# ╔═╡ 776c8d57-15c9-45e1-8ddb-c70afe4240a2



# ╔═╡ 4f558c30-f5ba-4b15-a72a-da16939c0bf5


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



# ╔═╡ 4059589a-bea7-4822-b8ad-d0eda5b7152c
begin
	DX3, FX3, TL3, xsdf =  fit_data2(trnz, niter=10, thr=0.0, rmax=120000, fplot=1)
	xsdf
end

# ╔═╡ ea5d4c44-c9cb-4603-8ad5-0d6070779e5d
plot_trx(DX3,FX3, nx=10, ny=5, figsize=(1500, 2500))

# ╔═╡ afd8b98c-aee7-44bb-96ac-3e7acd7b2f00
let
	htl, htlp = hist1d(xsdf.nstep .-1, "number of steps", 10, 0, 10, datap=false)
	htl
	xf0 = htl.weights[1]/sum(htl.weights)
	xf1 = sum(htl.weights[2:4])/sum(htl.weights)
	md"""
	- Fraction of zero steps =$(xf0)
	- Fraction of one/three steps =$(xf1)
	- number of mol trajectories $(length(unique(xsdf.nmol)))
	"""
end

# ╔═╡ 5fe515e5-23af-4e04-a1da-944c6f464062
plot_sdf(xsdf, nbins=(10,20,20,20))

# ╔═╡ 94562765-2067-4111-b343-206f404cc599
sdf1 = filter(row -> row.stepHeight > 0, xsdf)

# ╔═╡ 7ea511fb-c91a-4b2b-9914-3a4fb29540eb
begin
	sdf2 = filter(row -> row.nstep .-1 > 0, sdf1)
	renumber_nmol!(sdf2)
end

# ╔═╡ 52d7b937-1ffd-40cb-9dc5-6c25fcce2de0
plot_sdf(sdf2, nbins=(10,20,20,20))

# ╔═╡ 3e9fd7f8-56c0-4b3d-a422-7a8faa7a8245
plot_traces(sdf2, collect(1:47), lyt=(10,5), size=(1500, 1500))

# ╔═╡ 84023e48-6534-42a8-a4f4-23d242d2ca88
sdf2

# ╔═╡ e6da6f93-3123-4d52-ad47-8812847fbef8
msthdf = combine(groupby(sdf2, :nmol), :stepHeight => minimum => :min_stepHeight)

# ╔═╡ cc02c3cd-cc1f-4084-ac53-1f274571843b
begin
	hmsh = hist1d(msthdf.min_stepHeight, 100, 0.0, 200.0)
	plot(step_hist(msthdf.min_stepHeight;
              nbins=100,
              xlim=(0.0, 200.0),
              logy=false,
              xlabel="step height",
              ylabel=" #events "), size=(400,200))
end

# ╔═╡ 9aa87760-ac29-49a9-93a7-8993330e683b
begin
	hh = get_histo1d(hmsh)
	xf1 = sum(hh.weights[1])/sum(hh.weights)
	
	md"""
	- Fraction of step height in first three bins =$(xf1)
	- cutoff at $(hh.centers[3])
	"""
	
end

# ╔═╡ 5b4f4962-9118-4c85-b005-08c9f46324a2
begin
	thrl = 10
	sdf2.stepHeight .= ifelse.(sdf2.stepHeight .< thrl, 1.0, sdf2.stepHeight)
	#sdf3 = filter(row -> row.stepHeight >9, sdf2)
	sdf3 = transform(groupby(sdf2, :nmol), nrow => :nstep)
end

# ╔═╡ 12af547f-37c6-4954-b512-08b369b00966
plot_traces(sdf3, collect(1:45), lyt=(9,5), size=(1500, 1500))

# ╔═╡ c5edbc73-07ff-4f5e-ae3b-ee9dc4ffa98b
plot_sdf(sdf3, nbins=(10,20,20,20))

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
	
end

# ╔═╡ e92082d9-7b34-4d27-b46b-f7e197c1b2ad
mxsthdf = combine(groupby(sdf3, :nmol), :stepHeight => maximum => :max_stepHeight)

# ╔═╡ 309c89fa-a70d-43fa-adc2-6eacd10369f8
begin
	plot(step_hist(mxsthdf.max_stepHeight;
              nbins=10,
              xlim=(0.0, 3000.0),
              logy=false,
              xlabel="# max step height",
              ylabel=" #entries "), size=(400,200))
end

# ╔═╡ 883efc84-9357-4673-b80b-242999b7b56d
begin
	hmsx = hist1d(mxsthdf.max_stepHeight, 10, 0.0, 3000.0)
	hhx = get_histo1d(hmsx)
	x1 = hhx.weights[1]/sum(hhx.weights)
	md"""
	Fraction of height below $(hhx.edges[2]) is $(x1)
	"""
end

# ╔═╡ 2df3df4f-60b0-4c81-8c17-59b92b9d1a74
sdff = unify_molecules_by_decay(sdf3, threshold=1.0)

# ╔═╡ 76808ab6-6887-4fcb-92ec-92d0233c0943
count_mol(sdff, cmol="nmol")

# ╔═╡ 28cda987-0ad5-4a34-b3e1-06f08ea81c71
count_mol(sdff, cmol="nmolx")

# ╔═╡ 39aac501-6d96-4ab2-a5bc-ca8b2376a8a1
plot_traces_by_nmolx(sdff::DataFrame, 2; lyt=(4, 4), size=(1200, 800))

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

# ╔═╡ Cell order:
# ╠═9292025e-0d7d-11f0-365e-f1724fc39b3c
# ╠═5b3b5f83-058a-4582-9482-9b2fc5892692
# ╠═f8573e77-38d8-467c-a123-9a83f73e8970
# ╠═19476521-78bd-44d6-a302-6e971f9fc095
# ╠═cba7fc75-8363-4ffa-b5a6-6e7d34363813
# ╠═583a9aee-08eb-4f5a-95ef-d0087eb98cbc
# ╠═39b011c6-f511-42dd-befc-eaf3fd17ea1a
# ╠═b5f399bf-5713-4f26-afb0-2d5771dbbc6f
# ╠═11e43f25-3fa8-4f2d-bba5-1773a4989178
# ╠═b3b16805-b5b8-4782-a49c-15029b3a749d
# ╠═8e8cc911-5a3c-4ec5-9498-cd60a09fd276
# ╠═b7508873-d85d-4502-8845-b779feb4666c
# ╠═84678b4f-847b-4883-ab0d-3086468e30b8
# ╠═c8dfaaf4-1430-419d-b048-93dccc92885b
# ╠═77fa8f9e-2b6f-4743-a4ad-247c9e08ef54
# ╠═bfc6bdfa-9e79-4be2-b2f9-c264ea27cfb4
# ╠═cfbf4c54-78ee-46e6-ace2-0877e04543a5
# ╠═6b54a654-7850-4a27-8ec4-2dcbd3566d82
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═709ca0ac-aaba-440c-a3c7-39e8dd0fee2e
# ╠═de42a1a5-aa7e-49a6-98b8-d27b3ee16849
# ╠═3a8e7300-d209-4e50-ace8-c4b841d71f42
# ╠═a7542de4-5e78-4c55-91b7-d68c4e2f7985
# ╠═08dd1429-4c2d-4def-8282-e7abe469f318
# ╠═4059589a-bea7-4822-b8ad-d0eda5b7152c
# ╠═ea5d4c44-c9cb-4603-8ad5-0d6070779e5d
# ╠═afd8b98c-aee7-44bb-96ac-3e7acd7b2f00
# ╠═5fe515e5-23af-4e04-a1da-944c6f464062
# ╠═8e8fd656-17df-44d0-9633-474155d7dc7f
# ╠═94562765-2067-4111-b343-206f404cc599
# ╠═b40f0395-9252-4381-9ac9-32e7d9084a12
# ╠═7ea511fb-c91a-4b2b-9914-3a4fb29540eb
# ╠═52d7b937-1ffd-40cb-9dc5-6c25fcce2de0
# ╠═3e9fd7f8-56c0-4b3d-a422-7a8faa7a8245
# ╠═84023e48-6534-42a8-a4f4-23d242d2ca88
# ╠═614c7468-34a4-4d39-8820-4d6aef99bd83
# ╠═e6da6f93-3123-4d52-ad47-8812847fbef8
# ╠═cc02c3cd-cc1f-4084-ac53-1f274571843b
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═9aa87760-ac29-49a9-93a7-8993330e683b
# ╠═5b4f4962-9118-4c85-b005-08c9f46324a2
# ╠═12af547f-37c6-4954-b512-08b369b00966
# ╠═c5edbc73-07ff-4f5e-ae3b-ee9dc4ffa98b
# ╠═4014329c-e371-4ef8-a77e-738a668bf790
# ╠═e92082d9-7b34-4d27-b46b-f7e197c1b2ad
# ╠═309c89fa-a70d-43fa-adc2-6eacd10369f8
# ╠═883efc84-9357-4673-b80b-242999b7b56d
# ╠═7f3d36a1-ac31-4cbb-b994-cbaf946d25fd
# ╠═a3c2e539-eda7-4784-b90c-43fcf6982347
# ╠═9b37fd0a-99eb-492a-a271-f584b043ef89
# ╠═2df3df4f-60b0-4c81-8c17-59b92b9d1a74
# ╠═83df15d6-5493-48ef-a8f4-29b8fc375da3
# ╠═76808ab6-6887-4fcb-92ec-92d0233c0943
# ╠═28cda987-0ad5-4a34-b3e1-06f08ea81c71
# ╠═39aac501-6d96-4ab2-a5bc-ca8b2376a8a1
# ╠═1f908e9b-855a-4ecc-a8d7-8ef92b502c54
# ╠═baad60ca-3e4f-480c-a11a-c5b1d83e1eb6
# ╠═418f05eb-9b20-418f-9a35-783eac94e501
# ╠═9893562b-800d-4384-aa91-d45e10de1cb6
# ╠═b46f0c62-45f4-44ca-946c-bc9447306c45
# ╠═776c8d57-15c9-45e1-8ddb-c70afe4240a2
# ╠═4f558c30-f5ba-4b15-a72a-da16939c0bf5
# ╠═8b5448b2-73b2-4c9b-9f64-398a7a7a7720
# ╠═f41f04da-5666-4c28-97c4-16eba8e6e7cf
