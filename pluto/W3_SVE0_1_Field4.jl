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
	using LabStepAnalysis
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
	using Images, FileIO, ImageIO
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

# ╔═╡ a94ab132-2949-4292-94d3-46db64809749
names(LabStepAnalysis)

# ╔═╡ b5f399bf-5713-4f26-afb0-2d5771dbbc6f
names(JStepFinder)

# ╔═╡ 11e43f25-3fa8-4f2d-bba5-1773a4989178
names(histos)

# ╔═╡ b3b16805-b5b8-4782-a49c-15029b3a749d
names(BoldLab)

# ╔═╡ 8e8cc911-5a3c-4ec5-9498-cd60a09fd276
set_log_level(SimpleLogger.ERROR)

# ╔═╡ b7508873-d85d-4502-8845-b779feb4666c
debug("test")

# ╔═╡ e30919fe-95fe-47fb-b765-38730b57f526
begin
	root = "/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign"
	week ="Week3"
	day = "09_Wednesday"
	sample = "W3_SVE0_1"
	field = "Field4"
	md"""
	### Data description
	
	- 09/04/2025
	- power = 130 μW
	- 6 fields
	- Coverslip treated with O3 and then functionalized with NAPH3, solvent used was ACN Uvasol
	- Evaporated
	"""
		
end

# ╔═╡ 534904a0-92f2-435d-9a6d-15d7d0a53560
begin
	folder = joinpath(root, week, "Data", day, sample, field)
	
	readdir(folder)
end

# ╔═╡ a9ce63dd-5c84-4928-8f59-5907540624fb
begin
	nn = 1
md"""
- Reading frame $nn.
"""
end

# ╔═╡ 6d8c5e31-64bd-4431-b363-b3880f6479fe
md"""
- Regularize frame, to eliminate outlayers 
"""

# ╔═╡ 3e002e3b-d59a-499d-8460-44e3acb20bef
md"""
- ci gives the cartesian indices of spots above threshold
"""

# ╔═╡ 92d1a34c-9ca0-4f0f-a289-515b04f1e6d6
begin
	imxt = load_tif_stack_int16(folder, pedestal=0)
	regularize_stack!(imxt, nsigma=10)
	md"""
	- Load the full stack: this is 3D matrix $(typeof(imxt))
	"""
end

# ╔═╡ 2fe0ccbc-0b38-415c-995e-987017bcb499
imx1, CI = regularize_img(imxt[:,:,1])

# ╔═╡ 57212cd7-e45d-44a5-85f9-0da89007806f
imxt

# ╔═╡ af3e3cde-368d-49ad-a315-e83e9414a263
plot_frames(imxt; nscale=40)

# ╔═╡ 783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
begin
	totalI, meanI, stdI = get_stats(imxt; bck=0.0)
	plot_stats(totalI, meanI, stdI)
end

# ╔═╡ 57bd56be-5729-44e8-aba3-783a78c714d2
TRZ = build_traces(imxt; nsigma=10)

# ╔═╡ 6b54a654-7850-4a27-8ec4-2dcbd3566d82
#hsum2d, hmean2d = traces_h2d(xsum, xmean, xstd; bins=50)

# ╔═╡ 74deedfc-a344-4373-9fc9-db22e83d48ac
md"""
#### Get traces from the imst matrix
"""

# ╔═╡ 3a8e7300-d209-4e50-ace8-c4b841d71f42
md"""
#### Plot traces
"""

# ╔═╡ a7542de4-5e78-4c55-91b7-d68c4e2f7985
begin
	i = 1
	j = 2
	md"""
	- Fitting trace in pixel ($(i), $(j))
	- Intensity in pixel = $(sum(TRZ[i,j]))
	"""
end

# ╔═╡ d1c82dbe-8493-464d-bdac-f505657678d6
begin
	
	dataX, FitX, S_curve, bshot = fit_traces(TRZ, i, j; niter=10)
	sth, stt, stl = getsteps(FitX)

	md"""
	- Fit results
	- bes shot =$(bshot)
	- Step heights: = $(vect_to_fstr(sth, "%.2f"))
	- Step times =$(vect_to_fstr(stt, "%.2f"))
	- Segment lengths = $(vect_to_fstr(stl, "%.2f"))
	"""
	
end

# ╔═╡ 6e2b4709-a8df-47ad-9567-bcf39e020bc8
size(TRZ)

# ╔═╡ 15a5e0fc-9416-455e-ba93-fb8ea97a75fe
plotfit(dataX, FitX, S_curve, i, j)

# ╔═╡ 2b424306-f327-4c7f-9e50-0a4cf9d44de9
typeof(FitX)

# ╔═╡ 193ee1a5-3803-4f3a-9c4c-57ddee590c25
typeof(dataX)

# ╔═╡ 08dd1429-4c2d-4def-8282-e7abe469f318
md"""
### Fit data
"""

# ╔═╡ 85a3b5cd-9130-475f-acdf-c0eb20ee8a6a
begin
	DX, FX, TL, sdf, npixel, nfailed = fit_data3(TRZ; ped=1600.0, tt = 0, niter=10, thr=0.0, fplot=2500)
	md"""
	- Fitted $(npixel) pixels 
	- failed fit fraction =$(nfailed/npixel)
	"""
end

# ╔═╡ 63518d63-7480-42a0-875a-7d5834e28f1b



# ╔═╡ ea5d4c44-c9cb-4603-8ad5-0d6070779e5d
plot_trx(DX,FX, nx=10, ny=5, figsize=(1500, 2500))

# ╔═╡ afd8b98c-aee7-44bb-96ac-3e7acd7b2f00
let
	htl = hist1d(sdf.nstep .-1, 10)
	xf0 = htl.weights[1]/sum(htl.weights)
	xf1 = sum(htl.weights[2:4])/sum(htl.weights)
	md"""
	- Fraction of zero steps =$(xf0)
	- Fraction of one/three steps =$(xf1)
	- number of mol trajectories $(length(unique(sdf.nmol)))
	"""
end

# ╔═╡ aaa4feca-d839-4736-bfae-8eec83bd9c05
plot_sdf(sdf::DataFrame; nbins=(10,100,100,100),  
		 xlim=((0.0,10.0),(-10.0, 500.0), (0.0, 400.0), (0.0, 400.0)))

# ╔═╡ 3e9fd7f8-56c0-4b3d-a422-7a8faa7a8245
plot_traces(sdf, collect(1:9), lyt=(3,3), size=(1500, 1500))

# ╔═╡ 7f415bac-05d6-4771-8f4b-4bc552fadf4d
sdf

# ╔═╡ 84023e48-6534-42a8-a4f4-23d242d2ca88
sdf2 = sdf[sdf.stepLength .> 10, :]

# ╔═╡ 47064a15-9ad8-4250-9608-6f1165187db3
plot_traces(sdf2, collect(1:9), lyt=(3,3), size=(1500, 1500))

# ╔═╡ 820beb62-fe16-494b-8e76-323ce30e46d0
#update_nstep!(sdf2)

# ╔═╡ 3b85b24d-34bb-4df7-a0c3-cb89f7d90b84
begin

vnm = sdf2.nmol 

# Get unique values
unique_nm = unique(vnm)

# Get count of each element
counts_dict = countmap(vnm)

# Convert dictionary to count vector in the same order as unique_vals
vns = [counts_dict[val] for val in unique_nm]
	md"""
	- Data set has $(length(unique_nm)) molecules
	"""
end

# ╔═╡ 867ab7f9-c5d3-4840-9773-eeb8c33183db
plot_sdf(sdf2::DataFrame; nbins=(10,50,50,50),  
		 xlim=((0.0,10.0),(-10.0, 500.0), (0.0, 400.0), (0.0, 400.0)),
		 log=(false,false,false,false)
		)

# ╔═╡ cd6db8a7-c328-4453-94da-14d9678a4c66
begin
	plot(step_hist(vns;
              nbins=10,
              xlim=(0.0, 10.0),
              logy=false,
              xlabel="# steps",
              ylabel=" #entries "), size=(400,200))
end

# ╔═╡ 614c7468-34a4-4d39-8820-4d6aef99bd83
md"""
- Compute min step height for each story. This is done to define a baseline
"""

# ╔═╡ e6da6f93-3123-4d52-ad47-8812847fbef8
msthdf = combine(groupby(sdf2, :nmol), :stepHeight => minimum => :min_stepHeight)

# ╔═╡ 309c89fa-a70d-43fa-adc2-6eacd10369f8
begin
	plot(step_hist(msthdf.min_stepHeight;
              nbins=100,
              xlim=(0.0, 300.0),
              logy=false,
              xlabel="# min step height",
              ylabel=" #entries "), size=(400,200))
end

# ╔═╡ f6100b17-a7f0-4a03-9e18-a073619b4db3
mxsthdf = combine(groupby(sdf2, :nmol), :stepHeight => maximum => :max_stepHeight)

# ╔═╡ 9fba1c88-0b18-4c1f-a687-c4b452a8adac
begin
	plot(step_hist(mxsthdf.max_stepHeight;
              nbins=100,
              xlim=(10.0, 300.0),
              logy=false,
              xlabel="# max step height",
              ylabel=" #entries "), size=(400,200))
end

# ╔═╡ 81181098-822a-48f9-b776-680735de6430
begin
	hsth = hist1d(sdf2.stepHeight, 100)
	typeof(hsth)
end

# ╔═╡ 1be93848-5757-48e6-8d5b-638cb11c4a61
md"""
## Functions
"""

# ╔═╡ f0aab374-2e76-47e0-a62c-a5ec21129767


# ╔═╡ 1b968f57-798f-471a-9203-437f373cb65f
function update_nstep!(df::DataFrame)
    @assert "nmol" in names(df) "Missing column 'nmol'"
    @assert "nstep" in names(df) "Missing column 'nstep'"

    gdf = groupby(df, :nmol)
    for g in gdf
        g[!, :nstep] .= nrow(g)  # assign group size to all rows in group
    end
    return df
end

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

# ╔═╡ baad60ca-3e4f-480c-a11a-c5b1d83e1eb6
begin
	"""
	Find max of a frame in a stack and return value of the maximum and indices
	"""
	function find_max(imst::AbstractArray{T, 3}; tf=1) where {T<:Real}
		find_max(imst[:, :, tf])
	end
	
	
	"""
	Find max of a matrix and return value of the maximum and indices
	"""
	function find_max(frame::AbstractArray{T, 2}) where {T<:Real}
		
		# Find the maximum value and its linear index
		max_val = maximum(frame)
		idx = argmax(frame)
		
		# Convert linear index to (i, j)
		(i, j) = Tuple(CartesianIndices(frame)[idx])
		i,j,max_val
	end
	
	
	"""
	Find max of a vector and return maximum and index
	"""
	function find_max(frame::AbstractVector{T}) where {T<:Real}
		
		# Find the maximum value and its linear index
		max_val = maximum(frame)
		idx = argmax(frame)
	    idx, max_val
	end
end

# ╔═╡ 4b7e4552-ffb9-4f9a-a25f-f45178b2fbd3
function stats(sdf, vns)
	hnst = hist1d(vns, 10)
	xf1 = hnst.weights[1]/sum(hnst.weights)
	xf2 = sum(hnst.weights[1:2])/sum(hnst.weights)
	xf3 = sum(hnst.weights[1:3])/sum(hnst.weights)
	xf4 = sum(hnst.weights[4:end])/sum(hnst.weights)
	xmnst = sum(hnst.weights .* hnst.centers) ./sum(hnst.weights)

	hsth = hist1d(sdf.stepHeight, 100)
	idx, max_val = find_max(hsth.weights)
	mxsth = hsth.centers[idx]

	hstl = hist1d(sdf.stepLength, 100)
	idx, max_val = find_max(hstl.weights)
	mxstl = hstl.centers[idx]
	xf1, xf2, xf3, xf4, xmnst, mxsth, mxstl	
	
	
end

# ╔═╡ 2d9edd87-8ff1-4d90-a440-4a702a8c1144
begin
	xf1, xf2, xf3, xf4, xmnst, mxsth, mxstl = stats(sdf2, vns)
	md"""
	#### Statistics
	- sample = $(sample)
	- field = $(field)
	- Fitted $(npixel) pixels 
	- Pixel fraction with no signal =$(nfailed/npixel)
	- Fraction with 0 steps = $(xf1)
	- Fraction with (0,1) steps = $(xf2)
	- Fraction with (0,1,2) steps = $(xf3)
	- Fraction with >3  steps = $(xf4)
	- Mean number of steps =$(xmnst)
	- Max value of step height = $(mxsth)
	- Max value of step length = $(mxstl)
	"""
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

# ╔═╡ Cell order:
# ╠═9292025e-0d7d-11f0-365e-f1724fc39b3c
# ╠═5b3b5f83-058a-4582-9482-9b2fc5892692
# ╠═f8573e77-38d8-467c-a123-9a83f73e8970
# ╠═19476521-78bd-44d6-a302-6e971f9fc095
# ╠═cba7fc75-8363-4ffa-b5a6-6e7d34363813
# ╠═583a9aee-08eb-4f5a-95ef-d0087eb98cbc
# ╠═39b011c6-f511-42dd-befc-eaf3fd17ea1a
# ╠═a94ab132-2949-4292-94d3-46db64809749
# ╠═b5f399bf-5713-4f26-afb0-2d5771dbbc6f
# ╠═11e43f25-3fa8-4f2d-bba5-1773a4989178
# ╠═b3b16805-b5b8-4782-a49c-15029b3a749d
# ╠═8e8cc911-5a3c-4ec5-9498-cd60a09fd276
# ╠═b7508873-d85d-4502-8845-b779feb4666c
# ╠═e30919fe-95fe-47fb-b765-38730b57f526
# ╠═534904a0-92f2-435d-9a6d-15d7d0a53560
# ╠═a9ce63dd-5c84-4928-8f59-5907540624fb
# ╠═6d8c5e31-64bd-4431-b363-b3880f6479fe
# ╠═3e002e3b-d59a-499d-8460-44e3acb20bef
# ╠═92d1a34c-9ca0-4f0f-a289-515b04f1e6d6
# ╠═2fe0ccbc-0b38-415c-995e-987017bcb499
# ╠═57212cd7-e45d-44a5-85f9-0da89007806f
# ╠═af3e3cde-368d-49ad-a315-e83e9414a263
# ╠═783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
# ╠═57bd56be-5729-44e8-aba3-783a78c714d2
# ╠═6b54a654-7850-4a27-8ec4-2dcbd3566d82
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═3a8e7300-d209-4e50-ace8-c4b841d71f42
# ╠═a7542de4-5e78-4c55-91b7-d68c4e2f7985
# ╠═d1c82dbe-8493-464d-bdac-f505657678d6
# ╠═6e2b4709-a8df-47ad-9567-bcf39e020bc8
# ╠═15a5e0fc-9416-455e-ba93-fb8ea97a75fe
# ╠═2b424306-f327-4c7f-9e50-0a4cf9d44de9
# ╠═193ee1a5-3803-4f3a-9c4c-57ddee590c25
# ╠═08dd1429-4c2d-4def-8282-e7abe469f318
# ╠═85a3b5cd-9130-475f-acdf-c0eb20ee8a6a
# ╠═63518d63-7480-42a0-875a-7d5834e28f1b
# ╠═ea5d4c44-c9cb-4603-8ad5-0d6070779e5d
# ╠═afd8b98c-aee7-44bb-96ac-3e7acd7b2f00
# ╠═aaa4feca-d839-4736-bfae-8eec83bd9c05
# ╠═3e9fd7f8-56c0-4b3d-a422-7a8faa7a8245
# ╠═7f415bac-05d6-4771-8f4b-4bc552fadf4d
# ╠═84023e48-6534-42a8-a4f4-23d242d2ca88
# ╠═47064a15-9ad8-4250-9608-6f1165187db3
# ╠═820beb62-fe16-494b-8e76-323ce30e46d0
# ╠═3b85b24d-34bb-4df7-a0c3-cb89f7d90b84
# ╠═867ab7f9-c5d3-4840-9773-eeb8c33183db
# ╠═cd6db8a7-c328-4453-94da-14d9678a4c66
# ╠═614c7468-34a4-4d39-8820-4d6aef99bd83
# ╠═e6da6f93-3123-4d52-ad47-8812847fbef8
# ╠═309c89fa-a70d-43fa-adc2-6eacd10369f8
# ╠═f6100b17-a7f0-4a03-9e18-a073619b4db3
# ╠═9fba1c88-0b18-4c1f-a687-c4b452a8adac
# ╠═81181098-822a-48f9-b776-680735de6430
# ╠═2d9edd87-8ff1-4d90-a440-4a702a8c1144
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═4b7e4552-ffb9-4f9a-a25f-f45178b2fbd3
# ╠═f0aab374-2e76-47e0-a62c-a5ec21129767
# ╠═1b968f57-798f-471a-9203-437f373cb65f
# ╠═7f3d36a1-ac31-4cbb-b994-cbaf946d25fd
# ╠═a3c2e539-eda7-4784-b90c-43fcf6982347
# ╠═9b37fd0a-99eb-492a-a271-f584b043ef89
# ╠═83df15d6-5493-48ef-a8f4-29b8fc375da3
# ╠═baad60ca-3e4f-480c-a11a-c5b1d83e1eb6
# ╠═418f05eb-9b20-418f-9a35-783eac94e501
