### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

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
	using StepPreprocessing 
	using histos
	import Measures
	using NPZ
	using Unitful
end

# ╔═╡ 5b3b5f83-058a-4582-9482-9b2fc5892692
begin
	using PlutoUI
	using Observables
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
	using Images, FileIO, ImageIO, ImageView, ImageSegmentation, ImageMorphology
	using SparseArrays
end

# ╔═╡ b9ef4450-d619-4e67-9ff3-9c8045b3863d
using LinearAlgebra

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

# ╔═╡ 7142e579-224c-474d-966f-461f8ce82e3a
names(StepPreprocessing)

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

# ╔═╡ f724b5fc-b506-4442-b6b0-92b8fc9ad16b
	list_subfolders(dir) = [basename(e) for e in readdir(dir; join=true) if isdir(e) && !startswith(basename(e), ".")]


# ╔═╡ 55941a45-56f4-48c8-b823-880bdecaca29
md"""
### Select Data
"""

# ╔═╡ 5d65761c-c47f-49d2-b4e8-3b758bc96e8b
md"""
2. Top level data set maps
"""

# ╔═╡ e76e1ca2-c736-42e7-8ab7-1f7bd09086e0
begin
    BASEDIRS = Dict("MonteCarlo" => "/Users/jjgomezcadenas/Projects/BoldLab/pluto/npy" ,"Data"=>"/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign"				    
    )
    base_directory(label) = BASEDIRS[label]
end

# ╔═╡ fe676ab1-4ec5-4c3b-beb1-62c68f4f87ee
begin
  
md"""
3. Folder scanner
"""
end

# ╔═╡ 53555a56-6602-4eb1-84a7-a64f4d8cb41e
md"""
4. Choose Data/MC
"""

# ╔═╡ abb8b3bf-6d9f-4505-b13f-34de69460c51
@bind dataset_label Select(collect(keys(BASEDIRS)))


# ╔═╡ 5f931b55-8578-4958-b6dc-24f4cfb3c011
md"""
7. Change root dir if needed
"""

# ╔═╡ d155780e-74ed-48f4-98d9-45f23ff20e3f
begin
	root_dir = base_directory(dataset_label)
	@bind path_input TextField(120, root_dir)
end

# ╔═╡ d4d1fcb8-ff53-4f49-bd5e-47e12672d8dc
md"""
10. Folder picker
"""

# ╔═╡ 8e1005ba-ecff-4b97-be2c-3683529fa3c5
begin
   subdirs, npys, _ = scan_level(root_dir)

    if !isempty(npys)
		casemc = true
		casedata = false
		@bind file_mc Select(npys) 
	else
		casemc=false
		casedata = true
		md"""
		- Select week 
		"""
	end

end


# ╔═╡ 096d3fd6-38e6-4da8-8a81-6339c3624ac5
if casemc ==true
	pedx = 0.0
	path_mc = joinpath(root_dir, file_mc)
else
	pedx=1600.0
	@bind folder_week Select(subdirs) 
	
end

# ╔═╡ 802b9639-4766-4209-899a-5a6fc1cf2d24
if casedata == true	
	path_week = joinpath(root_dir, folder_week, "Data")
	subdirs2, _, _ = scan_level(path_week)
	@bind folder_day Select(subdirs2) 
end

# ╔═╡ 6c9cb955-f2b3-41fb-85af-52873837e127
if casedata == true
	path_day = joinpath(path_week, folder_day)
	subdirs3, _, _ = scan_level(path_day)
	@bind folder_scan Select(subdirs3) 
end

# ╔═╡ 49c7448d-824e-4604-9c90-c28e45f232d4
if casedata == true
	path_scan = joinpath(path_day, folder_scan)
	subdirs4, _, _ = scan_level(path_scan)
	@bind folder_field Select(subdirs4) 
end

# ╔═╡ 7732c78d-25be-4648-a386-4d455b8bd0d5
if casedata == true
	path_tiff = joinpath(path_scan, folder_field)
end

# ╔═╡ 92d1a34c-9ca0-4f0f-a289-515b04f1e6d6
if casedata
	nsigma = 5
	imxt = load_tif_stack_int16(path_tiff, pedestal=0.0)
	mask = regularize_stack!(imxt, nsigma=nsigma)
	md"""
	- Reading Data: Load the stack:  $(length(mask[mask .>0])) pixels masked with signal > $(nsigma) times the mean of the first frame
	"""
else
	imxt = npzread(path_mc)
	md"""
	- Reading MC: Image size = $(size(imxt))
	"""
end

# ╔═╡ c7cc7bfa-2ae3-44ef-8b98-5e5afbefcf2e
path_mc

# ╔═╡ fa4b76f1-f3e1-49ef-9054-404bfa942786
path_mc_tif = replace(path_mc, "/npy/" => "/tif/", r"\.npy$" => ".tif")

# ╔═╡ 2fe0ccbc-0b38-415c-995e-987017bcb499
#img = load(path_mc_tif)


# ╔═╡ b93360f2-5890-4138-b7e8-53a9d6f539e4
#imshow(img[:,:,1])

# ╔═╡ 968648d8-54f2-4485-8209-8c22d4b63a8a
#img2 = load(path_mc_tif) |> channelview  # ensures 2D array

# ╔═╡ d096680e-a09c-42cf-998e-b3304c509932
#plot_frames(img2; nscale=20)

# ╔═╡ 11a04a68-91d0-4595-9296-144b3e21f478
begin
	σ = 10.0 #average over large number of pixels 
	nlf = 5  # number of frames to average
	background = compute_background_from_stack(imxt; σ= σ, nlf=nlf)
	histogram(vec(background))
end

# ╔═╡ d2a8345c-9c44-43d7-a9ec-6ee4ab58674d
typeof(background)

# ╔═╡ b0a4c0be-12b9-416c-a45c-5fe35fbbd168
begin
	imgbsub = subtract_background_from_stack(imxt, background)
	vim = vec(imgbsub)
	histogram(vim[vim.>50])
end

# ╔═╡ 89354633-e220-45f2-989e-c5575acd2988
plot_frames(imgbsub; nscale=20)

# ╔═╡ f4ff7938-f16e-4159-95e8-cb98c59a9d80
begin
	σ_denoise = 1.0
	#k = Kernel.gaussian(σ_denoise)

# Apply Gaussian denoising to each frame
	#img_denoised_stack = imfilter.(eachslice(imgbsub, dims=3), Ref(k))
	#imgd = cat(img_denoised_stack...; dims=3)
	imgd = denoise_stack(imgbsub; σ=σ_denoise)
	plot_frames(imgd, nscale=20)
end

# ╔═╡ fd39234f-6c62-42e3-bde9-9a9e565fa519
begin
	vi2m = vec(imgd[:,:,2])
	histogram(vi2m[vi2m.>30])
end

# ╔═╡ f9a54f94-b109-4ae0-abe1-f44fdedc0d25
begin
	threshold = 30.0  
	min_area = 10
end

# ╔═╡ 389be056-e3d5-45fc-b80c-19f7cd4d5433
begin
	filtered_masks, region_stats = filter_regions_in_stack(imgd; i_thr=15, min_area=10)
	mask_stack = cat(filtered_masks...; dims=3)
	t = 1  # any frame
	frame = imgd[:, :, t]
	centroids = [r[:centroid] for r in region_stats[t]]
	
	heatmap(frame, color=:grays, aspect_ratio=1)
	scatter!([c[1] for c in centroids], [c[2] for c in centroids], markersize=4, color=:red)
end

# ╔═╡ e7cb1f63-130c-4e75-af5d-779fc1a4c755
"""
    detect_local_maxima(frame::AbstractMatrix{<:Real}; 
                        threshold::Real=0.0, 
                        dx::Int=0, dy::Int=0) 
        -> DataFrame

Detect local maxima in a 2D image, excluding a border of `dx` and `dy` pixels from the search area.

# Arguments
- `frame`: 2D image matrix (e.g., one frame from an image stack).
- `threshold`: Minimum intensity a peak must exceed (default: 0.0).
- `dx`: Margin to exclude on the left and right edges (columns).
- `dy`: Margin to exclude on the top and bottom edges (rows).

# Returns
- A `DataFrame` with columns:
    - `i`: row index of each peak (vertical coordinate)
    - `j`: column index (horizontal)
    - `intensity`: value at the peak

# Notes
- Padding uses `Pad(:replicate)` to preserve edge values, but maxima near the border are excluded.
"""
function detect_local_maxima(frame::AbstractMatrix{<:Real}; threshold=0.0, dx=0, dy=0)
    is_max = mapwindow(x -> x[5] == maximum(x), frame, (3, 3); 
					   border=Pad(:replicate))
    candidates = findall(is_max .& (frame .> threshold))

    i_vals = Int[]
    j_vals = Int[]
    intensities = Float64[]

    for I in candidates
        i, j = Tuple(I)

        # Skip peaks near the edge
        if i ≤ dy || j ≤ dx || i > size(frame, 1) - dy || j > size(frame, 2) - dx
            continue
        end

        push!(i_vals, i)
        push!(j_vals, j)
        push!(intensities, frame[I])
    end

    return DataFrame(i = i_vals, j = j_vals, intensity = intensities)
end

# ╔═╡ 62da10f6-0d5c-41c0-a985-d15c946f5b84
begin
    @bind show_peaks CheckBox(true)
end

# ╔═╡ 12960d51-135c-47cd-ab86-f2ab5bacef08
@bind nframe NumberField(0:199, default=1)

# ╔═╡ d34b9105-09d4-4876-842d-bcf74249cca9
@bind pthr NumberField(0:0.1:200.0, default=30.0)

# ╔═╡ 4d0c63a1-89e0-4d89-846c-e1501dbc2696
begin
	peaks = detect_local_maxima(imgd[:, :, nframe]; threshold=pthr, dx=0, dy=0)
	md"""
	- found $(size(peaks)[1]) molecule candidates 
	- in frame $(nframe) 
	- with thr = $(pthr)
	"""
end

# ╔═╡ 60b16840-2530-41df-bbc1-12b1066e9829
size(peaks)

# ╔═╡ 7927a18c-6942-42bd-ac6b-2ff720b155d0
histogram(peaks.intensity, bins=20)

# ╔═╡ 80068e3e-e9c4-47c5-b6cf-e10ff2de19ea
begin
	# Plot the grayscale image
	p1 = heatmap(imgd[:, :, nframe]; 
			color = cgrad(:grays, rev = true),
			#color = :grays, 
			aspect_ratio = 1, title = "Detected Peaks")

	if show_peaks
		scatter!(p1, 
		    peaks.j,               # x-axis (columns)
		    peaks.i,               # y-axis (rows)
		    zcolor = peaks.intensity,
		    colorbar = true,
		    marker = (:circle, 3),
		    c = :viridis,
		    label = "Peaks"
		)
	else
		p1
	end
end

# ╔═╡ 783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
begin
	totalI, meanI, stdI = get_stats(imxt; bck=0.0)
	plot_stats(totalI, meanI, stdI)
end

# ╔═╡ 57692746-065d-4d75-8a41-7ffafd69550e
md"""
- total intensity = $(Float64(sum(totalI)))
"""

# ╔═╡ 74deedfc-a344-4373-9fc9-db22e83d48ac
md"""
#### Get traces from the imst matrix
"""

# ╔═╡ 1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
TRZS = build_sparse_traces(imxt, peaks)

# ╔═╡ d92c4660-4950-4adf-aceb-bc94663411c6
function select_trace(TRZS, df::DataFrame, row::Int) 
	i = df.i[row]
	j = df.j[row]
    trace = TRZS[i, j]  # safe, only defined where needed
	i,j, trace
end

# ╔═╡ 3a8e7300-d209-4e50-ace8-c4b841d71f42
md"""
#### Plot traces
"""

# ╔═╡ 291c9e90-a915-4b35-a134-745ef253a72a
function plot_traces(TRZS, peaks; ftrz=1, ltrz=9,  figsize=(1500,1500))
	function pltd(tz, i, j)
		plot(1:length(tz), tz, 
		label="Trace =($(i),$(i))", color=:gray, lw=1,
		xlabel="time steps", ylabel="Intensity", title="", 
		legend=:topright, grid=true)
	end
	PP =[]
	ntrz = ltrz - ftrz + 1
	
	if ntrz <= 9 
		ly = (3,3)
	elseif ntrz <= 16 
		ly = (4,4)
	elseif  ntrz <= 25 
		ly = (5,5)
	end

	for it in 1:ntrz
		i,j, tz = select_trace(TRZS, peaks, it)
		
		push!(PP, pltd(tz, i, j))
	end
	plot(PP..., layout=ly, size=figsize)
end

# ╔═╡ af233a9a-2a3b-4b23-a982-c76d4a4c16f2
plot_traces(TRZS, peaks; ftrz=1, ltrz=25,  figsize=(1500,1500))

# ╔═╡ 7360e3df-c583-4894-86a6-5654b50a389c
plot_traces(TRZS, peaks; ftrz=26, ltrz=50,  figsize=(1500,1500))

# ╔═╡ 1c6507e3-821f-4e07-b41e-b3c106df3671
@bind ntrz NumberField(0:100, default=1)

# ╔═╡ 6f3a242a-6dff-432a-b30e-1b7ee1d234fc
i,j, tz = select_trace(TRZS, peaks, ntrz)

# ╔═╡ e04c9c53-7ae2-474b-87bc-97bd58f458fa
@bind thr NumberField(0.0:0.1:1.1, default=0.5)

# ╔═╡ 86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
@bind niter NumberField(0:10, default=3)

# ╔═╡ d1c82dbe-8493-464d-bdac-f505657678d6
begin
	dataX, FitX, S_curve, best_shot, iter, cc = fit_traces(tz, niter=niter, tresH=thr, sel="core")
	sth, stt, stl = getsteps(FitX)

	#- Step heights: = $(vect_to_fstr(sth, "%.2f"))
	#- Step times =$(vect_to_fstr(stt, "%.2f"))
	#- Segment lengths = $(vect_to_fstr(stl, "%.2f"))	

	md"""
	- Fit results
	- bes shot =$(best_shot)
	
	"""
	
end

# ╔═╡ 8aa50fbb-74a9-424d-9a84-e717247c65a9
plotfit(dataX, FitX, S_curve, i,j )

# ╔═╡ 1b250af1-ffe5-488c-9eee-c1ca41085901
typeof(S_curve)

# ╔═╡ 77f5c1c3-11bb-4976-af44-1b488f08de6b
md"""
## Fits to data
"""

# ╔═╡ 3115a147-d41b-4ab6-9ad9-0f3b30fb4504
md"""
- Set the vaue of threshold
"""

# ╔═╡ f175508e-f0ea-48a9-ba72-d2e21215de1d
md"""
- Set the vaue of iterations (number of steps to be sought in the data)
"""

# ╔═╡ 5042becb-29c6-447e-84ad-a965a9961992
begin
	dfs, zI, zJ, zDX, zFX, zSC  =find_fit_candidates2(TRZS, peaks;  sel="core", ped=0.0, niter=niter, thr=thr)
	dfs
end

# ╔═╡ 01391609-5034-4dbd-8319-10ac82126cfc
length(zDX)

# ╔═╡ 264fcbb0-ec40-42b3-ab05-01b2171da6f2
begin
	md"""
	- number of fitted molecules = $(length(unique(dfs.nmol)))
	- threshold = $(thr)
	"""
end

# ╔═╡ 1be93848-5757-48e6-8d5b-638cb11c4a61
md"""
## Functions
"""

# ╔═╡ efd033eb-bde9-4f4e-85f0-777125555edd
function plotsc(S_curve, II, JJ)
	plt2 = plot(1:length(S_curve), S_curve, 
	marker=:circle, label="S-curve",
	xlabel="Iteration", ylabel="S-value (MSE ratio)", 
	title="Goodness of Fit (S-curve), Trace =($(II),$(JJ))", 
	grid=true, legend=:topright)
end

# ╔═╡ d9477e2f-2bad-4243-a11f-393f5b3405c7
function plot_sc(VSC, VI,VJ; plotsel="3x3", figsize=(1500,1500))
	PP =[]
	jl = 9
	ly = (3,3)
	
	if plotsel == "4x4"
		jl = 16
		ly = (4,4)
	elseif plotsel == "5x5"
		jl = 25
		ly = (5,5)
	end

	jl = min(jl, length(VSC))
	for i in 1:jl
		push!(PP, plotsc(VSC[i], VI[i], VJ[i]))
	end
	plot(PP..., layout=ly, size=figsize)
#plotfit(dataX, FitX, S_curve, i, j)
end

# ╔═╡ 64ed321d-3050-465e-8c4f-fea826265779
plot_sc(zSC, zI,zJ; plotsel="4x4", figsize=(1500,1500))

# ╔═╡ fa0da6ef-68ba-4db2-aec0-26fd3a280b5a
function pltf(dataX,  FitX1, II, JJ)  
	
	plt1 = plot(1:length(dataX), dataX, 
	label="Noisy Signal", color=:gray, lw=1,
	xlabel="time steps", ylabel="Intensity", title="Trace =($(II),$(JJ))", 
	legend=:topright, grid=true)

	plot!(plt1, 1:length(FitX1), FitX1, 
	label="Fit1", color=:blue, lw=2)
	plt1
end

# ╔═╡ 862d9a0c-ecc6-4178-80f4-74fc71a79f38
function plot_fits(VDX, VFX, VI, VJ; plotsel="3x3", figsize=(1500,1500))
	PP =[]
	jl = 9
	ly = (3,3)
	
	if plotsel == "4x4"
		jl = 16
		ly = (4,4)
	elseif plotsel == "5x5"
		jl = 25
		ly = (5,5)
	end

	jl = min(jl, length(VDX))
	for i in 1:jl
		push!(PP, pltf(VDX[i], VFX[i], VI[i], VJ[i]))
	end
	plot(PP..., layout=ly, size=figsize)
#plotfit(dataX, FitX, S_curve, i, j)
end

# ╔═╡ 1015d8f6-d32c-4920-b688-ab1848cd9c5a
plot_fits(zDX, zFX, zI, zJ; plotsel="4x4")

# ╔═╡ 57100463-320e-4768-9249-4d38b5b4f6b4


# ╔═╡ d4a11884-b670-4ceb-9444-36d3bcc6e8a7
function plotsdf(sdf)

		
	# take only one row per (i,j), keeping the first nstep
	df_pix = unique(sdf, [:i, :j, :nstep])

	
	hnstep, pnstep = step_hist(df_pix.nstep;
              nbins=10,
              xlim=(0.0, 10.0),
              logy=false,
              xlabel=" # steps",
              ylabel=" #entries ")

	hbshot, pbshot = step_hist(df_pix.bestShot;
              nbins=10,
              xlim=(0.0, 10.0),
              logy=false,
              xlabel=" # best shot",
              ylabel=" #entries ")
	
	hstptime, pstptime = step_hist(sdf.stepTime;
              nbins=50,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" step time",
              ylabel=" #entries ")

	hstph, pstph = step_hist(sdf.stepHeight;
              nbins=50,
              xlim=(0.0, 500.0),
              logy=false,
              xlabel=" # step height",
              ylabel=" #entries ")

	hstphx, pstphx = step_hist(df_pix.stepHeightMax;
              nbins=50,
              xlim=(0.0, 3500.0),
              logy=false,
              xlabel=" # step height max",
              ylabel=" #entries ")

	hstplx, pstplx = step_hist(sdf.stepLength;
              nbins=50,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" # step length",
              ylabel=" #entries ")
	
	plot(pstphx, size=(400,200))
	
	plot(pbshot, pnstep, pstptime, pstph, pstphx, pstplx, 
		 layout=(3, 2), size=(800,600))
	
	
end

# ╔═╡ 0a8baaba-06a6-4d44-bbc7-437e740be7b2
plotsdf(dfs)

# ╔═╡ 4b7e4552-ffb9-4f9a-a25f-f45178b2fbd3
function stats(sdf,)
	df_pix = unique(sdf, [:i, :j, :nstep])
	hnst, _ = step_hist(df_pix.nstep;
              nbins=20,
              xlim=(0.0, 10.0),
              logy=false,
              xlabel=" # steps",
              ylabel=" #entries ")
	

	hsth, _ = step_hist(sdf.stepHeight;
              nbins=100,
              xlim=(0.0, 500.0),
              logy=false,
              xlabel=" # step height",
              ylabel=" #entries ")

	hstl, _ = step_hist(sdf.stepLength;
              nbins=100,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" # step length",
              ylabel=" #entries ")
	
	xmnst = sum(hnst.weights .* hnst.centers) ./sum(hnst.weights)

	hsth = hist1d(sdf.stepHeight, 100)
	idx, max_sth = find_max(hsth.weights)
	mxsth = hsth.centers[idx]

	idx, max_stl = find_max(hstl.weights)
	mxstl = hstl.centers[idx]
	hnst, xmnst, mxsth, max_sth, mxstl,max_stl 	
	
	
end

# ╔═╡ 81181098-822a-48f9-b776-680735de6430
function save_stats(sdf2, npixels, sample, field)
	hnst, xmnst, mxsth, max_sth, mxstl, max_stl = stats(sdf2)

	# Construct output string
	outstr = """
	#### Statistics
	- sample = $(sample)
	- field = $(field)
	- total intensity = $(Float64(sum(totalI)))
	- Fitted $(npixels) pixels 
	
	- #### steps:
	- edges = $(vect_to_fstr(hnst.edges, "%.2f"))
	- weights = $(vect_to_fstr(hnst.weights, "%.2f"))
	- Mean number of steps = $(xmnst)
	- Max of step height at = $(mxsth)
	- Max of step height value = $(max_sth)
	- Max of step length at = $(mxstl)
	- Max of step length value = $(max_stl)
	"""
	file =string("step_stats_sample_", sample, "_field_", field, ".md")
	# Save to file
	write(file, outstr)

	# Return markdown for Pluto display
	outstr
end

# ╔═╡ f0aab374-2e76-47e0-a62c-a5ec21129767


# ╔═╡ a3c2e539-eda7-4784-b90c-43fcf6982347
md"""
- Look for molecules with the same decay time
"""

# ╔═╡ 83df15d6-5493-48ef-a8f4-29b8fc375da3
count_mol(df; cmol="nmol") = length(unique(df[:, cmol]))

# ╔═╡ Cell order:
# ╠═9292025e-0d7d-11f0-365e-f1724fc39b3c
# ╠═5b3b5f83-058a-4582-9482-9b2fc5892692
# ╠═b9ef4450-d619-4e67-9ff3-9c8045b3863d
# ╠═f8573e77-38d8-467c-a123-9a83f73e8970
# ╠═19476521-78bd-44d6-a302-6e971f9fc095
# ╠═cba7fc75-8363-4ffa-b5a6-6e7d34363813
# ╠═583a9aee-08eb-4f5a-95ef-d0087eb98cbc
# ╠═39b011c6-f511-42dd-befc-eaf3fd17ea1a
# ╠═a94ab132-2949-4292-94d3-46db64809749
# ╠═7142e579-224c-474d-966f-461f8ce82e3a
# ╠═b5f399bf-5713-4f26-afb0-2d5771dbbc6f
# ╠═11e43f25-3fa8-4f2d-bba5-1773a4989178
# ╠═b3b16805-b5b8-4782-a49c-15029b3a749d
# ╠═8e8cc911-5a3c-4ec5-9498-cd60a09fd276
# ╠═b7508873-d85d-4502-8845-b779feb4666c
# ╠═f724b5fc-b506-4442-b6b0-92b8fc9ad16b
# ╠═55941a45-56f4-48c8-b823-880bdecaca29
# ╠═5d65761c-c47f-49d2-b4e8-3b758bc96e8b
# ╠═e76e1ca2-c736-42e7-8ab7-1f7bd09086e0
# ╠═fe676ab1-4ec5-4c3b-beb1-62c68f4f87ee
# ╠═53555a56-6602-4eb1-84a7-a64f4d8cb41e
# ╠═abb8b3bf-6d9f-4505-b13f-34de69460c51
# ╠═5f931b55-8578-4958-b6dc-24f4cfb3c011
# ╠═d155780e-74ed-48f4-98d9-45f23ff20e3f
# ╠═d4d1fcb8-ff53-4f49-bd5e-47e12672d8dc
# ╠═8e1005ba-ecff-4b97-be2c-3683529fa3c5
# ╠═096d3fd6-38e6-4da8-8a81-6339c3624ac5
# ╠═802b9639-4766-4209-899a-5a6fc1cf2d24
# ╠═6c9cb955-f2b3-41fb-85af-52873837e127
# ╠═49c7448d-824e-4604-9c90-c28e45f232d4
# ╠═7732c78d-25be-4648-a386-4d455b8bd0d5
# ╠═92d1a34c-9ca0-4f0f-a289-515b04f1e6d6
# ╠═c7cc7bfa-2ae3-44ef-8b98-5e5afbefcf2e
# ╠═fa4b76f1-f3e1-49ef-9054-404bfa942786
# ╠═2fe0ccbc-0b38-415c-995e-987017bcb499
# ╠═b93360f2-5890-4138-b7e8-53a9d6f539e4
# ╠═968648d8-54f2-4485-8209-8c22d4b63a8a
# ╠═d096680e-a09c-42cf-998e-b3304c509932
# ╠═11a04a68-91d0-4595-9296-144b3e21f478
# ╠═d2a8345c-9c44-43d7-a9ec-6ee4ab58674d
# ╠═b0a4c0be-12b9-416c-a45c-5fe35fbbd168
# ╠═89354633-e220-45f2-989e-c5575acd2988
# ╠═f4ff7938-f16e-4159-95e8-cb98c59a9d80
# ╠═fd39234f-6c62-42e3-bde9-9a9e565fa519
# ╠═f9a54f94-b109-4ae0-abe1-f44fdedc0d25
# ╠═389be056-e3d5-45fc-b80c-19f7cd4d5433
# ╠═e7cb1f63-130c-4e75-af5d-779fc1a4c755
# ╠═4d0c63a1-89e0-4d89-846c-e1501dbc2696
# ╠═60b16840-2530-41df-bbc1-12b1066e9829
# ╠═7927a18c-6942-42bd-ac6b-2ff720b155d0
# ╠═62da10f6-0d5c-41c0-a985-d15c946f5b84
# ╠═12960d51-135c-47cd-ab86-f2ab5bacef08
# ╠═d34b9105-09d4-4876-842d-bcf74249cca9
# ╠═80068e3e-e9c4-47c5-b6cf-e10ff2de19ea
# ╠═783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
# ╠═57692746-065d-4d75-8a41-7ffafd69550e
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
# ╠═d92c4660-4950-4adf-aceb-bc94663411c6
# ╠═6f3a242a-6dff-432a-b30e-1b7ee1d234fc
# ╠═3a8e7300-d209-4e50-ace8-c4b841d71f42
# ╠═af233a9a-2a3b-4b23-a982-c76d4a4c16f2
# ╠═7360e3df-c583-4894-86a6-5654b50a389c
# ╠═1c6507e3-821f-4e07-b41e-b3c106df3671
# ╠═e04c9c53-7ae2-474b-87bc-97bd58f458fa
# ╠═86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
# ╟─d1c82dbe-8493-464d-bdac-f505657678d6
# ╠═8aa50fbb-74a9-424d-9a84-e717247c65a9
# ╠═1b250af1-ffe5-488c-9eee-c1ca41085901
# ╠═77f5c1c3-11bb-4976-af44-1b488f08de6b
# ╠═3115a147-d41b-4ab6-9ad9-0f3b30fb4504
# ╠═f175508e-f0ea-48a9-ba72-d2e21215de1d
# ╠═5042becb-29c6-447e-84ad-a965a9961992
# ╠═1015d8f6-d32c-4920-b688-ab1848cd9c5a
# ╠═64ed321d-3050-465e-8c4f-fea826265779
# ╠═01391609-5034-4dbd-8319-10ac82126cfc
# ╠═264fcbb0-ec40-42b3-ab05-01b2171da6f2
# ╠═0a8baaba-06a6-4d44-bbc7-437e740be7b2
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═291c9e90-a915-4b35-a134-745ef253a72a
# ╠═862d9a0c-ecc6-4178-80f4-74fc71a79f38
# ╠═d9477e2f-2bad-4243-a11f-393f5b3405c7
# ╠═efd033eb-bde9-4f4e-85f0-777125555edd
# ╠═fa0da6ef-68ba-4db2-aec0-26fd3a280b5a
# ╠═57100463-320e-4768-9249-4d38b5b4f6b4
# ╠═d4a11884-b670-4ceb-9444-36d3bcc6e8a7
# ╠═81181098-822a-48f9-b776-680735de6430
# ╠═4b7e4552-ffb9-4f9a-a25f-f45178b2fbd3
# ╠═f0aab374-2e76-47e0-a62c-a5ec21129767
# ╠═a3c2e539-eda7-4784-b90c-43fcf6982347
# ╠═83df15d6-5493-48ef-a8f4-29b8fc375da3
