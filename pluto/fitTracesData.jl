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
	using JPELT
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
	using Changepoints
	using Distributions 
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

# ╔═╡ 23531712-6da1-458a-aaff-36acca0b79f9
names(JPELT)

# ╔═╡ b3b16805-b5b8-4782-a49c-15029b3a749d
names(BoldLab)

# ╔═╡ 8e8cc911-5a3c-4ec5-9498-cd60a09fd276
set_log_level(SimpleLogger.WARN)

# ╔═╡ b7508873-d85d-4502-8845-b779feb4666c
debug("test")

# ╔═╡ 9cad012b-4fb3-49e1-a00b-5814991c2b3a
begin
	jbl = ENV["JBoldLab"]
	naph3f, naph3c = fm_from_csv(joinpath(jbl, "data"), 
									  "naph3_emission.csv",
									  "naph3_abs.csv")
end

# ╔═╡ 66f96c28-0883-4ef6-a777-9a8980922d60
ActiveFilters = Dict("Filter3"=>[500.0, 520.0],
					 "Filter4"=>[524.0, 544.0],
			  "Filter5"=>[554.0, 568.0],
			  "Filter6"=>[576.0, 596.0],
			  "Filter7"=>[605.0, 625.0],
			  "Filter8"=>[633.0, 647.0])

# ╔═╡ 555f27e7-a820-4830-b47a-7bb50458cdc5
plot_molecule_emission(naph3f, naph3c, ActiveFilters)

# ╔═╡ 82e619dd-27d1-41b2-81e5-d102b972b3d8
compute_filter_coverage(naph3f, naph3c, ActiveFilters, λ_i=500.0, λ_f=650.0)

# ╔═╡ f724b5fc-b506-4442-b6b0-92b8fc9ad16b
	list_subfolders(dir) = [basename(e) for e in readdir(dir; join=true) if isdir(e) && !startswith(basename(e), ".")]


# ╔═╡ 95804d28-2720-4a5a-8003-05b1bc86f00c
begin
    BASEDIRS = Dict("CampaignApril2025" =>"/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign"	    )
    base_directory(label) = BASEDIRS[label]
	root_dir =  base_directory("CampaignApril2025")
end

# ╔═╡ 89ac3496-c36a-4472-a21c-bf4fca5e79d6
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

# ╔═╡ 6c7fb057-3c32-46c1-8a1d-6d2dae824a02
@bind folder_week Select(subdirs)

# ╔═╡ 626abdf3-401b-417d-97f2-23982570891f
begin
	if occursin("Week1", folder_week) ||  occursin("Week2", folder_week) ||  occursin("Week3", folder_week)
		path_week = joinpath(root_dir, folder_week, "Data")
	else
		path_week = joinpath(root_dir, folder_week)
	end
	subdirs2, _, _ = scan_level(path_week)
	@bind folder_day Select(subdirs2) 
end

# ╔═╡ a841c738-f6e3-4cc8-8ac0-f7d328d4677b
subdirs2

# ╔═╡ b7903ec3-c815-4574-9018-9c4202cf95ef
begin
	path_day = joinpath(path_week, folder_day)
	subdirs3, _, _ = scan_level(path_day)
	@bind folder_scan Select(subdirs3) 
end

# ╔═╡ bbf32af5-9175-4c9e-b71e-66e00a93b8f9
subdirs3

# ╔═╡ d5d957ad-89f3-4875-9ea2-91a0409d55ee
begin
	path_scan = joinpath(path_day, folder_scan)
	subdirs4, _, _ = scan_level(path_scan)
	if length(subdirs4) >0
		@bind folder_field Select(subdirs4) 
	else
		imf = path_scan
		folder_field = "" 
	end
	
end

# ╔═╡ d825fc83-e398-4319-a71f-ea67b97cf08f



# ╔═╡ e0040358-fab4-489d-8e4f-fa7128841ded
function find_filter(imf, Filters)
	
	ff = split(imf, "/")[end]
	
	i = findfirst("Filter", ff)
	suffix = i === nothing ? "" : ff[i.start:end]
	suffix, Filters[suffix]
end

# ╔═╡ f8ea2fe9-20af-4283-924f-d13729c4eca0
begin
	if length(folder_field) == 0 
		flt, fltv = find_filter(imf, ActiveFilters)
	elseif occursin("Filter", folder_field)
		flt, fltv = find_filter(folder_field, ActiveFilters)
	end
	md"""
		- Analyzing Fiter ? $(flt)
		- Filter band (in nm) = $(vect_to_fstr(fltv, "%.1f")) 
		"""
end

# ╔═╡ 2c50b5b4-a227-423a-bd4b-812fb4cd425c
begin
	if length(folder_field) > 0
		path_tiff = joinpath(path_scan, folder_field)
	else 
		path_tiff =imf
	end
	imxt = load_tif_stack_int16(path_tiff, pedestal=0.0)
	plot_frames(imxt; nscale=20)
end

# ╔═╡ 797aaf3e-f284-4f62-bb1a-5514190ca422
function png_from_path(path)
	
	parts = split(path, '/')
	folder1 = parts[end - 1]  # W3_SVE0_1
	folder2 = parts[end]      # Field1
	string("stepAnalysis","/",folder1,"_", folder2, ".png")
	#join(["stepAnalysis", folder1 * "_" * folder2], '/')
end

# ╔═╡ 4e4252bc-3b6e-4507-8817-dfe9558b51b6
md"""
- Show global statistics
"""

# ╔═╡ 783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
begin
	totalI, meanI, stdI = get_stats(imxt; bck=0.0)
	plot_stats(totalI, meanI, stdI)
end

# ╔═╡ 11a04a68-91d0-4595-9296-144b3e21f478
begin
	σ = 10.0 #average over large number of pixels 
	nlf = 5  # number of frames to average
	background = compute_background_from_stack(imxt; σ= σ, nlf=nlf)
	hbk, pbk = step_hist(vec(background);
              nbins=40,
              xlim=(1500.0, 3500.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="residual intensity background")
	#histogram(vec(background))
	plot(pbk, size=(600, 300))
end

# ╔═╡ eb6b49be-ec65-447d-afe6-2c7580bd6693
md"""
- Measure the average background.
- Average over the last $(nlf) frames
"""

# ╔═╡ d2a8345c-9c44-43d7-a9ec-6ee4ab58674d
md"""
- substrack background from stack
"""

# ╔═╡ b0a4c0be-12b9-416c-a45c-5fe35fbbd168
begin
	imgbsub = subtract_background_from_stack(imxt, background)
	vim = vec(imgbsub)
	icut = 5.0

	hibs, pibs = step_hist(vim;
              nbins=50,
              xlim=(0.0, 100.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="IBS")
	hibs2, pibs2 = step_hist(vim[vim.>icut];
              nbins=50,
              xlim=(0.0, 350.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="IBS with cut = $(icut) ")
	#histogram(vec(background))
	plot(pibs, pibs2, size=(800, 400))
	#histogram(vim[vim.>50])
end

# ╔═╡ 7c267c22-6d75-46a4-a459-6a03d74da0cb
md"""
- plot some frames background subtracted
"""

# ╔═╡ 89354633-e220-45f2-989e-c5575acd2988
plot_frames(imgbsub; nscale=20)

# ╔═╡ 83b89e38-10e8-4ce1-9639-29b3afb38a64
md"""
- Denoise the stack and plot some frames
"""

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

	hibsf2, pibsf2 = step_hist(vi2m[vi2m.>icut];
              nbins=50,
              xlim=(0.0, 1000.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="IBSD with cut = $(icut) ")
	#histogram(vec(background))
	plot(pibsf2, size=(600, 300))
	#histogram(vi2m[vi2m.>30])
end

# ╔═╡ 62da10f6-0d5c-41c0-a985-d15c946f5b84
begin
    @bind show_peaks CheckBox(true)
end

# ╔═╡ 12960d51-135c-47cd-ab86-f2ab5bacef08
@bind nframe NumberField(0:199, default=1)

# ╔═╡ d34b9105-09d4-4876-842d-bcf74249cca9
@bind pthr NumberField(0:0.1:200.0, default=50.0)

# ╔═╡ 74deedfc-a344-4373-9fc9-db22e83d48ac
md"""
#### Get traces from the imst matrix
"""

# ╔═╡ 3a8e7300-d209-4e50-ace8-c4b841d71f42
md"""
#### Plot traces
"""

# ╔═╡ 7360e3df-c583-4894-86a6-5654b50a389c
#plot_traces(TRZS, peaks; ftrz=26, ltrz=50,  figsize=(1500,1500))

# ╔═╡ 22fc945b-532f-40d9-8e72-ba3a5a58f95a
md"""
- Select traces for fit illustrations 
"""

# ╔═╡ d92c4660-4950-4adf-aceb-bc94663411c6
function select_trace(TRZS, df::DataFrame, row::Int) 
	i = df.i[row]
	j = df.j[row]
    trace = TRZS[i, j]  # safe, only defined where needed
	i,j, trace
end

# ╔═╡ 51062b57-15c3-4288-bd14-74afb1bd08d2
md"""
## Autostepfinder
- We call autostepfinder with thr=0, and a guess on the number of interactions. The resulting S-curve can be used to refine the actual interactions. 
"""

# ╔═╡ 567d39c7-f9ae-453a-8e9a-25672b399264
md"- select the traze to fit"

# ╔═╡ 1c6507e3-821f-4e07-b41e-b3c106df3671
@bind ntrz NumberField(0:100, default=1)

# ╔═╡ 415f6577-d2e6-4f7d-b1ee-5e1695ce98e3
md"""
- In standard fit, keep thr=0, used as parameter only to show that increasing thr values ruins the fit often
"""

# ╔═╡ e04c9c53-7ae2-474b-87bc-97bd58f458fa
@bind thr NumberField(0.0:0.1:1.1, default=0.0)

# ╔═╡ 3ad05178-0420-41af-8b6d-19994428bc5e
md"- Guess on the number of interactions. Check the value of the S-Curve. The optimal number of interactions should be set to the peak of the S-Curve +1. If the data is not noisy, the best fit will find the minimum number of steps by itself, but noisy data may result in higher steps"

# ╔═╡ 86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
@bind niter NumberField(0:10, default=10)

# ╔═╡ 7929e575-bddf-4da5-ba21-6fdcc71f81da
begin
	
	md"""
## PELT
"""
end

# ╔═╡ 77f5c1c3-11bb-4976-af44-1b488f08de6b
md"""
## Fits to data
"""

# ╔═╡ 3115a147-d41b-4ab6-9ad9-0f3b30fb4504
begin
    @bind fitASF10 CheckBox(false)
end

# ╔═╡ 816cc8ba-db93-49ff-ba0a-0e42138c1ca7
begin
    @bind fitASF3 CheckBox(false)
end

# ╔═╡ 5349482c-ba2a-4a85-ae2c-d6b20c8ae6e1
begin
    @bind fitPELT CheckBox(false)
end

# ╔═╡ f175508e-f0ea-48a9-ba72-d2e21215de1d
md"""
- Set the vaue of iterations (number of steps to be sought in the data)
"""

# ╔═╡ 25b059da-37e0-4498-b7cf-4d87ac26e126
if fitASF3
	niter2 = 3
md"""
- There is no clear recipe to manipulate the algorithm parameters. One can, for example, reduce the number of iterations, knowing that the data has few steps. As an example let us take niter = $(niter2): 
"""
end

# ╔═╡ 60628fe1-6060-4d1f-b7c8-e1d9bd4032a5
if fitPELT
	md"""
	## Fits with PELT
	"""
end

# ╔═╡ b9e5ce1f-4a86-4377-9d50-6732fae9516c
readdir("stepAnalysis")

# ╔═╡ 76f1b23b-b6c4-43a0-9e53-3b2ff553f167
md"""
## PELT code
"""

# ╔═╡ 1be93848-5757-48e6-8d5b-638cb11c4a61
md"""
## Functions
"""

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

# ╔═╡ 4d0c63a1-89e0-4d89-846c-e1501dbc2696
begin
	#threshold = 10.0  
	#min_area = 10
	peaks = detect_local_maxima(imgd[:, :, nframe]; threshold=pthr, dx=0, dy=0)
	md"""
	- Search for local maxima
	- found $(size(peaks)[1]) molecule candidates 
	- in frame $(nframe) 
	- with thr = $(pthr)
	"""
end

# ╔═╡ 60b16840-2530-41df-bbc1-12b1066e9829
peaks

# ╔═╡ 7927a18c-6942-42bd-ac6b-2ff720b155d0
begin
hisp, pisp = step_hist(peaks.intensity;
              nbins=20,
              xlim=(0.0, 1000.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="Intensity of selected peaks with cut = $(pthr) ")
	
	plot(pisp, size=(600, 300))
	
#histogram(peaks.intensity, bins=20)
end

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

# ╔═╡ 57692746-065d-4d75-8a41-7ffafd69550e
md"""
- total intensity = $(Float64(sum(totalI)))
- number of detected peaks = $(size(peaks)[1])
"""

# ╔═╡ 1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
TRZS = build_sparse_traces(imgd, peaks)

# ╔═╡ 02ba78f2-1de5-46d1-a8d6-d6895b2d8771
peaks

# ╔═╡ 6f3a242a-6dff-432a-b30e-1b7ee1d234fc
i,j, tz = select_trace(TRZS, peaks, ntrz)

# ╔═╡ d1c82dbe-8493-464d-bdac-f505657678d6
begin
	dataX, FitX, S_curve, best_shot, iter, cc = fit_traces(tz, niter=niter, tresH=thr, sel="core")
	sth, stt, stl = getsteps(FitX)
	

	md"""
	- Fit results
	- best shot =$(best_shot)
	- number of steps = $(length(stl))
	- Step heights: = $(vect_to_fstr(sth, "%.d"))
	- Step times =$(vect_to_fstr(stt, "%.d"))
	"""
	
end

# ╔═╡ 8aa50fbb-74a9-424d-9a84-e717247c65a9
plotfit(dataX, FitX, S_curve, i,j )

# ╔═╡ ca42733e-3343-4676-a0a8-e347bdd600cf
begin
	cro = pelt_fit(tz)
	pfit, pen, cost = steps_from_pelt(tz, cro; ic=1)
	pfit2, pen2, cost2 = steps_from_pelt(tz, cro; ic=2)
	if pen > 0 && pen2 > 0
		plotfit2(tz, pfit, pfit2; figsize=(1000, 800))
	elseif pen > 0 
		plotfit2(tz, pfit, pfit; figsize=(1000, 800))
	else
		md"""
		- PELT did not find a good fit
		"""
	end
end

# ╔═╡ 5042becb-29c6-447e-84ad-a965a9961992
if fitASF10
	dfs, zI, zJ, zDX, zFX, zSC  =find_fit_candidates2(TRZS, peaks;  sel="core", ped=0.0, niter=niter, thr=thr)
	dfs
end

# ╔═╡ dc5bf3d6-9bac-42a9-b914-ab45e744a11d
if fitASF3
	dfs2, zI2, zJ2, zDX2, zFX2, zSC2  =find_fit_candidates2(TRZS, peaks;  sel="core", ped=0.0, niter=niter2, thr=thr)
	dfs2
end

# ╔═╡ b7b1a6bc-aa85-4daa-a77b-e69267d2c633
if fitPELT
	pdf, pI, pJ, pDX, pFX, pSC  =find_fit_pelt(TRZS, peaks)
	pdf
end

# ╔═╡ 7cfbc9b3-965b-45ca-bc33-f3a5f57b2ea3
if fitPELT
	hstph, pstph = step_hist(pdf.stepHeight;
	              nbins=25,
	              xlim=(0.0, 200.0),
	              logy=false,
	              xlabel=" # step height",
	              ylabel=" #entries ")
		plot(pstph, size=(600,300))
end


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
plot_traces(TRZS, peaks; ftrz=1, ltrz=6,  figsize=(1500,1500))

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
if fitASF10
	plot_sc(zSC, zI,zJ; plotsel="4x4", figsize=(1500,1500))
end

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
if fitASF10
	plot_fits(zDX, zFX, zI, zJ; plotsel="5x5")
end

# ╔═╡ c0919a79-785f-4cae-b576-f1cece2f680f
if fitASF3
	plot_fits(zDX2, zFX2, zI2, zJ2; plotsel="5x5")
end

# ╔═╡ 3410b51e-ed17-4970-8290-cd98f11df26e
if fitPELT
	plot_fits(pDX, pFX, pI, pJ; plotsel="5x5")
end

# ╔═╡ d4a11884-b670-4ceb-9444-36d3bcc6e8a7
function plotsdf_pelt(sdf; ppath="")

		
	# take only one row per (i,j), keeping the first nstep
	df_pix = unique(sdf, [:i, :j, :nstep])

	
	hnstep, pnstep = step_hist(df_pix.nstep;
              nbins=10,
              xlim=(0.0, 10.0),
              logy=false,
              xlabel=" # steps",
              ylabel=" #entries ")

	
	hstptime, pstptime = step_hist(sdf.stepTime;
              nbins=25,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" step time",
              ylabel=" #entries ")

	hstph, pstph = step_hist(sdf.stepHeight;
              nbins=25,
              xlim=(0.0, 250.0),
              logy=false,
              xlabel=" # step height",
              ylabel=" #entries ")


	hstplx, pstplx = step_hist(sdf.stepLength;
              nbins=50,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" # step length",
              ylabel=" #entries ")
	
	plt = plot(pnstep, pstptime, pstph, pstplx, 
		 layout=(2, 2), size=(800,600))
	
	savefig(plt, ppath)
	plt
	
	
end

# ╔═╡ 7ddb353f-ef75-4ffd-9273-b02a3321474c
if fitPELT
	ppng = png_from_path(path_tiff)
	plotsdf_pelt(pdf; ppath=ppng)
end

# ╔═╡ d9d972f7-e196-4232-a062-c963c9147801
md"""
- Saving plots to $(ppng)
"""

# ╔═╡ f9c4fe44-b596-4d54-a07a-47a9ca76c108
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
              nbins=25,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" step time",
              ylabel=" #entries ")

	hstph, pstph = step_hist(sdf.stepHeight;
              nbins=25,
              xlim=(0.0, 250.0),
              logy=false,
              xlabel=" # step height",
              ylabel=" #entries ")

	hstphx, pstphx = step_hist(df_pix.stepHeightMin;
              nbins=25,
              xlim=(0.0, 250.0),
              logy=false,
              xlabel=" # step height min",
              ylabel=" #entries ")

	hstplx, pstplx = step_hist(sdf.stepLength;
              nbins=25,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" # step length",
              ylabel=" #entries ")
	
	plot(pstphx, size=(400,200))
	
	plot(pbshot, pnstep, pstptime, pstph, pstphx, pstplx, 
		 layout=(3, 2), size=(800,600))
	
	
end

# ╔═╡ 0a8baaba-06a6-4d44-bbc7-437e740be7b2
if fitASF10
	plotsdf(dfs)
end

# ╔═╡ 931b5be3-d6cf-49cc-82e8-050776406a5c
if fitASF3
	plotsdf(dfs2)
end

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

# ╔═╡ 264fcbb0-ec40-42b3-ab05-01b2171da6f2
if fitASF10
	md"""
	### The ASF algorithm 
	
	- number of fitted molecules = $(count_mol(dfs))
	- threshold = $(thr)
	- ninter = $(niter)

	The algorithm succeeds almost always to fit the data, but finds too many steps and the results are hard to intepret, in particular the flat tale of sept times above 100 time units. The step height at 200 corresponds with the baseline, also showing that many steps are fake. 
	"""
end

# ╔═╡ 43c3f8da-d787-4e46-89b4-433ddca30843
if fitASF3
	md"""
	
	- number of fitted molecules = $(count_mol(dfs2))
	- threshold = $(thr)
	- ninter = $(niter2)

	The distribution in step time is better (not the distribution in step length), but the efficiency of the algorithm is smaller (as expected). 
	"""
end

# ╔═╡ 300fcee3-17cd-44f4-bd76-6c8bc03b30a4
if fitPELT
md"""
### PELT
- number of fitted molecules = $(count_mol(pdf))
- The algorithm can be run without parameters. The efficiency is higher than ASF. The algorithm finds always 2 steps as expected in the data and the expected distributions in step hegith and time. 
"""
end

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
# ╠═23531712-6da1-458a-aaff-36acca0b79f9
# ╠═b3b16805-b5b8-4782-a49c-15029b3a749d
# ╠═8e8cc911-5a3c-4ec5-9498-cd60a09fd276
# ╠═b7508873-d85d-4502-8845-b779feb4666c
# ╠═9cad012b-4fb3-49e1-a00b-5814991c2b3a
# ╠═66f96c28-0883-4ef6-a777-9a8980922d60
# ╠═555f27e7-a820-4830-b47a-7bb50458cdc5
# ╠═82e619dd-27d1-41b2-81e5-d102b972b3d8
# ╠═f724b5fc-b506-4442-b6b0-92b8fc9ad16b
# ╠═95804d28-2720-4a5a-8003-05b1bc86f00c
# ╠═89ac3496-c36a-4472-a21c-bf4fca5e79d6
# ╠═6c7fb057-3c32-46c1-8a1d-6d2dae824a02
# ╠═626abdf3-401b-417d-97f2-23982570891f
# ╠═a841c738-f6e3-4cc8-8ac0-f7d328d4677b
# ╠═b7903ec3-c815-4574-9018-9c4202cf95ef
# ╠═bbf32af5-9175-4c9e-b71e-66e00a93b8f9
# ╠═d5d957ad-89f3-4875-9ea2-91a0409d55ee
# ╠═d825fc83-e398-4319-a71f-ea67b97cf08f
# ╠═e0040358-fab4-489d-8e4f-fa7128841ded
# ╠═f8ea2fe9-20af-4283-924f-d13729c4eca0
# ╠═2c50b5b4-a227-423a-bd4b-812fb4cd425c
# ╠═797aaf3e-f284-4f62-bb1a-5514190ca422
# ╠═4e4252bc-3b6e-4507-8817-dfe9558b51b6
# ╠═783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
# ╠═eb6b49be-ec65-447d-afe6-2c7580bd6693
# ╠═11a04a68-91d0-4595-9296-144b3e21f478
# ╠═d2a8345c-9c44-43d7-a9ec-6ee4ab58674d
# ╠═b0a4c0be-12b9-416c-a45c-5fe35fbbd168
# ╠═7c267c22-6d75-46a4-a459-6a03d74da0cb
# ╠═89354633-e220-45f2-989e-c5575acd2988
# ╠═83b89e38-10e8-4ce1-9639-29b3afb38a64
# ╠═f4ff7938-f16e-4159-95e8-cb98c59a9d80
# ╠═fd39234f-6c62-42e3-bde9-9a9e565fa519
# ╠═4d0c63a1-89e0-4d89-846c-e1501dbc2696
# ╠═60b16840-2530-41df-bbc1-12b1066e9829
# ╠═7927a18c-6942-42bd-ac6b-2ff720b155d0
# ╠═62da10f6-0d5c-41c0-a985-d15c946f5b84
# ╠═12960d51-135c-47cd-ab86-f2ab5bacef08
# ╠═d34b9105-09d4-4876-842d-bcf74249cca9
# ╠═80068e3e-e9c4-47c5-b6cf-e10ff2de19ea
# ╠═57692746-065d-4d75-8a41-7ffafd69550e
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
# ╠═02ba78f2-1de5-46d1-a8d6-d6895b2d8771
# ╠═3a8e7300-d209-4e50-ace8-c4b841d71f42
# ╠═af233a9a-2a3b-4b23-a982-c76d4a4c16f2
# ╠═7360e3df-c583-4894-86a6-5654b50a389c
# ╠═22fc945b-532f-40d9-8e72-ba3a5a58f95a
# ╠═d92c4660-4950-4adf-aceb-bc94663411c6
# ╠═6f3a242a-6dff-432a-b30e-1b7ee1d234fc
# ╠═51062b57-15c3-4288-bd14-74afb1bd08d2
# ╟─567d39c7-f9ae-453a-8e9a-25672b399264
# ╠═1c6507e3-821f-4e07-b41e-b3c106df3671
# ╟─415f6577-d2e6-4f7d-b1ee-5e1695ce98e3
# ╠═e04c9c53-7ae2-474b-87bc-97bd58f458fa
# ╠═3ad05178-0420-41af-8b6d-19994428bc5e
# ╠═86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
# ╠═d1c82dbe-8493-464d-bdac-f505657678d6
# ╠═8aa50fbb-74a9-424d-9a84-e717247c65a9
# ╠═7929e575-bddf-4da5-ba21-6fdcc71f81da
# ╠═ca42733e-3343-4676-a0a8-e347bdd600cf
# ╠═77f5c1c3-11bb-4976-af44-1b488f08de6b
# ╠═3115a147-d41b-4ab6-9ad9-0f3b30fb4504
# ╠═816cc8ba-db93-49ff-ba0a-0e42138c1ca7
# ╠═5349482c-ba2a-4a85-ae2c-d6b20c8ae6e1
# ╠═f175508e-f0ea-48a9-ba72-d2e21215de1d
# ╠═5042becb-29c6-447e-84ad-a965a9961992
# ╠═1015d8f6-d32c-4920-b688-ab1848cd9c5a
# ╠═64ed321d-3050-465e-8c4f-fea826265779
# ╠═264fcbb0-ec40-42b3-ab05-01b2171da6f2
# ╠═0a8baaba-06a6-4d44-bbc7-437e740be7b2
# ╠═25b059da-37e0-4498-b7cf-4d87ac26e126
# ╠═dc5bf3d6-9bac-42a9-b914-ab45e744a11d
# ╠═c0919a79-785f-4cae-b576-f1cece2f680f
# ╠═931b5be3-d6cf-49cc-82e8-050776406a5c
# ╠═43c3f8da-d787-4e46-89b4-433ddca30843
# ╠═60628fe1-6060-4d1f-b7c8-e1d9bd4032a5
# ╠═b7b1a6bc-aa85-4daa-a77b-e69267d2c633
# ╠═3410b51e-ed17-4970-8290-cd98f11df26e
# ╠═7ddb353f-ef75-4ffd-9273-b02a3321474c
# ╠═d9d972f7-e196-4232-a062-c963c9147801
# ╠═b9e5ce1f-4a86-4377-9d50-6732fae9516c
# ╠═7cfbc9b3-965b-45ca-bc33-f3a5f57b2ea3
# ╠═300fcee3-17cd-44f4-bd76-6c8bc03b30a4
# ╠═76f1b23b-b6c4-43a0-9e53-3b2ff553f167
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═e7cb1f63-130c-4e75-af5d-779fc1a4c755
# ╠═291c9e90-a915-4b35-a134-745ef253a72a
# ╠═862d9a0c-ecc6-4178-80f4-74fc71a79f38
# ╠═d9477e2f-2bad-4243-a11f-393f5b3405c7
# ╠═efd033eb-bde9-4f4e-85f0-777125555edd
# ╠═fa0da6ef-68ba-4db2-aec0-26fd3a280b5a
# ╠═d4a11884-b670-4ceb-9444-36d3bcc6e8a7
# ╠═f9c4fe44-b596-4d54-a07a-47a9ca76c108
# ╠═81181098-822a-48f9-b776-680735de6430
# ╠═4b7e4552-ffb9-4f9a-a25f-f45178b2fbd3
# ╠═f0aab374-2e76-47e0-a62c-a5ec21129767
# ╠═a3c2e539-eda7-4784-b90c-43fcf6982347
# ╠═83df15d6-5493-48ef-a8f4-29b8fc375da3
