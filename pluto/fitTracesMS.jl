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
	#using JStepFinder
	#using StepAnalysis
	#using LabStepAnalysis
	#using StepPreprocessing
	#using JPELT
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

# ╔═╡ c5d3b3b0-d5f2-4102-8bc7-e6af113e8e63
md"""
- Trace Analysis (MS stands for Mikel Setup)
"""

# ╔═╡ cba7fc75-8363-4ffa-b5a6-6e7d34363813
import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W

# ╔═╡ 62925080-87f3-4c84-8fc4-1b9ae36ec7c5
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
                :(eval(x) = $(Expr(:core, :eval))($name, x)),
                :(include(x) = $(Expr(:top, :include))($name, x)),
                :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
                :(include($path))))
    m
end

# ╔═╡ 6586349e-37f9-4cd2-8513-d05558ec3d6b
begin
	jn = ingredients(string(ENV["JBoldLab"],"/src/BoldLab.jl"))
end

# ╔═╡ 583a9aee-08eb-4f5a-95ef-d0087eb98cbc
names(SimpleLogger)

# ╔═╡ 39b011c6-f511-42dd-befc-eaf3fd17ea1a
#names(StepAnalysis)

# ╔═╡ a94ab132-2949-4292-94d3-46db64809749
#names(LabStepAnalysis)

# ╔═╡ 7142e579-224c-474d-966f-461f8ce82e3a
#names(StepPreprocessing)

# ╔═╡ b5f399bf-5713-4f26-afb0-2d5771dbbc6f
#names(JStepFinder)

# ╔═╡ 11e43f25-3fa8-4f2d-bba5-1773a4989178
names(histos)

# ╔═╡ 23531712-6da1-458a-aaff-36acca0b79f9
#names(JPELT)

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
			  "Filter8"=>[633.0, 647.0],
			   "Filter12"=>[500.0, 800.0])

# ╔═╡ 967ef731-1cee-4091-8f98-4a7f080058be
DoubleFilters = Dict("Filter1"=>[473.0, 552.0],
					 "Filter2"=>[600.0, 750.0])

# ╔═╡ 555f27e7-a820-4830-b47a-7bb50458cdc5
plot_molecule_emission(naph3f, naph3c, DoubleFilters)

# ╔═╡ 82e619dd-27d1-41b2-81e5-d102b972b3d8
compute_filter_coverage(naph3f, naph3c, ActiveFilters, λ_i=500.0, λ_f=650.0)

# ╔═╡ f724b5fc-b506-4442-b6b0-92b8fc9ad16b
	list_subfolders(dir) = [basename(e) for e in readdir(dir; join=true) if isdir(e) && !startswith(basename(e), ".")]


# ╔═╡ 95804d28-2720-4a5a-8003-05b1bc86f00c
begin
    BASEDIRS = Dict("CampaignApril2025" =>"/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign",
	"MikelSeptember2025" =>"/Users/jjgomezcadenas/BoldLab/BOLD/Data1/ZFFL179_NAPH3_SIL_5e-7")
    base_directory(label) = BASEDIRS[label]
	root_dir =  base_directory("MikelSeptember2025")
end

# ╔═╡ 89ac3496-c36a-4472-a21c-bf4fca5e79d6
begin
   phys_dirs, _, _ = scan_level(root_dir)

	md"""
	- Select Physics dir (FREE or EVAP)
	"""
end


# ╔═╡ e3ce6b68-4186-41a5-a3f1-a458b78d8821
begin
	@bind phys_folder Select(phys_dirs, default="FREE") 
	path_samples = joinpath(root_dir, phys_folder)
   sample_dirs, _, _ = scan_level(path_samples)

	md"""
	- Select sample dir (A1, A2...)
	"""
end

# ╔═╡ 62c9f533-a314-4c80-b700-8b83bc3c8f83
@bind sample_folder Select(sample_dirs) 

# ╔═╡ c8770df9-727b-4095-b45c-d9d10877227e
begin
	path_files = joinpath(root_dir, phys_folder, sample_folder)
    _, _, all_files = scan_level(path_files)
	file_numbers = sort(unique([parse(Int, match(r"Image(\d+)_", f).captures[1]) for f in all_files]))
	md"""
	- Select field(an integer number)
	"""
end

# ╔═╡ 04601060-8310-4e76-a2d0-8e258d2d061f
@bind field Select(file_numbers) 

# ╔═╡ 9c6db6db-b5b7-49e6-b275-970a29ba69cb
begin
	filtered_files = filter(f -> occursin("Image$(field)_", f), all_files)
	tif_files = [joinpath(path_files, ff) for ff in filtered_files]
	imxt = load_tif_stack_int16(tif_files, pedestal=0.0)
	igreen = imxt[1:399, :, :]
	ired = imxt[400:800, :, :]
end

# ╔═╡ 91b5254a-b21a-4243-a838-811ddc90d515
plot_frames(igreen; colormap=:viridis, nscale=30)

# ╔═╡ 5a63ed41-9233-4654-a37c-d00013c35197
plot_frames(ired; colormap=:hot, nscale=30)

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
- Show global statistics for green channel 
"""

# ╔═╡ 783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
begin
	totalIg, meanIg, stdIg = get_stats(igreen; bck=0.0)
	plot_stats(totalIg, meanIg, stdIg)
end

# ╔═╡ 7e40b817-a1d3-45c6-b651-136031aad7b3
md"""
- Show global statistics for red channel 
"""

# ╔═╡ ee8175fb-f6dd-46fc-baed-ae82209af0ed
begin
	totalIr, meanIr, stdIr = get_stats(ired; bck=0.0)
	plot_stats(totalIr, meanIr, stdIr)
end

# ╔═╡ 11a04a68-91d0-4595-9296-144b3e21f478
begin
	σ = 10.0 #average over large number of pixels 
	nlf = 5  # number of frames to average
	gbackground = compute_background_from_stack(igreen; σ= σ, nlf=nlf)
	rbackground = compute_background_from_stack(ired; σ= σ, nlf=nlf)
	ghbk, gpbk = step_hist(vec(gbackground);
              nbins=40,
              xlim=(400.0, 600.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="residual intensity background GREEN")
	rhbk, rpbk = step_hist(vec(rbackground);
              nbins=40,
              xlim=(400.0, 600.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="residual intensity background RED")
	#histogram(vec(background))
	plot(gpbk, rpbk, layout=(1,2), size=(800, 400))
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
	gimg = subtract_background_from_stack(igreen, gbackground)
	vim = vec(gimg)
	icut = 5.0

	hibs, pibs = step_hist(vim;
              nbins=50,
              xlim=(0.0, 100.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="GREEN: IBS")
	hibs2, pibs2 = step_hist(vim[vim.>icut];
              nbins=50,
              xlim=(0.0, 350.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="GREEN: IBS with cut = $(icut) ")
	#histogram(vec(background))
	plot(pibs, pibs2, size=(800, 400))
	#histogram(vim[vim.>50])
end

# ╔═╡ 69734f52-a953-44e3-81e6-a8f911c4f4d0
begin
	rimg = subtract_background_from_stack(ired, rbackground)
	vim2 = vec(rimg)
	icut2 = 5.0

	rhibs, rpibs = step_hist(vim2;
              nbins=50,
              xlim=(0.0, 100.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="RED: IBS")
	rhibs2, rpibs2 = step_hist(vim2[vim2.>icut2];
              nbins=50,
              xlim=(0.0, 350.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="RED: IBS with cut = $(icut) ")
	#histogram(vec(background))
	plot(rpibs, rpibs2, size=(800, 400))
	#histogram(vim[vim.>50])
end

# ╔═╡ 7c267c22-6d75-46a4-a459-6a03d74da0cb
md"""
- plot some frames background subtracted
"""

# ╔═╡ 89354633-e220-45f2-989e-c5575acd2988
plot_frames(gimg; colormap=:viridis, nscale=30)

# ╔═╡ 5fcf28e8-0467-4f6f-b2e7-676693e8355e
plot_frames(rimg; colormap=:hot, nscale=30)

# ╔═╡ 83b89e38-10e8-4ce1-9639-29b3afb38a64
md"""
- Denoise the stack and plot some frames
"""

# ╔═╡ f4ff7938-f16e-4159-95e8-cb98c59a9d80
begin
	σ_denoise = 1.0
	gimgd = denoise_stack(gimg; σ=σ_denoise)
	plot_frames(gimgd; colormap=:viridis, nscale=30)
end

# ╔═╡ 156d523f-1b80-4cbd-912c-ed21a1327aff
begin
	rimgd = denoise_stack(rimg; σ=σ_denoise)
	plot_frames(rimgd; colormap=:hot, nscale=30)
end

# ╔═╡ fd39234f-6c62-42e3-bde9-9a9e565fa519
begin
	vi2m = vec(gimgd[:,:,2])

	hibsf2, pibsf2 = step_hist(vi2m[vi2m.>icut];
              nbins=50,
              xlim=(0.0, 2500.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="GREEN: IBSD with cut = $(icut) ")
	
	vi2r = vec(rimgd[:,:,2])

	rhibsf2, rpibsf2 = step_hist(vi2r[vi2r.>icut2];
              nbins=50,
              xlim=(0.0, 350.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="RED: IBSD with cut = $(icut2) ")
	plot(pibsf2, rpibsf2, size=(600, 300))
	#histogram(vi2m[vi2m.>30])
end

# ╔═╡ 62da10f6-0d5c-41c0-a985-d15c946f5b84
begin
    @bind show_peaks CheckBox(true)
end

# ╔═╡ 12960d51-135c-47cd-ab86-f2ab5bacef08
@bind nframe NumberField(0:199, default=1)

# ╔═╡ d34b9105-09d4-4876-842d-bcf74249cca9
@bind pthr NumberField(0:1.0:200.0, default=50.0)

# ╔═╡ 4d0c63a1-89e0-4d89-846c-e1501dbc2696
begin
	gpeaks = detect_local_maxima(gimgd[:, :, nframe]; threshold=pthr, dx=0, dy=0)
	md"""
	- GREEN: Search for local maxima in frame $(nframe)
	- found $(size(gpeaks)[1]) molecule candidates 
	- with thr = $(pthr)
	"""
end

# ╔═╡ a74eff14-b6b4-4e80-bfc7-5a1db6a32aba
begin
	rpeaks = detect_local_maxima(rimgd[:, :, nframe]; threshold=pthr, dx=0, dy=0)
	md"""
	- RED: Search for local maxima in frame $(nframe)
	- found $(size(rpeaks)[1]) molecule candidates 
	- with thr = $(pthr)
	"""
end

# ╔═╡ 7927a18c-6942-42bd-ac6b-2ff720b155d0
begin
	ghisp, gpisp = step_hist(gpeaks.intensity;
	              nbins=20,
	              xlim=(0.0, 2500.0),
	              logy=false,
	              xlabel=" intensity",
	              ylabel=" #entries ",
				  title="GREEN:  Int of sel peaks; cut = $(pthr) ")

	rhisp, rpisp = step_hist(rpeaks.intensity;
              nbins=20,
              xlim=(0.0, 2500.0),
              logy=false,
              xlabel=" intensity",
              ylabel=" #entries ",
			  title="RED: Int of sel peaks; cut = $(pthr) ")
	
	plot(gpisp, rpisp, size=(800, 400))
	
#histogram(peaks.intensity, bins=20)
end

# ╔═╡ 1333c10b-c751-414b-9d71-df9f6499172b
md"""
- Maximum intensity of green channel: $(maximum(gpeaks.intensity))
- Total intensity of green channel = $(Float64(sum(totalIg)))
- number of detected peaks in green channel = $(size(gpeaks)[1])
-
- Maximum intensity of red channel: $(maximum(rpeaks.intensity))
- Total intensity of red channel = $(Float64(sum(totalIr)))
- number of detected peaks in red channel = $(size(rpeaks)[1])
"""

# ╔═╡ 74deedfc-a344-4373-9fc9-db22e83d48ac
md"""
#### Get traces from the imst matrix
"""

# ╔═╡ 1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
gTRZS = build_sparse_traces(gimgd, gpeaks)

# ╔═╡ 574e8729-7610-47c0-be0a-767f4af15c7c
rTRZS = build_sparse_traces(rimgd, rpeaks)

# ╔═╡ 77f5c1c3-11bb-4976-af44-1b488f08de6b
md"""
# Fits to data
"""

# ╔═╡ b24005f5-e019-4f41-9193-ebf7e5d69dac
@bind thr NumberField(0.0:0.1:1.1, default=0.5)

# ╔═╡ 381bb69d-1b27-4403-9f95-bfde03fc1c15
@bind niter NumberField(0:100, default=10)

# ╔═╡ 56e7b887-3304-47e4-8c17-8268db12a46b
md"""
## Autostepfinder
- Autostepfinder depend on two parameters. 
	- threshold (thr)
	- number of interactions (niter)
- Default values: thr=$(thr), niter = $(niter) 
"""


# ╔═╡ 8d46c82c-6585-4d79-b5a9-2ee51d321e36
md"""
### Select traces for fit illustrations 
"""

# ╔═╡ a8b6fe7c-f328-4602-b7e2-240f53b78b64
@bind ntrz NumberField(0:1000, default=10)

# ╔═╡ c11c894f-cefa-4d50-9c27-19cb744ef7ef
begin
	#dataX, FitX, S_curve, best_shot, iter, cc = fit_traces(tz, niter=niter, tresH=thr, sel="core")
	#sth, stt, stl = getsteps(FitX)
	

	#md"""
	#- Fit results
	#- best shot =$(best_shot)
	#- number of steps = $(length(stl))
	#- Step heights: = $(vect_to_fstr(sth, "%.d"))
	#- Step times =$(vect_to_fstr(stt, "%.d"))
	#"""
	
end

# ╔═╡ 709275a8-af9f-4cea-9e9f-1c22c95cf48c
#md"""

#- Steps = $(length(afsFit.steps))
#- Heights = $(vect_to_fstr(afsFit.stepHeights, "%.d"))
#- Times = $(vect_to_fstr(afsFit.stepTimes, "%.d")) 
#"""

# ╔═╡ 073db006-22d0-4d27-8148-ef142b163ba2
#plotfit(dataX, FitX, S_curve, i,j )

# ╔═╡ 4b0eae99-e537-459d-a94b-9beb49a5ac4c
md"""
### ASF fit 
"""

# ╔═╡ 10162b29-ed32-431f-99d0-d302f3f4bd59
md"""
### PELT fit 
"""

# ╔═╡ 64154de3-6523-444b-aa5f-6ada2ae9a9e0
md"""
## Comparison between ASF and PELT/CROP
"""

# ╔═╡ 5349482c-ba2a-4a85-ae2c-d6b20c8ae6e1
begin
    @bind fitPELT CheckBox(false)
end

# ╔═╡ 664b3e4a-531b-4215-9909-242358d85d15
begin
    @bind fitEXP CheckBox(false)
end

# ╔═╡ cd836ae2-25dc-433d-8a06-fea4dcfbe861
begin
    @bind fitASF3 CheckBox(false)
end

# ╔═╡ cf6e8f34-4a1d-4ea3-8b7f-e569dac6c2c3
if fitEXP
	md"""
	## Fits with an exponential function 
	"""
end


# ╔═╡ 5c598f18-17dc-49cf-95b0-0c555fdb17c1
if fitEXP
	dfe, VI, VJ, VDX, VFX = fit_traces_exponential(gTRZS, gpeaks)
	dfe
end

# ╔═╡ ab210c22-68da-4a5b-b3e0-81d05da44c55
#if fitEXP
#	dfne = select_nexp(A_thresh = 100.0, tau_thresh = 35.0, chi2_thresh = 50.0)
#	md"""
#	- number of bad fits = $(size(dfne)[1])
#	"""
#end

# ╔═╡ 01f9827d-cef9-44db-a187-90918e70ea6f
#if fitEXP
#	plot_fits2(dfne, VDX, VFX, VI, VJ; plotsel="5x5")
#end

# ╔═╡ 9dc7b946-6213-475f-ae6d-ebead0844463
if fitPELT
	md"""
	## Fits with PELT GREEN
	"""
end

# ╔═╡ f8433824-dd11-4a8d-ae98-8d323926b72f
if fitASF3
	niter2 = 10
	md"""
	- The number of iterations defines the maximum number of steps that can be found in the data. A short number may underfit while a large number may overfit the data. Here we take: niter = $(niter2): 
"""
end

# ╔═╡ 2cf8cb42-61fe-4ea9-8047-7004511eaa19
if fitPELT
	#pdf, pI, pJ, pDX, pFX, pSC  =find_fit_pelt(gTRZS, gpeaks)
	pdf, pI, pJ, pDX, pFX  = jn.BoldLab.fit_pelt(gTRZS, gpeaks, ic=4) 
	pdf
end

# ╔═╡ 67ffc85b-8f2e-4f6c-8a80-3ba2c441bcc3
if fitPELT
	md"""
	## Fits with PELT RED
	"""
end

# ╔═╡ f14a9a85-414d-4b0e-8f84-569b74c5953c
#if fitPELT
#	rpdf, rpI, rpJ, rpDX, rpFX, rpSC  =find_fit_pelt(rTRZS, rpeaks)
#	pdf
#end

# ╔═╡ deddcc3c-f4ed-4e99-94db-bc59eb0446f4
#if fitPELT
#	plotsdf_pelt(rpdf, stept=200.0, steph=200.0, stepl=200.0)
#end

# ╔═╡ 7ddb353f-ef75-4ffd-9273-b02a3321474c
#if fitPELT
#	ppng = png_from_path(path_tiff)
#	plotsdf_pelt(pdf; ppath=ppng)
#end

# ╔═╡ d9d972f7-e196-4232-a062-c963c9147801
#md"""
#- Saving plots to $(ppng)
#"""

# ╔═╡ b9e5ce1f-4a86-4377-9d50-6732fae9516c
#readdir("stepAnalysis")

# ╔═╡ 1be93848-5757-48e6-8d5b-638cb11c4a61
md"""
## Functions
"""

# ╔═╡ de5a1406-000a-4acb-8e70-051c0aa76429
function select_trace(TRZS, df::DataFrame, row::Int) 
	i = df.i[row]
	j = df.j[row]
    trace = TRZS[i, j]  # safe, only defined where needed
	i,j, trace
end

# ╔═╡ 18688718-f034-463b-9f79-4997fe44c5b6
begin
	i,j, tz = select_trace(gTRZS, gpeaks, ntrz)
	md"""
	- selected trace = $(ntrz)
	"""
end

# ╔═╡ 56a220ee-adc1-42c2-b62d-6c520456c219
afsFit = jn.BoldLab.fit_asf(tz; niter=niter, tresH=thr)

# ╔═╡ c12cd5f1-6c56-400a-918a-ef0cfa842af0
md"""
#### Results from ASF: 
- number of steps = $(length(afsFit.steps))
- number of iterations = $(afsFit.iter)
- best shot = $(afsFit.best_shot)
- chi2/dof = $(afsFit.χ2_red)
"""

# ╔═╡ cc3653ee-d258-4928-a59c-22534d5d3d8c
p1 = jn.BoldLab.plot_asf_fit(tz, afsFit; show_residuals=true) 

# ╔═╡ 6c84ee94-ff04-4d85-8093-a1df36570db3
expFit = fit_traces_exponential(tz)

# ╔═╡ aeb6d87e-e115-4be0-9a55-8561592b8108
expFit

# ╔═╡ 00468149-d5ea-4ac7-962a-3b63e4566154
md"""
#### Results from Exponential Fit: 
- A0 = $(length(expFit.A0))
- tau = $(expFit.tau)
- chi2/dof = $(expFit.χ2_red)
"""

# ╔═╡ 219f2ef0-7bf4-4f15-a8fb-cd2893befd09
p2 = plot_exponential_fit(tz, expFit, show_residuals=true) 

# ╔═╡ 0b2ae59f-4eef-48cf-8daa-42f0b896c8cc
begin
	cro, si, sf = jn.BoldLab.pelt_fit(tz)
end

# ╔═╡ b1ddb27e-2a44-4f9d-b80a-23270ce30ffc
begin
	cost = cro["constrained"]
	steps = cro["number"]
	ielb = jn.BoldLab.elbow(steps, cost, method=:kneedle) 
	sc = scatter(steps, cost, xlabel="Number of Changepoints",
             ylabel="cost", 
             title="CROPS elbow", label="cost data")
	sc2=scatter!(sc,[steps[ielb]],[cost[ielb]], color=:red,
             marker=:circle,
             markersize=3, label="elbow")
	plot!(sc2, steps, cost, linewidth=1, label="")
	#jn.BoldLab.plot_crops_results(cro::Dict; title="CROPS Results")
end

# ╔═╡ e901d4c8-368c-4761-9b12-c1c88a16c0e5
md"""
#### PELT fit with ic =$(ielb)
"""

# ╔═╡ 236a3d2e-35e8-4a7f-9f81-2c160c719cac
length(cro)

# ╔═╡ d3f8806a-8d6c-4c0a-8220-d227ffa74bd2
begin
	pfit1  = jn.BoldLab.fit_pelt(tz, cro, ielb; si=si, sf=sf)
	pelt1 = jn.BoldLab.plot_pelt_fit(tz, pfit1, show_residuals=true) 
end

# ╔═╡ 5b8a48bf-8551-4dac-a58b-571e9c07f0ea
md"""
#### Results from PELT with IC=$(ielb): 
- number of steps = $(pfit1.nstep)
- ncost = $(pfit1.cost)
- penalty = $(pfit1.penalty)
- background region:  si = $(pfit1.si), sf = $(pfit1.sf)
- chi2/dof = $(pfit1.χ2_red)
"""

# ╔═╡ be2083ad-f167-419d-a674-4ecee527b7e8
begin
	pfit2  = jn.BoldLab.fit_pelt(tz, cro, 5; si=si, sf=sf)
	pelt2 = jn.BoldLab.plot_pelt_fit(tz, pfit2, show_residuals=true) 
end

# ╔═╡ 266ff2da-2901-4431-8268-26ec23ddc483
md"""
#### Results from PELT with IC=5: 
- number of steps = $(pfit2.nstep)
- ncost = $(pfit2.cost)
- penalty = $(pfit2.penalty)
- background region:  si = $(pfit2.si), sf = $(pfit2.sf)
- chi2/dof = $(pfit2.χ2_red)
"""

# ╔═╡ 1a1412b9-65b3-482e-b7c6-0da464635dec
function compare_asf_pelt(;ntx=100)
	
	ASFnsteps = Int[]
	ASFchi2   = Float64[]
	PLTnsteps = Int[]
	PLTchi2   = Float64[]
	PLTelb   = Int[]
	
	for ntrz in 1:ntx
		i,j, tz = select_trace(gTRZS, gpeaks, ntrz)
		
		afsFit = jn.BoldLab.fit_asf(tz; niter=niter, tresH=thr)
		if afsFit == nothing
			println("AFS fit failed for track number ", ntrz)
			continue
		end
		
		push!(ASFnsteps, length(afsFit.steps))
		push!(ASFchi2, afsFit.χ2_red)

		cro, si, sf = jn.BoldLab.pelt_fit(tz)

		if cro == nothing
			println("PELT fit failed to compute crops for track number ", ntrz)
			continue
		end
		
		cost = cro["constrained"]
		steps = cro["number"]
		ielb = jn.BoldLab.elbow(steps, cost, method=:kneedle)
		pfit1 = jn.BoldLab.fit_pelt(tz, cro, ielb; si=si, sf=sf)
		pfit1 = jn.BoldLab.fit_pelt(tz, cro, 5; si=si, sf=sf)

		if pfit1 == nothing
			println("PELT fit failed in fit_pelt for track number ", ntrz)
			continue
		end
		
		push!(PLTnsteps, pfit1.nstep)
		push!(PLTchi2, pfit1.χ2_red)
		push!(PLTelb, ielb)
	end
	return ASFnsteps, ASFchi2, PLTnsteps, PLTchi2, PLTelb
	
end

# ╔═╡ b23997a3-c47f-44d8-86a5-14d043eeb5a3
ASFnsteps, ASFchi2, PLTnsteps, PLTchi2, PLTelb = compare_asf_pelt(ntx=200)

# ╔═╡ c2438c47-31cd-4331-a2b2-35e96957a65d
begin
	hstepsASF, pstepsASF = step_hist(ASFnsteps;
	              nbins=20,
	              xlim=(0.0, 25.0),
	              logy=false,
	              xlabel="# steps",
	              ylabel="#entries ",
				  title="ASF nsteps ")
	
	hstepsPLT, pstepsPLT = step_hist(PLTnsteps;
	              nbins=20,
	              xlim=(0.0, 25.0),
	              logy=false,
	              xlabel="# steps",
	              ylabel="#entries ",
				  title="PLT nsteps ")
	
	hchi2ASF, pchi2ASF = step_hist(ASFchi2;
	              nbins=20,
	              xlim=(0.0, 200.0),
	              logy=false,
	              xlabel="chi2/dof",
	              ylabel="#entries ",
				  title="ASF chi2 ")
	
	hchi2PLT, pchi2PLT = step_hist(PLTchi2;
	              nbins=20,
	              xlim=(0.0, 200.0),
	              logy=false,
	              xlabel="chi2/dof",
	              ylabel="#entries ",
				  title="PLT chi2 ")

	plot(pstepsASF, pstepsPLT, pchi2ASF, pchi2PLT, layout=(2,2), size=(800, 400))
	
#histogram(peaks.intensity, bins=20)
end

# ╔═╡ 109e210b-5c0d-4ef7-97eb-61f48f339a8b
begin
	helbow, pelbow = step_hist(PLTelb;
	              nbins=10,
	              xlim=(0.0, 20.0),
	              logy=false,
	              xlabel="index elbow",
	              ylabel="#entries ",
				  title="Index Elbow ")
	plot(pelbow)
	
end

# ╔═╡ 80068e3e-e9c4-47c5-b6cf-e10ff2de19ea
function plot_peaks(imgd, peaks, nframe; colormap=:viridis)
	# Plot the grayscale image
	p1 = heatmap(gimgd[:, :, nframe]; 
			color = cgrad(:grays, rev = true),
			#color = :grays, 
			aspect_ratio = 1, title = "Detected Peaks")

	if show_peaks
		scatter!(p1, 
		    gpeaks.j,               # x-axis (columns)
		    gpeaks.i,               # y-axis (rows)
		    zcolor = gpeaks.intensity,
		    colorbar = true,
		    marker = (:circle, 3),
		    c = colormap,
		    label = "Peaks"
		)
	else
		p1
	end
end

# ╔═╡ eee7f7b1-b8fd-4f54-bbfb-122195d6bd5f
plot_peaks(gimgd, gpeaks, nframe; colormap=:viridis)

# ╔═╡ 920f1c33-2c46-4158-b88f-a9bef5590e72
plot_peaks(rimgd, rpeaks, nframe, colormap=:hot)

# ╔═╡ d25301bb-0fa2-4648-a76c-fc7d7218cbde
function select_nexp(; A_thresh = 100.0, tau_thresh = 35.0, chi2_thresh = 50.0)
	filter(row -> row.A0 > A_thresh || row.tau > tau_thresh || row.chi2_red > chi2_thresh, dfe)
end

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

# ╔═╡ c88c0cac-0ced-40c2-a472-a32e0e6ea20f
if fitEXP
	plot_fits(VDX, VFX, VI, VJ; plotsel="5x5")
end

# ╔═╡ 588add76-c9c8-4a21-a074-374e4f5e25ff
if fitPELT
	plot_fits(pDX, pFX, pI, pJ; plotsel="5x5")
end

# ╔═╡ d279c0a8-bdd4-408b-9132-399a9093f5a1
function plot_fits2(DF, VDX, VFX, VI, VJ; plotsel="5x5", figsize=(1500,1500))
    PP = []

    # Determine layout
    ly = (5,5)
    jl = 25
    if plotsel == "3x3"
        ly = (3,3); jl = 9
    elseif plotsel == "4x4"
        ly = (4,4); jl = 16
    end

    # Extract (i,j) pairs from DF
    target_pairs = Set([(row.i, row.j) for row in eachrow(DF)])

    # Match only entries in VDX/VFX that correspond to selected (i,j) in DF
    selected = [(k, i, j) for (k, (i,j)) in enumerate(zip(VI, VJ)) if (i, j) in target_pairs]

    # Limit to jl entries
    for (k, i, j) in Iterators.take(selected, jl)
        push!(PP, pltf(VDX[k], VFX[k], i, j))
    end

    return plot(PP..., layout=ly, size=figsize)
end

# ╔═╡ d4a11884-b670-4ceb-9444-36d3bcc6e8a7
function plotsdf_pelt(sdf; stept=200.0, steph=500.0, stepl=100.0, ppath="")

		
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
              xlim=(0.0, stept),
              logy=false,
              xlabel=" step time",
              ylabel=" #entries ")

	hstph, pstph = step_hist(sdf.stepHeight;
              nbins=25,
              xlim=(0.0, steph),
              logy=false,
              xlabel=" # step height",
              ylabel=" #entries ")


	hstplx, pstplx = step_hist(sdf.stepLength;
              nbins=50,
              xlim=(0.0, stepl),
              logy=false,
              xlabel=" # step length",
              ylabel=" #entries ")
	
	plt = plot(pnstep, pstptime, pstph, pstplx, 
		 layout=(2, 2), size=(800,600))
	
	#savefig(plt, ppath)
	plt
	
	
end

# ╔═╡ 7cfbc9b3-965b-45ca-bc33-f3a5f57b2ea3
if fitPELT
	plotsdf_pelt(pdf, stept=200.0, steph=200.0, stepl=200.0)
end


# ╔═╡ e1e3bc5d-9064-4794-91fd-23dd9f757846
function plotsdf_exp(sdf; amax=500.0, taumax=100.0, cmax=75.0, chi2max=200.0, ppath="")

	
	ha0, pa0 = step_hist(sdf.A0;
              nbins=25,
              xlim=(0.0, amax),
              logy=false,
              xlabel=" # A0",
              ylabel=" #entries ")

	
	htau, ptau = step_hist(sdf.tau;
              nbins=25,
              xlim=(0.0, taumax),
              logy=false,
              xlabel=" tau",
              ylabel=" #entries ")

	hc0, pc0 = step_hist(sdf.C0;
              nbins=25,
              xlim=(0.0, cmax),
              logy=false,
              xlabel=" # C0",
              ylabel=" #entries ")


	hchi2, pchi2 = step_hist(sdf.chi2_red;
              nbins=25,
              xlim=(0.0, chi2max),
              logy=false,
              xlabel=" # chi2",
              ylabel=" #entries ")
	
	plt = plot(pa0, ptau, pc0, pchi2, 
		 layout=(2, 2), size=(800,600))
	
	#savefig(plt, ppath)
	plt
	
	
end

# ╔═╡ d68e0b2c-25e4-4b79-be90-fb1c3149f33d
if fitEXP
	plotsdf_exp(dfe, amax=200.0, taumax=100.0, cmax=200.0, chi2max=200.0)
end

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

# ╔═╡ 300fcee3-17cd-44f4-bd76-6c8bc03b30a4
if fitPELT
md"""
### PELT, GREEN
- number of fitted molecules = $(count_mol(pdf))
"""
end

# ╔═╡ ab730d0c-0e11-4aa1-acb5-48dcb250f97b
if fitPELT
md"""
### PELT, RED
- number of fitted molecules = $(count_mol(rpdf))
"""
end

# ╔═╡ Cell order:
# ╠═9292025e-0d7d-11f0-365e-f1724fc39b3c
# ╠═c5d3b3b0-d5f2-4102-8bc7-e6af113e8e63
# ╠═5b3b5f83-058a-4582-9482-9b2fc5892692
# ╠═b9ef4450-d619-4e67-9ff3-9c8045b3863d
# ╠═f8573e77-38d8-467c-a123-9a83f73e8970
# ╠═19476521-78bd-44d6-a302-6e971f9fc095
# ╠═cba7fc75-8363-4ffa-b5a6-6e7d34363813
# ╠═62925080-87f3-4c84-8fc4-1b9ae36ec7c5
# ╠═6586349e-37f9-4cd2-8513-d05558ec3d6b
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
# ╠═967ef731-1cee-4091-8f98-4a7f080058be
# ╠═555f27e7-a820-4830-b47a-7bb50458cdc5
# ╠═82e619dd-27d1-41b2-81e5-d102b972b3d8
# ╠═f724b5fc-b506-4442-b6b0-92b8fc9ad16b
# ╠═95804d28-2720-4a5a-8003-05b1bc86f00c
# ╠═89ac3496-c36a-4472-a21c-bf4fca5e79d6
# ╠═e3ce6b68-4186-41a5-a3f1-a458b78d8821
# ╠═62c9f533-a314-4c80-b700-8b83bc3c8f83
# ╠═c8770df9-727b-4095-b45c-d9d10877227e
# ╠═04601060-8310-4e76-a2d0-8e258d2d061f
# ╠═9c6db6db-b5b7-49e6-b275-970a29ba69cb
# ╠═91b5254a-b21a-4243-a838-811ddc90d515
# ╠═5a63ed41-9233-4654-a37c-d00013c35197
# ╠═797aaf3e-f284-4f62-bb1a-5514190ca422
# ╠═4e4252bc-3b6e-4507-8817-dfe9558b51b6
# ╠═783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
# ╠═7e40b817-a1d3-45c6-b651-136031aad7b3
# ╠═ee8175fb-f6dd-46fc-baed-ae82209af0ed
# ╠═eb6b49be-ec65-447d-afe6-2c7580bd6693
# ╠═11a04a68-91d0-4595-9296-144b3e21f478
# ╠═d2a8345c-9c44-43d7-a9ec-6ee4ab58674d
# ╠═b0a4c0be-12b9-416c-a45c-5fe35fbbd168
# ╠═69734f52-a953-44e3-81e6-a8f911c4f4d0
# ╠═7c267c22-6d75-46a4-a459-6a03d74da0cb
# ╠═89354633-e220-45f2-989e-c5575acd2988
# ╠═5fcf28e8-0467-4f6f-b2e7-676693e8355e
# ╠═83b89e38-10e8-4ce1-9639-29b3afb38a64
# ╠═f4ff7938-f16e-4159-95e8-cb98c59a9d80
# ╠═156d523f-1b80-4cbd-912c-ed21a1327aff
# ╠═fd39234f-6c62-42e3-bde9-9a9e565fa519
# ╠═62da10f6-0d5c-41c0-a985-d15c946f5b84
# ╠═12960d51-135c-47cd-ab86-f2ab5bacef08
# ╠═d34b9105-09d4-4876-842d-bcf74249cca9
# ╠═4d0c63a1-89e0-4d89-846c-e1501dbc2696
# ╠═eee7f7b1-b8fd-4f54-bbfb-122195d6bd5f
# ╠═a74eff14-b6b4-4e80-bfc7-5a1db6a32aba
# ╠═920f1c33-2c46-4158-b88f-a9bef5590e72
# ╠═7927a18c-6942-42bd-ac6b-2ff720b155d0
# ╠═1333c10b-c751-414b-9d71-df9f6499172b
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
# ╠═574e8729-7610-47c0-be0a-767f4af15c7c
# ╠═77f5c1c3-11bb-4976-af44-1b488f08de6b
# ╠═56e7b887-3304-47e4-8c17-8268db12a46b
# ╠═b24005f5-e019-4f41-9193-ebf7e5d69dac
# ╠═381bb69d-1b27-4403-9f95-bfde03fc1c15
# ╠═8d46c82c-6585-4d79-b5a9-2ee51d321e36
# ╠═a8b6fe7c-f328-4602-b7e2-240f53b78b64
# ╠═18688718-f034-463b-9f79-4997fe44c5b6
# ╠═c11c894f-cefa-4d50-9c27-19cb744ef7ef
# ╠═709275a8-af9f-4cea-9e9f-1c22c95cf48c
# ╠═073db006-22d0-4d27-8148-ef142b163ba2
# ╠═4b0eae99-e537-459d-a94b-9beb49a5ac4c
# ╠═56a220ee-adc1-42c2-b62d-6c520456c219
# ╠═cc3653ee-d258-4928-a59c-22534d5d3d8c
# ╠═c12cd5f1-6c56-400a-918a-ef0cfa842af0
# ╠═6c84ee94-ff04-4d85-8093-a1df36570db3
# ╠═219f2ef0-7bf4-4f15-a8fb-cd2893befd09
# ╠═aeb6d87e-e115-4be0-9a55-8561592b8108
# ╠═00468149-d5ea-4ac7-962a-3b63e4566154
# ╠═10162b29-ed32-431f-99d0-d302f3f4bd59
# ╠═0b2ae59f-4eef-48cf-8daa-42f0b896c8cc
# ╠═b1ddb27e-2a44-4f9d-b80a-23270ce30ffc
# ╠═e901d4c8-368c-4761-9b12-c1c88a16c0e5
# ╠═236a3d2e-35e8-4a7f-9f81-2c160c719cac
# ╠═d3f8806a-8d6c-4c0a-8220-d227ffa74bd2
# ╠═5b8a48bf-8551-4dac-a58b-571e9c07f0ea
# ╠═be2083ad-f167-419d-a674-4ecee527b7e8
# ╠═266ff2da-2901-4431-8268-26ec23ddc483
# ╠═64154de3-6523-444b-aa5f-6ada2ae9a9e0
# ╠═1a1412b9-65b3-482e-b7c6-0da464635dec
# ╠═b23997a3-c47f-44d8-86a5-14d043eeb5a3
# ╠═c2438c47-31cd-4331-a2b2-35e96957a65d
# ╠═109e210b-5c0d-4ef7-97eb-61f48f339a8b
# ╠═5349482c-ba2a-4a85-ae2c-d6b20c8ae6e1
# ╠═664b3e4a-531b-4215-9909-242358d85d15
# ╠═cd836ae2-25dc-433d-8a06-fea4dcfbe861
# ╠═cf6e8f34-4a1d-4ea3-8b7f-e569dac6c2c3
# ╠═5c598f18-17dc-49cf-95b0-0c555fdb17c1
# ╠═c88c0cac-0ced-40c2-a472-a32e0e6ea20f
# ╠═d68e0b2c-25e4-4b79-be90-fb1c3149f33d
# ╠═ab210c22-68da-4a5b-b3e0-81d05da44c55
# ╠═01f9827d-cef9-44db-a187-90918e70ea6f
# ╠═9dc7b946-6213-475f-ae6d-ebead0844463
# ╠═f8433824-dd11-4a8d-ae98-8d323926b72f
# ╠═2cf8cb42-61fe-4ea9-8047-7004511eaa19
# ╠═588add76-c9c8-4a21-a074-374e4f5e25ff
# ╠═7cfbc9b3-965b-45ca-bc33-f3a5f57b2ea3
# ╠═300fcee3-17cd-44f4-bd76-6c8bc03b30a4
# ╠═67ffc85b-8f2e-4f6c-8a80-3ba2c441bcc3
# ╠═f14a9a85-414d-4b0e-8f84-569b74c5953c
# ╠═deddcc3c-f4ed-4e99-94db-bc59eb0446f4
# ╠═ab730d0c-0e11-4aa1-acb5-48dcb250f97b
# ╠═7ddb353f-ef75-4ffd-9273-b02a3321474c
# ╠═d9d972f7-e196-4232-a062-c963c9147801
# ╠═b9e5ce1f-4a86-4377-9d50-6732fae9516c
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═de5a1406-000a-4acb-8e70-051c0aa76429
# ╠═80068e3e-e9c4-47c5-b6cf-e10ff2de19ea
# ╠═d25301bb-0fa2-4648-a76c-fc7d7218cbde
# ╠═862d9a0c-ecc6-4178-80f4-74fc71a79f38
# ╠═d279c0a8-bdd4-408b-9132-399a9093f5a1
# ╠═d9477e2f-2bad-4243-a11f-393f5b3405c7
# ╠═efd033eb-bde9-4f4e-85f0-777125555edd
# ╠═fa0da6ef-68ba-4db2-aec0-26fd3a280b5a
# ╠═d4a11884-b670-4ceb-9444-36d3bcc6e8a7
# ╠═e1e3bc5d-9064-4794-91fd-23dd9f757846
# ╠═f9c4fe44-b596-4d54-a07a-47a9ca76c108
# ╠═81181098-822a-48f9-b776-680735de6430
# ╠═4b7e4552-ffb9-4f9a-a25f-f45178b2fbd3
# ╠═f0aab374-2e76-47e0-a62c-a5ec21129767
# ╠═a3c2e539-eda7-4784-b90c-43fcf6982347
# ╠═83df15d6-5493-48ef-a8f4-29b8fc375da3
