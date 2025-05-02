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
	using Images, FileIO, ImageIO
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


# ╔═╡ a0c20599-dfbf-4614-86fb-0c4a768f33dc
md"""
# Multiple pass versus single pass

- Results show that multiple pass does not improve fit in the Monte Carlo 
- To reproduce results, be sure to select "MC" and "n3d_mc_1e3_noise.py"
"""

# ╔═╡ 550cb0b0-69b5-4933-a494-36044e4beb9d
md"""
- The MC data corresponds to 10^3 molecules distributed in ~10^4 pixels, so the density of molecules expected is small.
"""

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
    BASEDIRS = Dict("Data"=>"/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign", 
					"MonteCarlo" => "/Users/jjgomezcadenas/Projects/BoldLab/pluto/npy"     
    )
    base_directory(label) = BASEDIRS[label]
end

# ╔═╡ fe676ab1-4ec5-4c3b-beb1-62c68f4f87ee
begin
    function scan_level(dir::AbstractString)
        entries = readdir(dir; join=true, sort=true)
        vis     = filter(e -> !startswith(basename(e), "."), entries)

        subdirs = filter(isdir, vis)
        npys    = filter(e -> endswith(e, ".npy"), vis)
        tiffs   = filter(e -> endswith(lowercase(e), ".tif")  ||
                               endswith(lowercase(e), ".tiff"), vis)

        return basename.(subdirs), basename.(npys), basename.(tiffs)
    end
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

# ╔═╡ 2fe0ccbc-0b38-415c-995e-987017bcb499
#imx1, CI = regularize_img(imxt[:,:,1])

# ╔═╡ af3e3cde-368d-49ad-a315-e83e9414a263
plot_frames(imxt; nscale=20)

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

# ╔═╡ 053dd59b-f6cb-47b2-9bfe-ae0c7f2a514b
typeof(imxt)

# ╔═╡ 57bd56be-5729-44e8-aba3-783a78c714d2
begin
	TRZ = build_traces(imxt)
	md"""
	- Compute traces:  $(size(TRZ)) 
	"""
end

# ╔═╡ 6b54a654-7850-4a27-8ec4-2dcbd3566d82
#hsum2d, hmean2d = traces_h2d(xsum, xmean, xstd; bins=50)

# ╔═╡ 7b89b967-0b2e-4f8d-8b8f-c6597a20a638
TRZ

# ╔═╡ 3a8e7300-d209-4e50-ace8-c4b841d71f42
md"""
#### Plot traces
"""

# ╔═╡ a7542de4-5e78-4c55-91b7-d68c4e2f7985
begin
	i = 100
	j = 200
	nx,mx = size(TRZ)
	if i > nx
		warn("i = $(i) is larger, than $(nx) setting i = $(nx)")
		i = nx 
	end

	if j > mx
		warn("j = $(j) is larger, than $(mx) setting j = $(mx)")
		j = mx
	end
	
	md"""
	- Fitting trace in pixel ($(i), $(j))
	- Intensity in pixel = $(sum(TRZ[i,j]))
	"""
end

# ╔═╡ d1c82dbe-8493-464d-bdac-f505657678d6
begin
	
	dataX, FitX, S_curve, bshot = fit_traces(TRZ, i, j; niter=3, thr=0.5)
	sth, stt, stl = getsteps(FitX)

	

	md"""
	- Fit results
	- bes shot =$(bshot)
	- Step heights: = $(vect_to_fstr(sth, "%.2f"))
	- Step times =$(vect_to_fstr(stt, "%.2f"))
	- Segment lengths = $(vect_to_fstr(stl, "%.2f"))
	"""
	
end

# ╔═╡ 8aa50fbb-74a9-424d-9a84-e717247c65a9
plotfit(dataX, FitX, S_curve, i,j )

# ╔═╡ 77f5c1c3-11bb-4976-af44-1b488f08de6b
md"""
## Fits to data
"""

# ╔═╡ 3115a147-d41b-4ab6-9ad9-0f3b30fb4504
md"""
- Set the vaue of threshold
"""

# ╔═╡ e04c9c53-7ae2-474b-87bc-97bd58f458fa
@bind thr NumberField(0.0:0.1:1.1, default=0.5)

# ╔═╡ f175508e-f0ea-48a9-ba72-d2e21215de1d
md"""
- Set the vaue of iterations (number of steps to be sought in the data)
"""

# ╔═╡ 86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
@bind niter NumberField(0:10, default=3)

# ╔═╡ 264fcbb0-ec40-42b3-ab05-01b2171da6f2
begin
	md"""
	- number of iterations (=steps to fit) = $(niter)
	- threshold = $(thr)
	"""
end

# ╔═╡ 55c4a936-659b-4adf-aaac-3817741e6732
pedx

# ╔═╡ 0c3d3a1e-96be-4654-8ff0-b8f03d6c5862
typeof(TRZ)

# ╔═╡ 08dd1429-4c2d-4def-8282-e7abe469f318
md"""
### Fit data
"""

# ╔═╡ 85a3b5cd-9130-475f-acdf-c0eb20ee8a6a
begin
	#DX, FX, TL, sdf, npixel, nfailed = fit_data3(TRZ; ped=1600.0, tt = 0, niter=10, thr=0.05, fplot=2500)

	asdf, anpixel, angfit, anfailed, aII, aJJ, aDX, aFX, aITER, aCC = fit_data3(TRZ; ped=pedx, tt = 0, niter=niter, thr=thr, stdf=0.0, sel="auto")
	md"""
	- Fitted $(anpixel) pixels 
	- failed fit fraction =$(anfailed/anpixel)
	- Found $(length(nonzeros(aDX))) good fits
	"""
end

# ╔═╡ 5a13e1e1-8679-4c6c-85d6-081e5abb9e57
begin
	#DX, FX, TL, sdf, npixel, nfailed = fit_data3(TRZ; ped=1600.0, tt = 0, niter=10, thr=0.05, fplot=2500)

	bsdf, bnpixel, bngfit, bnfailed, bII, bJJ, bDX, bFX, bITER, bCC = fit_data3(TRZ; ped=pedx, tt = 0, niter=niter, thr=thr, stdf=0.0, sel="base")
	md"""
	- Fitted $(bnpixel) pixels 
	- failed fit fraction =$(bnfailed/bnpixel)
	- Found $(length(nonzeros(bDX))) good fits
	"""
end

# ╔═╡ 177f250e-3bae-4f2b-a980-ca5357fc81bf
begin
	hiter, piter = step_hist(aITER;
	              nbins=12,
	              xlim=(-12.0, 12.0),
	              logy=false,
	              xlabel=" iterations",
	              ylabel=" #entries ")
	plot(piter, size=(600,300))
end

# ╔═╡ 2664248e-dfee-4c06-aa73-55f75c9f122d
begin
	hcc, pcc = step_hist(aCC;
	              nbins=12,
	              xlim=(-5.0, 5.0),
	              logy=false,
	              xlabel=" CC",
	              ylabel=" #entries ")
	plot(pcc, size=(600,300))
end

# ╔═╡ 3ed32697-10d4-4b76-8540-d4431e5a3fc0
unique(asdf.nstep)

# ╔═╡ a40752df-cfa7-40fa-84b0-8e380e69a924
unique(bsdf.nstep)

# ╔═╡ 0db1ec72-9f82-4e6c-928e-f7ef7ced0966
 plot_trx(aDX, aFX; ni=1, nf=16, layout=(4,4), figsize=(1500, 1500))

# ╔═╡ ea5d4c44-c9cb-4603-8ad5-0d6070779e5d
 plot_trx(bDX, bFX; ni=1, nf=16, layout=(4,4), figsize=(1500, 1500))

# ╔═╡ 7a68a368-5493-4090-ad6a-b9a7c3b10cb5
asdf

# ╔═╡ a760e4c4-e067-4385-b449-c2f897042839
anmx = unique(asdf.nmol)

# ╔═╡ 2b6e3ece-6a39-496b-ae78-a2fa6d2edd50
length(anmx)

# ╔═╡ c1ed3b0f-1c92-4a11-b424-042551275f8f
bnmx = unique(bsdf.nmol)

# ╔═╡ 26c4cfe9-1596-414a-aa31-a16158d8910b
length(bnmx)

# ╔═╡ 3e9fd7f8-56c0-4b3d-a422-7a8faa7a8245
plot_traces(asdf, anmx[1:16], lyt=(4,4), size=(1500, 1500))

# ╔═╡ c88f5014-7abb-4b50-8a52-2533cda749f5
plot_traces(bsdf, bnmx[1:16], lyt=(4,4), size=(1500, 1500))

# ╔═╡ 0dab9ef6-0cdb-4983-b87a-bf5d41ae9661
mstl = maximum(bsdf.stepLength)

# ╔═╡ 16dbb6bc-5d8f-43a7-8063-caaba9c4e71b
md"""
- Cut away flat cases (nstep=1, stepLength=mstl)
"""

# ╔═╡ aa40b9e7-0992-4228-8e93-a7ad91e7f8c6
#sdf2 = filter(row -> 
#	    row.stepTime > 1 && row.stepTime < mstl-20  && row.stepLength >2 && row.stepLength <mstl-20,
#	  bsdf)
#	update_nstep!(sdf2)

# ╔═╡ 154277d2-cdfa-4f7f-8f07-63ddd4e90ae8


# ╔═╡ 1be93848-5757-48e6-8d5b-638cb11c4a61
md"""
## Functions
"""

# ╔═╡ cb7f11ec-e468-41f3-b9fc-c6f12a27d5b3
function get_vals_from_sparse(sm)
	rows = Int[]
	cols = Int[]
	vals = Float64[]

	for j in 1:size(sm, 2)
	    for idx in sm.colptr[j]:(sm.colptr[j+1]-1)
	        i = sm.rowval[idx]
	        v = sm.nzval[idx]
	        push!(rows, i)
	        push!(cols, j)
	        push!(vals, v)
	    end
	end
	return rows, cols, vals
end

# ╔═╡ facc11b5-f86a-4390-ac68-a3a74aa8d247
function reduced_chi_squared_per_element(DX::SparseMatrixCSC{<:AbstractVector, Int},
                                          FX::SparseMatrixCSC{<:AbstractVector, Int})

    @assert size(DX) == size(FX) "DX and FX must have the same size"

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for col in 1:size(DX, 2)
        for idx in DX.colptr[col]:(DX.colptr[col+1] - 1)
            row = DX.rowval[idx]
            dvec = DX[row, col]
            fvec = FX[row, col]

            @assert length(dvec) == length(fvec) "Vectors at ($row,$col) must be same length"

            chi2 = 0.0
            dof = 0

            for k in eachindex(dvec)
                d = dvec[k]
                f = fvec[k]
                if d > 0
                    chi2 += (d - f)^2 / d
					#chi2 += (d - f)^2 
                    dof += 1
                end
            end

            if dof > 0
                push!(rows, row)
                push!(cols, col)
                push!(vals, chi2 / dof)
            end
        end
    end

    return sparse(rows, cols, vals, size(DX, 1), size(DX, 2))
end

# ╔═╡ d0e0c953-6fd9-4f2b-b328-b344ccdbde44
begin 
	chi2_matrix = reduced_chi_squared_per_element(aDX, aFX)
	X2_dense = Array(chi2_matrix)
	Ic2, Jc2, Vc2 = get_vals_from_sparse(chi2_matrix)
	md"""
	### aFIT
	- mean chi2 = $(mean(Vc2))
	- std chi2 = $(std(Vc2))
	- max chi2 = $(maximum(Vc2))
	"""
end

# ╔═╡ 2f94d37c-f2a1-40b8-87cd-b7b97139d0f0
heatmap(X2_dense, colorbar_title="Reduced χ²", xlabel="j", ylabel="i", title="Reduced Chi-squared Map")

# ╔═╡ 7c0fd731-9a54-4eeb-a009-cb6b5fc438a4
length(Vc2)

# ╔═╡ 6b642fe5-5f55-4a9d-86cb-1695dca29ce2
begin
	hvc2, pvc2 = step_hist(Vc2;
	              nbins=20,
	              xlim=(0.0, 5.0),
	              logy=false,
	              xlabel=" chi2: FIT A",
	              ylabel=" #entries ")
	plot(pvc2, size=(600,300))
end

# ╔═╡ ffdadd86-0e63-4cdc-b770-b1693614f1e2
begin 
	bchi2_matrix = reduced_chi_squared_per_element(bDX, bFX)
	bX2_dense = Array(bchi2_matrix)
	Ic2b, Jc2b, Vc2b = get_vals_from_sparse(bchi2_matrix)
	md"""
	### bFIT
	- mean chi2 = $(mean(Vc2b))
	- std chi2 = $(std(Vc2b))
	- max chi2 = $(maximum(Vc2b))
	"""
end

# ╔═╡ 3eecf083-f042-4514-ad6c-5bf265aa0607
begin
	hvc2b, pvc2b = step_hist(Vc2b;
	              nbins=20,
	              xlim=(0.0, 5.0),
	              logy=false,
	              xlabel=" chi2: FIT B",
	              ylabel=" #entries ")
	plot(pvc2b, size=(600,300))
end

# ╔═╡ 69fa05aa-8208-414c-9b92-a4a1f7be1c68
function reduced_fit(DX::SparseMatrixCSC{<:AbstractVector, Int},
                                          FX::SparseMatrixCSC{<:AbstractVector, Int})

    @assert size(DX) == size(FX) "DX and FX must have the same size"

    rows = Int[]
    cols = Int[]
    dvals = Float64[]
	fvals = Float64[]

    for col in 1:size(DX, 2)
        for idx in DX.colptr[col]:(DX.colptr[col+1] - 1)
            row = DX.rowval[idx]
            dvec = DX[row, col]
            fvec = FX[row, col]

            @assert length(dvec) == length(fvec) "Vectors at ($row,$col) must be same length"

            fr = 0.0
            dr = 0.0

            for k in eachindex(dvec)
                dr  += dvec[k]
                fr += fvec[k]
            end
			
            push!(rows, row)
            push!(cols, col)
            push!(dvals, dr)
			push!(fvals, fr)
            
        end
    end

    return sparse(rows, cols, dvals, size(DX, 1), size(DX, 2)), sparse(rows, cols, fvals, size(FX, 1), size(FX, 2))
end

# ╔═╡ 8f6609a3-b24d-404e-90c7-0ced26436f0e
begin 
	rDX, rFX = reduced_fit(aDX, aFX)
	rdDX = Array(rDX)
	rdFX = Array(rFX)
	#Ix2, Jx2, Vx2 = get_vals_from_sparse(rdDX)
	#Iy2, Jy2, Vy2 = get_vals_from_sparse(rdFX)
	hmp1 = heatmap(rdDX, colorbar_title="DX", xlabel="j", ylabel="i", title="Reduced DX Map")
	hmp2 = heatmap(rdFX, colorbar_title="FX", xlabel="j", ylabel="i", title="Reduced FX Map")
	plot(hmp1, hmp2, layout=(2,1), size=(900, 900))
end

# ╔═╡ cb56c6c6-3d92-43af-8556-232f7387dc97
function cut_length(df; maxl=150)

	function is_bad_molecule(gdf)
    	any(gdf.stepLength .> maxl)
	end

	grouped = groupby(df, :nmol)
	
	# Keep only groups that do not satisfy the bad condition
	good_groups = filter(gdf -> !is_bad_molecule(gdf), grouped)
	
	# Reassemble the cleaned DataFrame
	vcat(good_groups...)
end


# ╔═╡ 7898bcc1-651d-4a61-aff5-d8f058bf74b4
function chi_squared(DX::AbstractVector, FX::AbstractVector, σ::AbstractVector)
    sum(((DX .- FX) ./ σ).^2)
end

# ╔═╡ d7f05829-b372-4ad0-afc4-645f1a37193f
function find_fit_candidates(trz, nc, sel; ped=1600, niter=5, thr=0.5)
	n, m = size(trz)
	I = Int[]
	J = Int[]
	BS = []
    DX = Vector{Vector{Float32}}()
    FX = Vector{Vector{Float32}}()
	ITER = Int[]
	CC = Float64[]
	ng = 0
	for i in 1:n
	    for j in 1:m
	        dataX, FitX, S_curve, best_shot, iter, cc = fit_traces(trz, i, j; 
	                                                     niter=niter,
															   thr=thr, sel=sel)
	
	        if best_shot >0
				ng+=1
				push!(I,i)
				push!(BS, best_shot)
				push!(J,j)
				push!(DX, dataX .-ped)
				push!(FX, FitX .- ped)
				push!(ITER, iter)
				push!(CC, cc)
				if ng > nc-1
					return I, J, BS, DX, FX, ITER, CC
				end
			end	
		end
	end
	
	I, J, BS, DX, FX, ITER, CC
end

# ╔═╡ 6e2b4709-a8df-47ad-9567-bcf39e020bc8
cI, cJ, cBS, cDX, cFX, cITER, cCC = find_fit_candidates(TRZ, 16,"base"; ped=pedx, niter=niter, thr=thr)

# ╔═╡ 0ca90db3-4700-4028-ae04-b29d397cdae7
axI, axJ, axBS, axDX, axFX, axITER, axCC = find_fit_candidates(TRZ, 16,"auto"; ped=pedx, niter=niter, thr=thr)

# ╔═╡ 15a5e0fc-9416-455e-ba93-fb8ea97a75fe
let
	PP =[]
	for i in 1:16
		push!(PP, plotfit2(cDX[i], cFX[i], axFX[i]))
	end
	plot(PP..., layout=(4, 4), size=(1500,1500))
#plotfit(dataX, FitX, S_curve, i, j)
end

# ╔═╡ df8347f2-be16-413f-9a93-e3baa9ed8639
function plot_sl(aFX, cFX)
	function get_histos(FX, lbl)
		nzfx = nonzeros(aFX)
		stdfx = [std(nzfx[i]) for i in 1:length(nzfx)]
		xtdfx = [mean(nzfx[i]) for i in 1:length(nzfx)]
		hx, phfx = step_hist(stdfx;
	              nbins=20,
	              xlim=(0.0, 10.0),
	              logy=false,
	              xlabel=" std step length ($(lbl))",
	              ylabel=" #entries ")
		h2x, phfx2 = step_hist(xtdfx;
	              nbins=20,
	              xlim=(0.0, 400.0),
	              logy=false,
	              xlabel=" mean step length ($(lbl))",
	              ylabel=" #entries ")
		return nzfx, stdfx, xtdfx, hx, phfx, h2x, phfx2
	end
	
	nzafx, stdafx, stdafx, hstdax, phstdax, hmeanax, phmeanax = get_histos(aFX, "aFX")
	nzcfx, stdcfx, stdcfx, hstdcx, phstdcx, hmeancx, phmeancx = get_histos(cFX, "cFX")
	plot(phstdax, phmeanax, phstdcx, phmeancx, 
		layout=(2,2), size=(1200, 1200))
end

# ╔═╡ 66b74aa7-8cd4-4e72-a9ec-9f6671e9c547
plot_sl(aFX, bFX)

# ╔═╡ f9dd1f85-da7f-41db-af12-d1b26ca36f99
function plot_step_time(asdf, csdf)
	function histos(df, lbl)
		hstx, pstx = step_hist(df.stepTime;
	              nbins=50,
	              xlim=(0.0, 400.0),
	              logy=false,
	              xlabel=" step time ($(lbl))",
	              ylabel=" #entries ")
		return hstx, pstx
	end
	ahstx, apstx = histos(asdf, "asdf")
	chstx, cpstx = histos(csdf, "csdf")
	plot(apstx, cpstx, layout =(1, 2), size=(800,400))
end

# ╔═╡ 6abd282d-aea2-4b9f-980c-9f9621f80442
plot_step_time(asdf, bsdf)

# ╔═╡ d4a11884-b670-4ceb-9444-36d3bcc6e8a7
function plotsdf(sdf)

		
	# take only one row per (i,j), keeping the first nstep
	df_pix = unique(sdf, [:i, :j, :nstep])

	
	hstpx, pstpx = step_hist(df_pix.nstep;
              nbins=10,
              xlim=(0.0, 10.0),
              logy=false,
              xlabel=" # steps",
              ylabel=" #entries ")
	
	hstx, pstx = step_hist(sdf.stepTime;
              nbins=50,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" step time",
              ylabel=" #entries ")

	hstphx, pstphx = step_hist(sdf.stepHeight;
              nbins=50,
              xlim=(0.0, 500.0),
              logy=false,
              xlabel=" # step height",
              ylabel=" #entries ")

	hstplx, pstplx = step_hist(sdf.stepLength;
              nbins=50,
              xlim=(0.0, 400.0),
              logy=false,
              xlabel=" # step length",
              ylabel=" #entries ")
	
	plot(pstphx, size=(400,200))
	
	plot(pstpx, pstx, pstphx, pstplx, layout=(2, 2), size=(800,400))
	
	
end

# ╔═╡ 0a8baaba-06a6-4d44-bbc7-437e740be7b2
plotsdf(asdf)

# ╔═╡ b0d9d17b-d421-42ea-91ae-cf530675a44f
plotsdf(bsdf)

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

# ╔═╡ 84023e48-6534-42a8-a4f4-23d242d2ca88
begin
	sdf2 = filter(row -> 
	    !(row.nstep >=3),
	  bsdf)
	update_nstep!(sdf2)
end

# ╔═╡ 907dccfa-0726-487f-91af-ca6f1c096d84
histogram(sdf2.nstep)

# ╔═╡ ed77f48f-0df0-4718-8862-673c01250dfe
begin
	nsdf2 = unique(sdf2.nmol)
	lsdf2 = length(nsdf2)

	md"""
	- number of fits = $(lsdf2)
	"""
end

# ╔═╡ 28a90470-6d12-4de7-9f43-f53681ca7427
plot_traces(sdf2, nsdf2[1:8], lyt=(3,3), size=(1500, 1500))

# ╔═╡ fe28cfc3-d0b2-4584-8abc-41b7bea8b1c1
plotsdf(sdf2)

# ╔═╡ e97a39c5-6d3a-4169-8143-4f02dece1e04
sdf3 = cut_length(sdf2; maxl=mstl-20)

# ╔═╡ 9980987c-ac9d-4786-b31e-54f9ee45d467
begin
	nsdf3 = unique(sdf3.nmol)
	lsdf3 = length(nsdf3)

	md"""
	- number of fits = $(lsdf3)
	"""
end

# ╔═╡ 75e778de-9c70-43bf-b5c4-a0d2c6b1a539
plot_traces(sdf3, nsdf3[1:16], lyt=(4,4), size=(1500, 1500))

# ╔═╡ 820beb62-fe16-494b-8e76-323ce30e46d0
plotsdf(sdf3)

# ╔═╡ 3b85b24d-34bb-4df7-a0c3-cb89f7d90b84
begin

	vnm = sdf3.nmol 
	
	# Get unique values
	unique_nm = unique(vnm)
	
	# Get count of each element
	counts_dict = countmap(vnm)
	
	# Convert dictionary to count vector in the same order as unique_vals
	vns = [counts_dict[val] for val in unique_nm]
	npixels = length(unique_nm)
	md"""
	- Data set has $(npixels) molecules
	"""
end

# ╔═╡ 34634f51-0ec5-4ef1-8b21-b18ee1c6ced0
if casedata
	ss = save_stats(sdf3, npixels, folder_scan, folder_field)
	Markdown.parse(ss)
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
# ╠═b9ef4450-d619-4e67-9ff3-9c8045b3863d
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
# ╠═f724b5fc-b506-4442-b6b0-92b8fc9ad16b
# ╠═a0c20599-dfbf-4614-86fb-0c4a768f33dc
# ╠═550cb0b0-69b5-4933-a494-36044e4beb9d
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
# ╠═2fe0ccbc-0b38-415c-995e-987017bcb499
# ╠═af3e3cde-368d-49ad-a315-e83e9414a263
# ╠═783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
# ╠═57692746-065d-4d75-8a41-7ffafd69550e
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═053dd59b-f6cb-47b2-9bfe-ae0c7f2a514b
# ╠═57bd56be-5729-44e8-aba3-783a78c714d2
# ╠═6b54a654-7850-4a27-8ec4-2dcbd3566d82
# ╠═7b89b967-0b2e-4f8d-8b8f-c6597a20a638
# ╠═3a8e7300-d209-4e50-ace8-c4b841d71f42
# ╠═a7542de4-5e78-4c55-91b7-d68c4e2f7985
# ╠═d1c82dbe-8493-464d-bdac-f505657678d6
# ╠═8aa50fbb-74a9-424d-9a84-e717247c65a9
# ╠═77f5c1c3-11bb-4976-af44-1b488f08de6b
# ╠═3115a147-d41b-4ab6-9ad9-0f3b30fb4504
# ╠═e04c9c53-7ae2-474b-87bc-97bd58f458fa
# ╠═f175508e-f0ea-48a9-ba72-d2e21215de1d
# ╠═86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
# ╠═264fcbb0-ec40-42b3-ab05-01b2171da6f2
# ╠═55c4a936-659b-4adf-aaac-3817741e6732
# ╠═0c3d3a1e-96be-4654-8ff0-b8f03d6c5862
# ╠═6e2b4709-a8df-47ad-9567-bcf39e020bc8
# ╠═0ca90db3-4700-4028-ae04-b29d397cdae7
# ╠═15a5e0fc-9416-455e-ba93-fb8ea97a75fe
# ╠═08dd1429-4c2d-4def-8282-e7abe469f318
# ╠═85a3b5cd-9130-475f-acdf-c0eb20ee8a6a
# ╠═8f6609a3-b24d-404e-90c7-0ced26436f0e
# ╠═d0e0c953-6fd9-4f2b-b328-b344ccdbde44
# ╠═2f94d37c-f2a1-40b8-87cd-b7b97139d0f0
# ╠═7c0fd731-9a54-4eeb-a009-cb6b5fc438a4
# ╠═6b642fe5-5f55-4a9d-86cb-1695dca29ce2
# ╠═ffdadd86-0e63-4cdc-b770-b1693614f1e2
# ╠═3eecf083-f042-4514-ad6c-5bf265aa0607
# ╠═5a13e1e1-8679-4c6c-85d6-081e5abb9e57
# ╠═66b74aa7-8cd4-4e72-a9ec-9f6671e9c547
# ╠═177f250e-3bae-4f2b-a980-ca5357fc81bf
# ╠═2664248e-dfee-4c06-aa73-55f75c9f122d
# ╠═6abd282d-aea2-4b9f-980c-9f9621f80442
# ╠═3ed32697-10d4-4b76-8540-d4431e5a3fc0
# ╠═a40752df-cfa7-40fa-84b0-8e380e69a924
# ╠═0a8baaba-06a6-4d44-bbc7-437e740be7b2
# ╠═b0d9d17b-d421-42ea-91ae-cf530675a44f
# ╠═0db1ec72-9f82-4e6c-928e-f7ef7ced0966
# ╠═ea5d4c44-c9cb-4603-8ad5-0d6070779e5d
# ╠═7a68a368-5493-4090-ad6a-b9a7c3b10cb5
# ╠═a760e4c4-e067-4385-b449-c2f897042839
# ╠═2b6e3ece-6a39-496b-ae78-a2fa6d2edd50
# ╠═c1ed3b0f-1c92-4a11-b424-042551275f8f
# ╠═26c4cfe9-1596-414a-aa31-a16158d8910b
# ╠═3e9fd7f8-56c0-4b3d-a422-7a8faa7a8245
# ╠═c88f5014-7abb-4b50-8a52-2533cda749f5
# ╠═0dab9ef6-0cdb-4983-b87a-bf5d41ae9661
# ╠═16dbb6bc-5d8f-43a7-8063-caaba9c4e71b
# ╠═84023e48-6534-42a8-a4f4-23d242d2ca88
# ╠═aa40b9e7-0992-4228-8e93-a7ad91e7f8c6
# ╠═907dccfa-0726-487f-91af-ca6f1c096d84
# ╠═ed77f48f-0df0-4718-8862-673c01250dfe
# ╠═28a90470-6d12-4de7-9f43-f53681ca7427
# ╠═fe28cfc3-d0b2-4584-8abc-41b7bea8b1c1
# ╠═e97a39c5-6d3a-4169-8143-4f02dece1e04
# ╠═9980987c-ac9d-4786-b31e-54f9ee45d467
# ╠═75e778de-9c70-43bf-b5c4-a0d2c6b1a539
# ╠═154277d2-cdfa-4f7f-8f07-63ddd4e90ae8
# ╠═820beb62-fe16-494b-8e76-323ce30e46d0
# ╠═3b85b24d-34bb-4df7-a0c3-cb89f7d90b84
# ╠═34634f51-0ec5-4ef1-8b21-b18ee1c6ced0
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═cb7f11ec-e468-41f3-b9fc-c6f12a27d5b3
# ╠═facc11b5-f86a-4390-ac68-a3a74aa8d247
# ╠═69fa05aa-8208-414c-9b92-a4a1f7be1c68
# ╠═cb56c6c6-3d92-43af-8556-232f7387dc97
# ╠═7898bcc1-651d-4a61-aff5-d8f058bf74b4
# ╠═d7f05829-b372-4ad0-afc4-645f1a37193f
# ╠═df8347f2-be16-413f-9a93-e3baa9ed8639
# ╠═f9dd1f85-da7f-41db-af12-d1b26ca36f99
# ╠═d4a11884-b670-4ceb-9444-36d3bcc6e8a7
# ╠═81181098-822a-48f9-b776-680735de6430
# ╠═4b7e4552-ffb9-4f9a-a25f-f45178b2fbd3
# ╠═f0aab374-2e76-47e0-a62c-a5ec21129767
# ╠═1b968f57-798f-471a-9203-437f373cb65f
# ╠═7f3d36a1-ac31-4cbb-b994-cbaf946d25fd
# ╠═a3c2e539-eda7-4784-b90c-43fcf6982347
# ╠═9b37fd0a-99eb-492a-a271-f584b043ef89
# ╠═83df15d6-5493-48ef-a8f4-29b8fc375da3
# ╠═418f05eb-9b20-418f-9a35-783eac94e501
