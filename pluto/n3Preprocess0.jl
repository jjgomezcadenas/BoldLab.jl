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

# ╔═╡ d68f1939-4a56-4848-884a-18fdc75bc0af
"""
    subtract_background_from_stack(img_stack::Array{<:Real,3}, background::Array{<:Real,2}) -> Array{Float64,3}

Subtract a 2D background image from each frame of a 3D image stack, clamping negative values to zero.

# Arguments
- `img_stack::Array{<:Real,3}`: A 3D array of shape (height, width, time), typically a time series of images.
- `background::Array{<:Real,2}`: A 2D image (height × width) representing the estimated background to subtract from each frame.

# Returns
- `Array{Float64,3}`: A new 3D array of the same shape as `img_stack`, where the background has been subtracted frame-wise and values are clamped to ≥ 0.


"""
subtract_background_from_stack(img2::Array{<:Real,3}, background::Array{<:Real,2}) = max.(img2 .- reshape(background, size(background)..., 1), 0.0)


# ╔═╡ d17b6d05-266d-4e2f-b38d-26e6566988ce
function denoise_stack(img2::Array{<:Real,3}; σ::Real=1.0)
	k = Kernel.gaussian(σ)
	# Apply Gaussian denoising to each frame
	img_denoised_stack = imfilter.(eachslice(img2, dims=3), Ref(k))
	imgd = cat(img_denoised_stack...; dims=3)
end


# ╔═╡ 0ed71506-c320-4373-b844-a569cff9804f


# ╔═╡ f9a54f94-b109-4ae0-abe1-f44fdedc0d25
begin
	threshold = 30.0  
	min_area = 10
end

# ╔═╡ e90dbe16-7f7f-4fa2-b0a5-9cf10579f319
"""
    filter_regions_in_stack(binary_stack::Vector{BitMatrix}; min_area::Int=5)

Process a stack of binary masks (one per frame), labeling connected regions
in each frame, filtering out small regions, and computing region properties.

# Arguments
- `binary_stack`: A vector of `BitMatrix`, one per time frame. Each mask contains `true` where a molecule is detected.
- `min_area`: Minimum number of pixels a region must have to be kept (default = 5).

# Returns
- `filtered_stack`: A vector of filtered `BitMatrix`, same length as `binary_stack`.
- `all_props`: A vector (per frame) of lists of region property `Dict`s. Each region's dictionary includes:
    - `:label` (region ID),
    - `:area` (in pixels),
    - `:centroid` (x, y coordinates).


"""
function filter_regions_in_stack(img2::Array{<:Real,3}; i_thr::Real = 30.0, min_area::Int=5)
	
	binary_stack = map(frame -> frame .> i_thr, eachslice(img2, dims=3))
	#ibst = cat(binary_stack...; dims=3)
    
    filtered_stack = BitMatrix[]  # output cleaned masks
	all_props = Vector{Vector{Dict{Symbol, Any}}}()  # per-frame region properties
	
	for bin_frame in binary_stack
	    # Label connected regions in the binary frame
	    labels = label_components(bin_frame, strel_diamond((3, 3)))
		
		#label = label_components(A, se; [bkg])
	    region_labels = setdiff(unique(labels), 0)  # ignore background (label 0)
	
	    # Initialize output mask and list of region properties
	    mask_filtered = falses(size(bin_frame))
	    props = Dict{Symbol, Any}[]
	
	    # Loop through each labeled region
	    for label in region_labels
	        inds = findall(==(label), labels)
	        area = length(inds)
	
	        # Keep only if area is above threshold
	        if area ≥ min_area
	            mask_filtered[inds] .= true  # mark region in filtered mask
	
	            # Compute centroid (mean x, y)
	            ys = [I[1] for I in inds]
	            xs = [I[2] for I in inds]
	            centroid = (mean(xs), mean(ys))
	
	            # Store region properties
	            push!(props, Dict(
	                :label => label,
	                :area => area,
	                :centroid => centroid,
					:coords => inds 
	            ))
	        end
	    end
	
	    push!(filtered_stack, mask_filtered)
	    push!(all_props, props)
	end
	
	return filtered_stack, all_props
end

# ╔═╡ a22becf2-db98-4abe-8f5a-c1dff911b4d8
begin
	
end

# ╔═╡ 8a9916b3-7f0c-43ba-a3de-0616b88060d9


# ╔═╡ 958edf9b-7bab-4da1-83ad-c85e8a04677a
function unique_i_j_across_t(df::DataFrame, t_range::UnitRange{Int})
    seen_coords = Set{Tuple{Int, Int}}()
    result_rows = DataFrame()

    for t_val in t_range
        subdf = filter(row -> row.t == t_val, df)
        keep_rows = similar(subdf, 0)  # empty sub-DataFrame with same schema

        for row in eachrow(subdf)
            coord = (row.i, row.j)
            if !(coord in seen_coords)
                push!(seen_coords, coord)
                push!(keep_rows, row)
            end
        end

        append!(result_rows, keep_rows)
    end

    return result_rows
end

# ╔═╡ c58aaada-7974-4ba6-b832-547dc3e383ce
"""
    build_sparse_traces(imst, df)

Construct a sparse matrix of traces for (i, j) pixel positions listed in a DataFrame.
Each entry contains a `Vector{T}` with the time trace at that pixel.

# Arguments
- `imst::Array{T,3}`: Image stack of shape (height, width, time)
- `df::DataFrame`: Must contain columns `i` and `j` (pixel coordinates)

# Returns
- `SparseMatrixCSC{Vector{T}}`: (height × width) matrix where only selected (i,j) entries are filled with traces
"""
function build_sparse_traces(imst::Array{T,3}, df::DataFrame) where {T<:Real}
    n, m, t = size(imst)
    
    # Sparse matrix of vector traces
    traces = spzeros(Vector{T}, n, m)
    
    for row in eachrow(df)
        i, j = row.i, row.j
        if 1 ≤ i ≤ n && 1 ≤ j ≤ m
            traces[i, j] = vec(imst[i, j, :])
        end
    end
    
    return traces
end

# ╔═╡ d8d064c2-9062-4a42-ac43-d78f1ff2ec13

"""
    detections_dataframe(region_stats, img_stack; px_size_nm, img_unit="photons")

Convert region statistics into a structured DataFrame including physical units and pixel indices.

# Arguments
- `region_stats`: Output from `filter_regions_in_stack`, containing per-frame region info.
- `img_stack`: 3D array (height, width, time), e.g. denoised fluorescence image stack.
- `px_size_nm`: Pixel size in nanometers (e.g. 130.0).
- `img_unit`: Optional string to label the intensity column (e.g. "photons").

# Returns
- A `DataFrame` with the following columns:
    - `frame`: Frame index (1-based)
    - `x [nm]`, `y [nm]`: Centroid in nanometers
    - `i`, `j`: Pixel indices (row, column) nearest to centroid
    - `area [nm²]`: Region area in nm²
    - `intensity`: Integrated intensity per region
"""
function detections_dataframe(
    region_stats::Vector{Vector{Dict{Symbol, Any}}},
    img_stack::Array{<:Real,3};
    px_size_nm::Real
)
    detections = []

    for (t, props) in enumerate(region_stats)
        img = img_stack[:, :, t]

        for r in props
            coords = r[:centroid]
            area_px = r[:area]
            inds = get(r, :coords, [])

            # Safe pixel intensity sum
            intensity = isempty(inds) ? 0.0 : sum(img[i] for i in inds)

            # Nearest integer indices for centroid
            i = round(Int, coords[2])  # row (y)
            j = round(Int, coords[1])  # column (x)

            push!(detections, (
                frame = t,
                x = coords[1] * px_size_nm,
                y = coords[2] * px_size_nm,
                i = i,
                j = j,
				t = t,
                area = area_px * px_size_nm^2,
                intensity = intensity
            ))
        end
    end

    df = DataFrame(detections)

    rename!(df, Dict(
    :x => "x_nm",
    :y => "y_nm",
    :i => "i",
    :j => "j",
	:t => "t",
    :area => "area_nm2",
    :intensity => "intensity"
))

    return df
end

# ╔═╡ b9f8446c-76e7-42d7-b8b6-6ec2db8b3857
function plot_frames_with_centroids(imst::AbstractArray{T,3}, region_stats::Vector{Vector{Dict{Symbol, Any}}}; 
	nscale::Int=20) where {T<:Real}

	FF = []

	
    for i in 1:9
        fn = (i-1) * nscale + i
		if fn > size(imst)[3]
			warn("requested frame = $(fn) is to large, set smaller nscale")
			fn = size(imst)[3] -i
			warn("set fn = $(fn)")
		end
		frame = imst[:, :, fn]
		centroids = [r[:centroid] for r in region_stats[fn]]
	
		p = heatmap(frame, 
					color=:grays, 
					colorbar=false,
					title="Frame $fn",
				 	titlefontsize=7,
					tickfontsize=6,
					guidefontsize=6,
					titlelocation=:left,
					aspect_ratio=:equal)
		p = scatter!(p, [c[1] for c in centroids], [c[2] for c in centroids], markersize=4, color=:red)
        push!(FF, p)
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

# ╔═╡ b332855d-7caa-483e-8d5a-d3121b1ff877


# ╔═╡ abe91220-5762-4932-b03d-925a35cd4588
function find_fit_candidates2(trzs, df; sel="core", ped=0.0, niter=5, thr=0.5)
	I = Int[]
	J = Int[]
    DX = Vector{Vector{Float32}}()
    FX = Vector{Vector{Float32}}()
	SC = Vector{Vector{Float32}}()
	
	ng = 0

	df2 = DataFrame(i=Int[], j=Int[], nmol=Int[], nstep=Int[],
                   stepHeight=Float32[], stepTime=Int[], stepLength=Int[])

	for row in eachrow(df)
		trace = trzs[row.i, row.j] 
    	dataX, FitX, S_curve, best_shot, iter, cc = fit_traces(trace, niter=niter, 														tresH=thr, sel=sel)
		if best_shot >0
			ng+=1
			push!(I,row.i)
			push!(J,row.j)
			push!(DX, dataX .-ped)
			push!(FX, FitX .- ped)
			push!(SC, S_curve)

			sth, stt, stl = getsteps(FitX)
			nsteps = length(sth)

			for k in 1:nsteps
                push!(df2, (row.i, row.j, ng, nsteps, sth[k] - ped, stt[k], stl[k]))
			end
            
		end
	end
	
	df2, I, J, DX, FX, SC
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

# ╔═╡ 08dd1429-4c2d-4def-8282-e7abe469f318
md"""
### Fit data
"""

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

# ╔═╡ 97ba452a-90f1-47d2-88d5-da40867f012b
"""
    compute_background_from_stack(img_stack::Array{<:Real,3}; σ::Real = 10.0, nlf::Int = 5)

Estimate the per-pixel background of a fluorescence image stack by averaging the last `nlf` frames
and applying a Gaussian blur with standard deviation `σ`.

# Arguments
- `img_stack::Array{<:Real,3}`: A 3D array of shape (height, width, time) representing the image sequence.
- `σ::Real = 10.0`: Standard deviation (in pixels) of the Gaussian filter used to smooth the background.
- `nlf::Int = 5`: Number of last frames to average for estimating the background.

# Returns
- `background::Matrix{Float64}`: A 2D matrix (same spatial size as a single frame) containing the smoothed background.

# Example
```julia
background = compute_background_from_stack(img_stack; σ=8.0, nlf=10)
img_subtracted = max.(img_stack .- reshape(background, size(background)..., 1), 0.0)
"""
function compute_background_from_stack(img2::Array{<:Real,3}; σ::Real = 10.0, nlf::Int=5)
	last_frames = img2[:,:, end-nlf+1:end]
	background3d = mean(last_frames, dims=3)
	avg_bg = dropdims(background3d, dims=3)
	background = imfilter(avg_bg, Kernel.gaussian(σ))
	
end

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
	vi2m = vec(imgd[:,:,20])
	histogram(vi2m[vi2m.>30])
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

# ╔═╡ 140b9b8f-74a7-46f2-9464-c8bfe4b66151
typeof(region_stats)

# ╔═╡ 5d1f6f1d-782d-4f8b-a232-5d4d728e0097
 plot_frames_with_centroids(imgd,region_stats; nscale=20)

# ╔═╡ 01fa4813-9bbc-4460-a07b-3a01c1bcf3e3
df = detections_dataframe(region_stats, imgd; px_size_nm=260.0)

# ╔═╡ 55629773-caef-4ca0-8023-ca8ffdf1f0eb
histogram(df.intensity)

# ╔═╡ a874cd51-8946-42fe-8d07-defc41e92a1d
histogram(df.t)

# ╔═╡ 816140c5-9485-4060-ac72-3f7fe9c33063
dft1 = df[df.t.==3,:]

# ╔═╡ fe85b459-ae6c-4f9e-b97b-c8034d84423c
df5f = unique_i_j_across_t(df, 1:5)

# ╔═╡ bd6434e3-d4f0-4a73-8f06-141a3aa73f88
histogram(df5f.intensity)

# ╔═╡ 1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
TRZS = build_sparse_traces(imxt, df5f)

# ╔═╡ d92c4660-4950-4adf-aceb-bc94663411c6
function select_trace(df::DataFrame, row::Int) 
	i = df.i[row]
	j = df.j[row]
    trace = TRZS[i, j]  # safe, only defined where needed
	i,j, trace
end

# ╔═╡ 6f3a242a-6dff-432a-b30e-1b7ee1d234fc
i,j, tz = select_trace(df5f, 1)

# ╔═╡ d1c82dbe-8493-464d-bdac-f505657678d6
begin
	dataX, FitX, S_curve, best_shot, iter, cc = fit_traces(tz, niter=5, tresH=0.0, sel="core")
	sth, stt, stl = getsteps(FitX)

	

	md"""
	- Fit results
	- bes shot =$(best_shot)
	- Step heights: = $(vect_to_fstr(sth, "%.2f"))
	- Step times =$(vect_to_fstr(stt, "%.2f"))
	- Segment lengths = $(vect_to_fstr(stl, "%.2f"))
	"""
	
end

# ╔═╡ 1b250af1-ffe5-488c-9eee-c1ca41085901
typeof(S_curve)

# ╔═╡ 8aa50fbb-74a9-424d-9a84-e717247c65a9
plotfit(dataX, FitX, S_curve, i,j )

# ╔═╡ 5042becb-29c6-447e-84ad-a965a9961992
dfs, zI, zJ, zDX, zFX, zSC  =find_fit_candidates2(TRZS, df5f;  sel="core", ped=0.0, niter=3, thr=0.25)

# ╔═╡ 1015d8f6-d32c-4920-b688-ab1848cd9c5a
plot_fits(zDX, zFX, zI, zJ; plotsel="4x4")

# ╔═╡ 64ed321d-3050-465e-8c4f-fea826265779
plot_sc(zSC, zI,zJ; plotsel="4x4", figsize=(1500,1500))

# ╔═╡ 01391609-5034-4dbd-8319-10ac82126cfc
length(zDX)

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

# ╔═╡ 34634f51-0ec5-4ef1-8b21-b18ee1c6ced0
if casedata
	ss = save_stats(sdf3, npixels, folder_scan, folder_field)
	Markdown.parse(ss)
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
# ╠═d68f1939-4a56-4848-884a-18fdc75bc0af
# ╠═b0a4c0be-12b9-416c-a45c-5fe35fbbd168
# ╠═89354633-e220-45f2-989e-c5575acd2988
# ╠═f4ff7938-f16e-4159-95e8-cb98c59a9d80
# ╠═d17b6d05-266d-4e2f-b38d-26e6566988ce
# ╠═0ed71506-c320-4373-b844-a569cff9804f
# ╠═fd39234f-6c62-42e3-bde9-9a9e565fa519
# ╠═f9a54f94-b109-4ae0-abe1-f44fdedc0d25
# ╠═e90dbe16-7f7f-4fa2-b0a5-9cf10579f319
# ╠═389be056-e3d5-45fc-b80c-19f7cd4d5433
# ╠═a22becf2-db98-4abe-8f5a-c1dff911b4d8
# ╠═5d1f6f1d-782d-4f8b-a232-5d4d728e0097
# ╠═8a9916b3-7f0c-43ba-a3de-0616b88060d9
# ╠═01fa4813-9bbc-4460-a07b-3a01c1bcf3e3
# ╠═55629773-caef-4ca0-8023-ca8ffdf1f0eb
# ╠═a874cd51-8946-42fe-8d07-defc41e92a1d
# ╠═816140c5-9485-4060-ac72-3f7fe9c33063
# ╠═fe85b459-ae6c-4f9e-b97b-c8034d84423c
# ╠═bd6434e3-d4f0-4a73-8f06-141a3aa73f88
# ╠═958edf9b-7bab-4da1-83ad-c85e8a04677a
# ╠═c58aaada-7974-4ba6-b832-547dc3e383ce
# ╠═d8d064c2-9062-4a42-ac43-d78f1ff2ec13
# ╠═b9f8446c-76e7-42d7-b8b6-6ec2db8b3857
# ╠═140b9b8f-74a7-46f2-9464-c8bfe4b66151
# ╠═783bbb4f-1bb1-4d57-881d-6a6f3c61e25b
# ╠═57692746-065d-4d75-8a41-7ffafd69550e
# ╠═74deedfc-a344-4373-9fc9-db22e83d48ac
# ╠═053dd59b-f6cb-47b2-9bfe-ae0c7f2a514b
# ╠═57bd56be-5729-44e8-aba3-783a78c714d2
# ╠═6b54a654-7850-4a27-8ec4-2dcbd3566d82
# ╠═7b89b967-0b2e-4f8d-8b8f-c6597a20a638
# ╠═1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
# ╠═d92c4660-4950-4adf-aceb-bc94663411c6
# ╠═6f3a242a-6dff-432a-b30e-1b7ee1d234fc
# ╠═3a8e7300-d209-4e50-ace8-c4b841d71f42
# ╠═b332855d-7caa-483e-8d5a-d3121b1ff877
# ╠═d1c82dbe-8493-464d-bdac-f505657678d6
# ╠═8aa50fbb-74a9-424d-9a84-e717247c65a9
# ╠═1b250af1-ffe5-488c-9eee-c1ca41085901
# ╠═abe91220-5762-4932-b03d-925a35cd4588
# ╠═fa0da6ef-68ba-4db2-aec0-26fd3a280b5a
# ╠═efd033eb-bde9-4f4e-85f0-777125555edd
# ╠═d9477e2f-2bad-4243-a11f-393f5b3405c7
# ╠═862d9a0c-ecc6-4178-80f4-74fc71a79f38
# ╠═5042becb-29c6-447e-84ad-a965a9961992
# ╠═1015d8f6-d32c-4920-b688-ab1848cd9c5a
# ╠═64ed321d-3050-465e-8c4f-fea826265779
# ╠═01391609-5034-4dbd-8319-10ac82126cfc
# ╠═0a8baaba-06a6-4d44-bbc7-437e740be7b2
# ╠═77f5c1c3-11bb-4976-af44-1b488f08de6b
# ╠═3115a147-d41b-4ab6-9ad9-0f3b30fb4504
# ╠═e04c9c53-7ae2-474b-87bc-97bd58f458fa
# ╠═f175508e-f0ea-48a9-ba72-d2e21215de1d
# ╠═86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
# ╠═264fcbb0-ec40-42b3-ab05-01b2171da6f2
# ╠═08dd1429-4c2d-4def-8282-e7abe469f318
# ╠═16dbb6bc-5d8f-43a7-8063-caaba9c4e71b
# ╠═aa40b9e7-0992-4228-8e93-a7ad91e7f8c6
# ╠═154277d2-cdfa-4f7f-8f07-63ddd4e90ae8
# ╠═34634f51-0ec5-4ef1-8b21-b18ee1c6ced0
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═97ba452a-90f1-47d2-88d5-da40867f012b
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
