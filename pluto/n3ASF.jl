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

# ╔═╡ 7929e575-bddf-4da5-ba21-6fdcc71f81da
using Changepoints

# ╔═╡ 79722d05-bcd3-4655-9459-ab697480d652
using ChangePointDetection

# ╔═╡ 6aa0d8ec-2dc3-4b2b-9d5b-dbaaba4b6449
using Distributions 

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
set_log_level(SimpleLogger.DEBUG)

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

# ╔═╡ 62da10f6-0d5c-41c0-a985-d15c946f5b84
begin
    @bind show_peaks CheckBox(true)
end

# ╔═╡ 12960d51-135c-47cd-ab86-f2ab5bacef08
@bind nframe NumberField(0:199, default=1)

# ╔═╡ d34b9105-09d4-4876-842d-bcf74249cca9
@bind pthr NumberField(0:0.1:200.0, default=30.0)

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

# ╔═╡ 1c6507e3-821f-4e07-b41e-b3c106df3671
@bind ntrz NumberField(0:100, default=1)

# ╔═╡ e04c9c53-7ae2-474b-87bc-97bd58f458fa
@bind thr NumberField(0.0:0.1:1.1, default=0.5)

# ╔═╡ 86e3c52c-9119-4d7d-bc7a-dc2f6cd63303
@bind niter NumberField(0:10, default=3)

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

# ╔═╡ 76f1b23b-b6c4-43a0-9e53-3b2ff553f167
md"""
## Changepoints
"""

# ╔═╡ 9f5ca7a0-e06d-4c23-9b13-253b998ed371
names(ChangePointDetection)

# ╔═╡ 2b81eca2-bc1d-4d08-9813-7032be53f6cd
begin
	n = 1000
	  λ = 100
	  μ, σp = Normal(0.0, 10.0), 1.0
	  # Samples changepoints from Normal distribution with changing mean
	  sample, cps = @changepoint_sampler n λ Normal(μ, σp)
	  # Run PELT on sample
	  pelt_cps, pelt_cost = @PELT sample Normal(:?, σp)
end

# ╔═╡ 1612f3e5-7fb2-4aac-81e1-a920dd796bb7
pelt_cps

# ╔═╡ 493d9175-bf44-40a2-94a5-2e05ba9d931e
sample

# ╔═╡ 64ea787c-b719-4dd1-b0f2-1dcb278922bf
cps

# ╔═╡ e4024366-61de-474e-891d-084cd5882df0
changepoint_plot(sample, pelt_cps)

# ╔═╡ c03e418c-c205-4e60-87d6-83c4c3fddb57
#detected_steps = kv_step_detection(tz; min_step_size=150.0)

# ╔═╡ a7b044e2-b219-449e-9759-f308641d719e


# ╔═╡ dece3304-7765-4c2f-89f3-143d14bf6e18
#fit_result = sic_fit_driver(tz, collect(1:length(tz)))

# ╔═╡ b86a9318-79d0-4169-8744-665cb6d54337
md"""
## KV code
"""

# ╔═╡ 4e451c59-c7de-4a55-9fcd-4d55cb06170b
function make_time_vector(tz::AbstractVector{T}, dt::Real) where {T<:Real}
    n = length(tz)
    return dt .* collect(0:n-1)  # time starts at 0, dt spacing
end

# ╔═╡ 14caa5f7-b3f6-4e70-963d-d8965640dbed
begin

# Struct to hold the step fit result
mutable struct StepFit{T <: Real}
    means::Vector{T}
    step_locations::Vector{Int}
    SIC::T
    chisq::T
end

"Initialize StepFit fields with zeros."
function initialize_fit_to_zeros(T::Type{<:Real}, n_dwells::Int)
    StepFit{T}(zeros(T, n_dwells), zeros(Int, n_dwells - 1), zero(T), zero(T))
end

"Set all values in an existing StepFit to zero (in-place reset)."
function set_fit_to_zeros!(fit::StepFit{T}) where {T}
    fill!(fit.means, zero(T))
    fill!(fit.step_locations, 0)
    fit.SIC = zero(T)
    fit.chisq = zero(T)
    return fit
end

"Get the mean of the dwell that a given index belongs to."
function get_previous_dwell_mean(prev_fit::StepFit{T}, step_indices::Vector{Int}, index::Int) where {T}
    loc = sum(step_indices[1:index-1])
    return prev_fit.means[loc + 1]
end

"Compute the no-step SIC model (mean of full trace)."
function no_step_sic(pos::AbstractVector{T}) where {T <: Real}
    n = length(pos)
    fit = initialize_fit_to_zeros(T, 1)
    fit.means[1] = mean(pos)
    fit.chisq = mean((pos .- fit.means[1]).^2)
    fit.SIC = T(2) * log(n) + T(n) * log(fit.chisq)
    return fit
end

"Try adding a step at every valid position and return the best resulting fit."
function add_step_sic2(pos::AbstractVector{T}, step_indices::Vector{Int},
                      n_dwell::Int, prev_fit::StepFit{T}) where {T <: Real}
    n = length(pos)
    win_fit = initialize_fit_to_zeros(T, n_dwell + 1)
    win_fit.SIC = T(-1)
    win_step_location = 0

    for i in 2:(n - 1)
        if step_indices[i] == 0
            step_indices[i] = 1  # trial step

			# get means of the dwells
            # compute left and right means
            left_vals = T[]
            j = i - 1
            while j ≥ 1 && step_indices[j] == 0
                push!(left_vals, pos[j])
                j -= 1
            end
            tmpleft = mean(left_vals)

			debug("i = $(i), left_vals =$(left_vals)")
			debug("tmpleft = $(tmpleft)")

            right_vals = T[]
            j = i + 1
            while j ≤ n && step_indices[j] == 0
                push!(right_vals, pos[j])
                j += 1
            end
            tmpright = mean(right_vals)

			debug("i = $(i), right_vals =$(right_vals)")
			debug("tmpright = $(tmpright)")

            prev_mean = get_previous_dwell_mean(prev_fit, step_indices, i)

            # compute updated chisq
            chisq = prev_fit.chisq * T(n)

			debug("prev_mean = $(prev_mean), chisq =$(chisq)")

			
            j = i - 1
            while j ≥ 1 && step_indices[j] == 0
                chisq -= (pos[j] - prev_mean)^2
                chisq += (pos[j] - tmpleft)^2
                j -= 1
            end

            j = i
            while j ≤ n && step_indices[j] == 0
                chisq -= (pos[j] - prev_mean)^2
                chisq += (pos[j] - tmpright)^2
                j += 1
            end

            chisq /= T(n)
            curr_SIC = T(n_dwell + 2) * log(n) + T(n) * log(chisq)

            if curr_SIC < win_fit.SIC || win_fit.SIC == T(-1)
                win_fit = initialize_fit_to_zeros(T, n_dwell + 1)

                # build updated means vector
                curr_means = T[]
                for k in 1:n_dwell
                    m = prev_fit.means[k]
                    if m != prev_mean
                        push!(curr_means, m)
                    else
                        push!(curr_means, tmpleft)
                        push!(curr_means, tmpright)
                    end
                end

                win_fit.means .= curr_means
                win_fit.chisq = chisq
                win_fit.SIC = curr_SIC
                win_step_location = i
            end

            step_indices[i] = 0  # reset trial
        elseif step_indices[i] != 1
            error("step_indices must contain only 0 or 1")
        end
    end

    final_steps = findall(==(1), step_indices)
    if win_step_location > 0
        push!(final_steps, win_step_location)
    end
    win_fit.step_locations .= sort!(final_steps)

    return win_fit
end


end

# ╔═╡ 16cbc179-3710-43f6-b2ed-7eb3cbb9109f
function add_step_sic(pos::Vector{T}, step_indices::Vector{Int},
                      n_dwell::Int, prev_fit::StepFit{T}) where {T <: Real}
    n = length(pos)

    curr_fit = initialize_fit_to_zeros(T, n_dwell + 1)
    win_fit = initialize_fit_to_zeros(T, n_dwell + 1)
    win_fit.SIC = T(-1)
    win_step_location = 0

    for i in 2:(n - 1)
        if step_indices[i] == 0
            step_indices[i] = 1  # trial step location

            # --- compute means to the left and right of step ---

			debug("i = $(i)")
            
			j = i - 1

			debug("left: j = $(j)")
			
            tmpleft_sum = zero(T)
            count_left = 0
            while j ≥ 1 && step_indices[j] == 0
                tmpleft_sum += pos[j]
                count_left += 1
                j -= 1
            end
            tmpleft = tmpleft_sum / count_left

			debug("count_left = $(count_left), tmpleft = $(tmpleft)")

            j = i + 1

			debug("right: j = $(j)")
			
            tmpright_sum = zero(T)
            count_right = 0
            while j ≤ n && step_indices[j] == 0
                tmpright_sum += pos[j]
                count_right += 1
                j += 1
            end
            tmpright = tmpright_sum / count_right

			debug("count_right = $(count_right), tmpright = $(tmpright)")

            # --- compute updated chisq ---
            prev_mean = get_previous_dwell_mean(prev_fit, step_indices, i)
            curr_fit.chisq = prev_fit.chisq * T(n)

			debug("prev_mean = $(prev_mean)")

            j = i - 1
            while j ≥ 1 && step_indices[j] == 0
                curr_fit.chisq -= (pos[j] - prev_mean)^2
                curr_fit.chisq += (pos[j] - tmpleft)^2
                j -= 1
            end

			debug("left: curr_fit.chisq = $(curr_fit.chisq)")

            j = i
            while j ≤ n && step_indices[j] == 0
                curr_fit.chisq -= (pos[j] - prev_mean)^2
                curr_fit.chisq += (pos[j] - tmpright)^2
                j += 1
            end

			debug("right: curr_fit.chisq = $(curr_fit.chisq)")

            curr_fit.chisq /= T(n)
            curr_fit.SIC = T(n_dwell + 2) * log(n) + T(n) * log(curr_fit.chisq)

			debug("curr_fit.chisq = $(curr_fit.chisq)")
			debug("curr_fit.SIC = $(curr_fit.SIC)")
			

            if curr_fit.SIC < win_fit.SIC || win_fit.SIC == T(-1)

				debug("win_fit.SIC = $(win_fit.SIC)")
                win_fit = initialize_fit_to_zeros(T, n_dwell + 1)
                win_step_location = i

				debug("win_step_location = $(win_step_location)")

                # update curr_fit.means from prev_fit, splitting prev_mean
                fill!(curr_fit.means, zero(T))
                dum = 1
                for k in 1:n_dwell
                    m = prev_fit.means[k]
                    if m ≠ prev_mean
                        curr_fit.means[dum] = m
                        dum += 1
                    else
                        curr_fit.means[dum] = tmpleft; dum += 1
                        curr_fit.means[dum] = tmpright; dum += 1
                    end
                end

                # copy curr_fit → win_fit
                win_fit.means .= curr_fit.means
                win_fit.SIC = curr_fit.SIC
                win_fit.chisq = curr_fit.chisq

				debug("copy curr_fit → win_fit")
				debug("win_fit.SIC =$(win_fit.SIC)")
            end

            # reset current trial
            curr_fit = initialize_fit_to_zeros(T, n_dwell + 1)
            step_indices[i] = 0

        elseif step_indices[i] == 1
            # do nothing
        else
            error("step_indices must only contain 0 or 1")
        end
    end

    # --- update step_locations in win_fit ---
    final_steps = Int[]
    for k in 1:n
        if step_indices[k] == 1
            push!(final_steps, k)
        end
        if k == win_step_location
            push!(final_steps, k)
        end
    end
    win_fit.step_locations .= sort!(final_steps)

    return win_fit
end

# ╔═╡ c7524cff-3c4b-4da5-a5b8-61e330881f32
"""
    sic_fit_driver(position::Vector{T}, time::Vector{T}) where {T<:Real}

Iteratively applies the SIC-based step-fitting algorithm to the input trace
until no improvement is found.

Returns a vector of (time, mean) tuples representing the fitted steps.
"""
function sic_fit_driver(position::Vector{T}, time::Vector{Int}) where {T<:Real}
    n = length(position)
    @assert length(time) == n "Time and position vectors must match"

    this_fit = no_step_sic(position)
    n_dwells = 1
    step_indices = zeros(Int, n)
	prev_fit = initialize_fit_to_zeros(T, n_dwells)

	debug("position=$(position[1:10])")
	debug("time=$(time[1:10])")
    while true
        debug("Iteration with n_dwells = $n_dwells")
        debug("Current SIC = $(this_fit.SIC)")

        # Copy this_fit to prev_fit
        prev_fit = initialize_fit_to_zeros(T, n_dwells)

        prev_fit.means .= this_fit.means
        if n_dwells > 1
            prev_fit.step_locations .= this_fit.step_locations[1:n_dwells-1]
        end
        prev_fit.SIC = this_fit.SIC
        prev_fit.chisq = this_fit.chisq

        # Reset and update step_indices based on prev_fit
        step_indices .= 0
        for s in prev_fit.step_locations
            if s > 0
                step_indices[s] = 1
            end
        end

        # Try to add a new step
        this_fit = add_step_sic(position, step_indices, n_dwells, prev_fit)

        debug("Proposed new SIC = $(this_fit.SIC)")

        # If no improvement, stop
        if this_fit.SIC ≥ prev_fit.SIC
            debug("No further SIC improvement. Halting.")
            break
        end

        n_dwells += 1
    end

    # Format result as list of (time, mean)
    result = [(time[1], prev_fit.means[1])]
    for i in 1:(n_dwells - 2)
        t_step = time[prev_fit.step_locations[i]]
        push!(result, (t_step, prev_fit.means[i]))
        push!(result, (t_step, prev_fit.means[i + 1]))
    end
    push!(result, (time[end], prev_fit.means[n_dwells - 1]))

    return result
end

# ╔═╡ 431fa5c2-8056-4e24-9776-e57deae9d24b
"""
    kv_step_detection(trace::Vector{<:Real}; min_step_size::Float64=0.0)

Detect step changes in a 1D time series using the Kalafut–Visscher algorithm
based on the Schwarz Information Criterion (SIC).

# Arguments
- `trace::Vector{<:Real}`: The input 1D signal.
- `min_step_size::Float64=0.0`: Minimum difference in segment means to accept a step.

# Returns
- A sorted `Vector{Int}` of step positions (indices where steps occur).
"""
function kv_step_detection(trace::Vector{<:Real}; min_step_size::Float64 = 0.0)
    N = length(trace)
    steps = Int[]
    segments = [(1, N)]

    while true
        best_sic = Inf
        best_split = nothing
        best_segment_idx = nothing

        for (seg_idx, (i_start, i_end)) in enumerate(segments)
            if i_end - i_start + 1 < 3
                continue
            end

            for s in (i_start + 1):(i_end - 1)
                left = trace[i_start:s]
                right = trace[s+1:i_end]

                mean_left = mean(left)
                mean_right = mean(right)
                mean_total = mean(trace[i_start:i_end])

                if abs(mean_left - mean_right) < min_step_size
                    continue
                end

                ssr_left = sum((left .- mean_left).^2)
                ssr_right = sum((right .- mean_right).^2)
                ssr_total = sum((trace[i_start:i_end] .- mean_total).^2)

                n = i_end - i_start + 1
                k0 = 1
                k1 = 2

                sic_nostep = n * log(ssr_total / n) + k0 * log(n)
                sic_step = length(left) * log(ssr_left / length(left)) +
                           length(right) * log(ssr_right / length(right)) +
                           k1 * log(n)

                if sic_step < sic_nostep && sic_step < best_sic
                    best_sic = sic_step
                    best_split = s
                    best_segment_idx = seg_idx
                end
            end
        end

        if best_split === nothing
            break
        end

        push!(steps, best_split)

        # Update segments
        i_start, i_end = segments[best_segment_idx]
        deleteat!(segments, best_segment_idx)

        insert!(segments, best_segment_idx, (best_split + 1, i_end))
        insert!(segments, best_segment_idx, (i_start, best_split))
    end

    return sort(steps)
end

# ╔═╡ 15435237-8a55-48c8-ab07-101dad94b86b
md"""
## Theory
"""

# ╔═╡ a1c1511c-d5f5-4334-b4e2-bdd86b8d941a
md"""
## `quickpbsa` and Step Detection Algorithms

[`quickpbsa`](https://github.com/jjgomezcadenas/quickpbsa) is a Julia toolkit for **photobleaching step analysis (PBSA)**.  
It detects and quantifies discrete steps in noisy fluorescence time traces — typically arising from single-molecule photobleaching.

Two key statistical tools are used:

- **Schwarz Information Criterion (SIC)** for model selection
- **Kalafut–Visscher algorithm** for step detection
"""

# ╔═╡ e6721d12-784a-4011-b5b5-92503194b6e3
md"""
### 📊 Schwarz Information Criterion (SIC / BIC)

The **Schwarz Information Criterion (SIC)**, also known as the **Bayesian Information Criterion (BIC)**, helps select the best model by balancing fit quality and model simplicity.

It is defined as:

```BIC = k * log(n) - 2 * log(L̂)```

- `k`: number of parameters (e.g., number of steps),
- `n`: number of data points in the trace,
- `L̂`: maximized likelihood of the model.

**Use in PBSA**:
SIC penalizes overly complex models to prevent overfitting. The model with the **lowest SIC/BIC value** is considered optimal.
"""

# ╔═╡ 5811fbff-30a9-4c6b-9f5f-15c35e3d3e96
md"""
### 🧬 Kalafut–Visscher (KV) Algorithm

The **Kalafut–Visscher algorithm** is a nonparametric method for detecting step-like changes in noisy time-series data (e.g., fluorescence intensity).

**Key steps**:
1. Propose candidate step locations.
2. Compare models **with** and **without** each step using the SIC.
3. Accept the step if it improves the SIC.
4. Iterate until no further steps improve the fit.

**Why it's useful**:
- No assumptions about noise or number of steps.
- Robust to noise and baseline fluctuations.
- Fully automated and reproducible.

🔗 [Original C implementation](https://github.com/knyquist/KV_SIC)
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

# ╔═╡ 1d9ae22b-6cb6-49c4-81e7-b4b740c893a7
TRZS = build_sparse_traces(imxt, peaks)

# ╔═╡ 6f3a242a-6dff-432a-b30e-1b7ee1d234fc
i,j, tz = select_trace(TRZS, peaks, ntrz)

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

# ╔═╡ 1b250af1-ffe5-488c-9eee-c1ca41085901
typeof(S_curve)

# ╔═╡ 8aa50fbb-74a9-424d-9a84-e717247c65a9
plotfit(dataX, FitX, S_curve, i,j )

# ╔═╡ 861ecdd4-de6a-4a85-9d7c-303db260f009
std(tz[1:10])

# ╔═╡ 92b7aaf4-11cd-4fea-bca0-6718e04ac105
std(tz[end-20:end])

# ╔═╡ 7a259526-2cd1-4b3d-9188-742d4be85666
costf = NormalMeanVarSegment(tz)

# ╔═╡ bd94fe1d-1aa7-40b6-9bb1-8a90a56b8b09
pcps, pcst = PELT(costf, length(tz))

# ╔═╡ 54c5dace-4d42-4ae1-b00e-25409df515c6
changepoint_plot(tz, pcps)

# ╔═╡ 59b02014-4ded-4aee-affb-20942ac84392
mean(tz[1:30])

# ╔═╡ d0e8c78a-77dd-403a-ae5a-23bc5524dbc7
std(tz[1:31])

# ╔═╡ a1b0eec7-7f29-412f-90c4-8be4c826966f
mean(tz[30:end])

# ╔═╡ 113a84f8-35d4-47ed-ab00-f7b572094b56
tzt = make_time_vector(tz, 1.0)

# ╔═╡ 02ba78f2-1de5-46d1-a8d6-d6895b2d8771
peaks

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

	hstphx, pstphx = step_hist(df_pix.stepHeightMin;
              nbins=50,
              xlim=(0.0, 500.0),
              logy=false,
              xlabel=" # step height min",
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
# ╠═02ba78f2-1de5-46d1-a8d6-d6895b2d8771
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
# ╠═76f1b23b-b6c4-43a0-9e53-3b2ff553f167
# ╠═7929e575-bddf-4da5-ba21-6fdcc71f81da
# ╠═79722d05-bcd3-4655-9459-ab697480d652
# ╠═9f5ca7a0-e06d-4c23-9b13-253b998ed371
# ╠═6aa0d8ec-2dc3-4b2b-9d5b-dbaaba4b6449
# ╠═2b81eca2-bc1d-4d08-9813-7032be53f6cd
# ╠═1612f3e5-7fb2-4aac-81e1-a920dd796bb7
# ╠═493d9175-bf44-40a2-94a5-2e05ba9d931e
# ╠═64ea787c-b719-4dd1-b0f2-1dcb278922bf
# ╠═e4024366-61de-474e-891d-084cd5882df0
# ╠═c03e418c-c205-4e60-87d6-83c4c3fddb57
# ╠═861ecdd4-de6a-4a85-9d7c-303db260f009
# ╠═92b7aaf4-11cd-4fea-bca0-6718e04ac105
# ╠═7a259526-2cd1-4b3d-9188-742d4be85666
# ╠═bd94fe1d-1aa7-40b6-9bb1-8a90a56b8b09
# ╠═54c5dace-4d42-4ae1-b00e-25409df515c6
# ╠═59b02014-4ded-4aee-affb-20942ac84392
# ╠═d0e8c78a-77dd-403a-ae5a-23bc5524dbc7
# ╠═a7b044e2-b219-449e-9759-f308641d719e
# ╠═a1b0eec7-7f29-412f-90c4-8be4c826966f
# ╠═113a84f8-35d4-47ed-ab00-f7b572094b56
# ╠═dece3304-7765-4c2f-89f3-143d14bf6e18
# ╠═b86a9318-79d0-4169-8744-665cb6d54337
# ╠═4e451c59-c7de-4a55-9fcd-4d55cb06170b
# ╠═16cbc179-3710-43f6-b2ed-7eb3cbb9109f
# ╠═c7524cff-3c4b-4da5-a5b8-61e330881f32
# ╠═14caa5f7-b3f6-4e70-963d-d8965640dbed
# ╠═431fa5c2-8056-4e24-9776-e57deae9d24b
# ╠═15435237-8a55-48c8-ab07-101dad94b86b
# ╠═a1c1511c-d5f5-4334-b4e2-bdd86b8d941a
# ╠═e6721d12-784a-4011-b5b5-92503194b6e3
# ╠═5811fbff-30a9-4c6b-9f5f-15c35e3d3e96
# ╠═1be93848-5757-48e6-8d5b-638cb11c4a61
# ╠═e7cb1f63-130c-4e75-af5d-779fc1a4c755
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
