module StepPreprocessing 
using Revise
using StatsPlots
using Statistics
using StatsBase
import Measures
using SparseArrays
using DataFrames
using Images, FileIO, ImageIO, ImageView, ImageSegmentation, ImageMorphology

using NPZ

export subtract_background_from_stack, compute_background_from_stack, denoise_stack
export filter_regions_in_stack, plot_frames_with_centroids, detections_dataframe
export build_sparse_traces, unique_i_j_across_t


# plt_ft_trace


include("histos.jl")
using .histos

include("SimpleLogger.jl")
using .SimpleLogger
global_log_level = ERROR



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

"""
function compute_background_from_stack(img2::Array{<:Real,3}; σ::Real = 10.0, nlf::Int=5)
	last_frames = img2[:,:, end-nlf+1:end]
	background3d = mean(last_frames, dims=3)
	avg_bg = dropdims(background3d, dims=3)
	background = imfilter(avg_bg, Kernel.gaussian(σ))
	
end


"""
    subtract_background_from_stack(img_stack::Array{<:Real,3}, background::Array{<:Real,2}) -> Array{Float64,3}

Subtract a 2D background image from each frame of a 3D image stack, clamping negative values to zero.

# Arguments
- `img_stack::Array{<:Real,3}`: A 3D array of shape (height, width, time), typically a time series of images.
- `background::Array{<:Real,2}`: A 2D image (height × width) representing the estimated background to subtract from each frame.

# Returns
- `Array{Float64,3}`: A new 3D array of the same shape as `img_stack`, where the background has been subtracted frame-wise and values are clamped to ≥ 0.


"""
subtract_background_from_stack(img2::Array{<:Real,3}, 
background::Array{<:Real,2}) = max.(img2 .- reshape(background, size(background)..., 1), 0.0)


"""
    denoise_stack(img2::Array{<:Real,3}; σ::Real=1.0)

Apply Gaussian denoising to each frame of a 3D image stack independently.

# Arguments
- `img2::Array{<:Real,3}`: A 3D array of shape (height, width, frames), representing an image time stack.
- `σ::Real=1.0`: Standard deviation of the Gaussian kernel used for denoising.

# Returns
- A 3D array of the same shape, where each frame has been filtered with a Gaussian kernel.

# Notes
- Filtering is applied independently to each 2D frame using `imfilter`.
- The kernel is constructed once and passed by reference for efficiency.

# Example

"""
function denoise_stack(img2::Array{<:Real,3}; σ::Real=1.0)
	k = Kernel.gaussian(σ)
	# Apply Gaussian denoising to each frame
	img_denoised_stack = imfilter.(eachslice(img2, dims=3), Ref(k))
	imgd = cat(img_denoised_stack...; dims=3)
end


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


"""
    unique_i_j_across_t(df::DataFrame, t_range::UnitRange{Int}) -> DataFrame

Select rows from a DataFrame by progressively scanning frames `t` in the given range
and keeping only the first occurrence of each unique `(i, j)` pixel coordinate.

# Arguments
- `df::DataFrame`: A DataFrame containing at least the columns `t`, `i`, and `j`, representing time and pixel coordinates.
- `t_range::UnitRange{Int}`: A range of time frames to scan (e.g., `1:5`).

# Returns
- A new `DataFrame` containing rows from the input `df` with no repeated `(i, j)` coordinates across the specified time frames.
  If a coordinate appears multiple times in different frames, only the **earliest** (lowest `t`) instance is kept.

# Behavior
- Rows are included in the output in the order of increasing `t`.
- Once a coordinate `(i, j)` is seen, it is excluded from all future frames in `t_range`.
"""
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

end  #end module 