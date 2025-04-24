module LabStepAnalysis
using Revise
using StatsPlots
using Statistics
using StatsBase
import Measures
using SparseArrays
using DataFrames
using DelimitedFiles
using Images, FileIO, ImageIO

using NPZ


export load_tif_stack_int16, get_tif_image, tif_to_matrix
export regularize_img, regularize_stack!
# plt_ft_trace
"""
    load_tif_stack_int16(dir::AbstractString) -> Array{Float32,3}

Given a directory `dir`, finds all TIFF files of the form `PREFIX_NNNNN.tif`, sorts them by
the integer suffix `NNNNN`, loads each as a Gray image, scales its pixel‐values into
the full `Int16` range, then casts to `Float32` and stacks into a 3‐D array
(height × width × n_frames).

# Arguments
- `dir::AbstractString`
    Path to the folder containing your TIFF frames.

# Returns
- `Array{Float32,3}`  
    A 3‑D array where `[:,:,k]` is frame k, and each pixel has been scaled
    to the full `Int16` range then converted to `Float32`.
"""
function load_tif_stack_int16(dir::AbstractString; pedestal=1600)
    # 1) get all .tif files with an underscore
    files = filter(f -> occursin('_', f) && endswith(f, ".tif"), readdir(dir))

    # 2) parse out the integer suffix after the underscore (before ".tif")
    numbered = [(fname,
                 parse(Int, splitext(split(fname, "_")[2])[1]))
                for fname in files]

    # 3) sort by that integer
    sort!(numbered, by = x -> x[2])

    # 4) load each, scale to Int16, then convert to Float32
    max16 = typemax(UInt16)
    frames = Vector{Array{Float32,2}}(undef, length(numbered))
    for (idx, (fname, _)) in enumerate(numbered)
        img_gray = load(joinpath(dir, fname))            # Gray{N0f16} array
        img_f     = Float32.(img_gray)                   # [0,1] floats
        img_i16   = UInt16.(round.(max16 .* img_f))  .- pedestal     # [0..32767] Int16
        frames[idx] = Float32.(img_i16)                  # cast to Float32
    end

    # 5) stack into a 3‑D array
    return cat(frames...; dims = 3)
end


"""
    get_image(folder::String, imglbl,  n::Int)

Loads an image from a given folder where filenames follow the pattern "Image1_<n>.ext".

# Arguments
- `folder`: Directory containing the image files.
- `n`: Integer index of the image to load.

# Returns
- A 2D or 3D array (`Matrix{Int32}` or `Array{Int32,3}`) representing the image.
"""
function get_tif_image(folder::String, imglbl, n::Int)
    nstr = string(n)
    files = readdir(folder)
    
    # Filter files starting with "Image"
    files = filter(f -> startswith(f, imglbl), files)
	
    filename = filter(f -> parse(Int, split(split(f, "_")[2], ".")[1]) == n, files)
	println(filename)

    if isempty(filename)
        Base.error("No image file found for index $n in folder $folder")
    end

    fpath = joinpath(folder, filename[end])
	println(fpath)
    img = load(fpath)  # Load the image using FileIO
	#arrf  = Float32.(channelview(img)[1,:,:])
	#scale = typemax(Int32)
	#return Int32.(round.(scale .* arrf))
    #return Float64.(Gray.(img)) 
	img
end


"""
    tif_to_matrix(path::AbstractString; T=UInt16)

Load a single‐channel TIFF and return a matrix of type `T`, rescaled from fixed‐point [0,1].
- `T=UInt8` → 0–255  
- `T=UInt16` → 0–65535  
- `T=Int32` → 0–typemax(Int32)
"""
function tif_to_matrix(img; T::Type=Int16)
    
    f     = Float32.(img)            # convert each Gray{N0f16}→Float32 in [0,1]
    scale = typemax(T)
    i     = round.(scale .* f)       # still a Matrix{Float32}
    return i                     # cast elementwise to Matrix{T}, same dims
end

"""
    regularize_stack!(imst::AbstractArray{T,3}; nsigma::Real=10) where {T<:Real}

In-place version of `regularize_stack`. Replaces "noisy" pixels across all frames
with the frame-wise mean, based on noise detected in the first frame.

# Arguments
- `imst`: A 3D array (n × m × t) of image data. This array is modified in-place.
- `nsigma`: Threshold multiplier for noise detection (default: 10).

# Returns
- The modified array `imst`, updated in-place.
"""
function regularize_stack!(imst::AbstractArray{T,3}; nsigma::Real=10) where {T<:Real}
    n, m, t = size(imst)

    frame1 = imst[:, :, 1]
    μ = mean(frame1)
    σ = std(frame1)
    threshold = μ + nsigma * σ
    mask = frame1 .> threshold

    for k in 1:t
        μk = mean(imst[:, :, k])
        for i in 1:n, j in 1:m
            if mask[i, j]
                imst[i, j, k] = μk
            end
        end
    end

    return mask
end


"""
Mask spots which are too bright (>nsigma sigma over the mean).
"""
function regularize_img(img::AbstractArray{T,2}; nsigma::Real=10) where {T<:Real}
	ximg = mean(img)
	simg = std(img)
	dimg =size(img)
	img2 = copy(img)
	I = Int[]
	J = Int[]
	for i in 1:dimg[1]
		for j in 1:dimg[2]
			if img[i,j] > ximg + nsigma * simg 
				push!(I,i)
				push!(J,j)
				img2[i,j] = ximg  
			end
		end
	end
	img2, CartesianIndex.(I, J)
end



end
# end module
