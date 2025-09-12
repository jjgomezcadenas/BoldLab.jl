module LabStepAnalysis
using Statistics
using SparseArrays
using DataFrames
using DelimitedFiles
using Images, FileIO, ImageIO


export load_tif_stack_int16, get_tif_image, tif_to_matrix
export regularize_img, regularize_stack!


# Internal helper function to process TIFF files
function _process_tif_files(file_paths::Vector{Tuple{String,Int}}; pedestal=1600)
    # Load each, scale to Int16, then convert to Float32
    max16 = typemax(UInt16)
    frames = Vector{Array{Float32,2}}(undef, length(file_paths))
    
    for (idx, (fpath, _)) in enumerate(file_paths)
        img_gray = load(fpath)                           # Gray{N0f16} array
        img_f     = Float32.(img_gray)                   # [0,1] floats
        img_i16   = UInt16.(round.(max16 .* img_f)) .- pedestal     # [0..65535] UInt16 minus pedestal
        frames[idx] = Float32.(img_i16)                  # cast to Float32
    end
    
    # Stack into a 3-D array
    return cat(frames...; dims = 3)
end

# Internal helper to parse and sort TIFF files by their numeric suffix
function _parse_and_sort_tifs(files_with_paths::Vector{Tuple{String,String}})
    numbered = Tuple{String,Int}[]
    
    for (fname, fpath) in files_with_paths
        if occursin('_', fname)
            # Extract number after last underscore
            parts = split(fname, "_")
            num_str = splitext(parts[end])[1]
            push!(numbered, (fpath, parse(Int, num_str)))
        end
    end
    
    # Sort by the parsed integer
    sort!(numbered, by = x -> x[2])
    return numbered
end
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
    # Get all .tif files with an underscore
    files = filter(f -> occursin('_', f) && endswith(f, ".tif"), readdir(dir))
    
    if isempty(files)
        error("No TIFF files with underscore found in directory: $dir")
    end
    
    # Create tuples of (filename, full_path)
    files_with_paths = [(fname, joinpath(dir, fname)) for fname in files]
    
    # Parse and sort
    numbered = _parse_and_sort_tifs(files_with_paths)
    
    # Process the files
    return _process_tif_files(numbered; pedestal=pedestal)
end

"""
    load_tif_stack_int16(filenames::Vector{String}, dir::AbstractString; pedestal=1600)

Load a stack of TIFF images from a list of filenames.

Files are loaded from the given directory, scaled from `UInt16` range 
to the full `Int16` range then converted to `Float32`.

# Arguments
- `filenames`: Vector of TIFF filenames to load
- `dir`: Directory path where the files are located
- `pedestal`: Value to subtract from pixel values (default: 1600)
"""
function load_tif_stack_int16(filenames::Vector{String}, dir::AbstractString; pedestal=1600)
    # Filter to ensure we only have .tif files
    tif_files = filter(f -> endswith(lowercase(f), ".tif") || endswith(lowercase(f), ".tiff"), filenames)
    
    if isempty(tif_files)
        error("No TIFF files found in the provided list")
    end
    
    # Create tuples of (filename, full_path)
    files_with_paths = [(fname, joinpath(dir, fname)) for fname in tif_files]
    
    # Parse and sort
    numbered = _parse_and_sort_tifs(files_with_paths)
    
    # Process the files
    return _process_tif_files(numbered; pedestal=pedestal)
end

"""
    load_tif_stack_int16(filenames::Vector{String}; pedestal=1600)

Load a stack of TIFF images from a list of full file paths.

Files are loaded directly from their full paths, scaled from `UInt16` range 
to the full `Int16` range then converted to `Float32`.

# Arguments
- `filenames`: Vector of full paths to TIFF files
- `pedestal`: Value to subtract from pixel values (default: 1600)
"""
function load_tif_stack_int16(filenames::Vector{String}; pedestal=1600)
    # Filter to ensure we only have .tif files
    tif_files = filter(f -> endswith(lowercase(f), ".tif") || endswith(lowercase(f), ".tiff"), filenames)
    
    if isempty(tif_files)
        error("No TIFF files found in the provided list")
    end
    
    # Create tuples of (basename, full_path) for parsing
    files_with_paths = [(basename(fpath), fpath) for fpath in tif_files]
    
    # Parse and sort
    numbered = _parse_and_sort_tifs(files_with_paths)
    
    # Process the files
    return _process_tif_files(numbered; pedestal=pedestal)
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
