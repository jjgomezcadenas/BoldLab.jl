using DataFrames
using Plots
using ImageFiltering
using Statistics
# ============================================================================
# Diffusion Functions for ITACA Detector Simulation
# ============================================================================

"""
    GridConfig

Configuration for histogram grid setup.
Fields:
- x_min, x_max, y_min, y_max: Grid boundaries
- nbins_x, nbins_y: Number of bins in each dimension
- dx, dy: Bin widths
"""
struct GridConfig
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    nbins_x::Int
    nbins_y::Int
    dx::Float64
    dy::Float64
end

"""
    extract_hits_and_electrons(df::DataFrame) -> (x_hits, y_hits, electrons)

Extract hit positions and electron counts from DataFrame.
If electrons column doesn't exist, assumes 1 electron per hit.
"""
function extract_hits_and_electrons(df::DataFrame)
    x_hits = Float64.(df.x)
    y_hits = Float64.(df.y)

    if hasproperty(df, :electrons)
        electrons = Int.(df.electrons)
    else
        @warn "DataFrame does not have 'electrons' column, using 1 electron per hit"
        electrons = ones(Int, length(x_hits))
    end

    return x_hits, y_hits, electrons
end

"""
    setup_grid_variable_bins(x_hits, y_hits; nbins=100, nsigma=3.0, sigma_mm=0.6) -> GridConfig

Setup grid with variable bin size based on data range and requested number of bins.
"""
function setup_grid_variable_bins(x_hits, y_hits; nbins=100, nsigma=3.0, sigma_mm=0.6)
    # Calculate histogram range with padding
    x_min = minimum(x_hits) - nsigma * sigma_mm
    x_max = maximum(x_hits) + nsigma * sigma_mm
    y_min = minimum(y_hits) - nsigma * sigma_mm
    y_max = maximum(y_hits) + nsigma * sigma_mm

    # Calculate bin widths
    dx = (x_max - x_min) / nbins
    dy = (y_max - y_min) / nbins

    return GridConfig(x_min, x_max, y_min, y_max, nbins, nbins, dx, dy)
end

"""
    setup_grid_fixed_pitch(x_hits, y_hits; pitch_mm=1.0, nsigma=3.0, sigma_mm=0.6) -> GridConfig

Setup grid with fixed pitch (bin size), aligning to exact pitch boundaries.
Useful for matching physical detector geometry (e.g., 1 mm ion collection grid).
"""
function setup_grid_fixed_pitch(x_hits, y_hits; pitch_mm=1.0, nsigma=3.0, sigma_mm=0.6)
    # Calculate data range with padding
    x_min_data = minimum(x_hits) - nsigma * sigma_mm
    x_max_data = maximum(x_hits) + nsigma * sigma_mm
    y_min_data = minimum(y_hits) - nsigma * sigma_mm
    y_max_data = maximum(y_hits) + nsigma * sigma_mm

    # Align grid to pitch_mm boundaries (round outward)
    x_min = floor(x_min_data / pitch_mm) * pitch_mm
    x_max = ceil(x_max_data / pitch_mm) * pitch_mm
    y_min = floor(y_min_data / pitch_mm) * pitch_mm
    y_max = ceil(y_max_data / pitch_mm) * pitch_mm

    # Calculate number of bins (exact multiples of pitch)
    nbins_x = Int(round((x_max - x_min) / pitch_mm))
    nbins_y = Int(round((y_max - y_min) / pitch_mm))

    # Verify pitch
    dx = (x_max - x_min) / nbins_x
    dy = (y_max - y_min) / nbins_y
    @assert abs(dx - pitch_mm) < 1e-10 "X bin width $(dx) != $(pitch_mm) mm"
    @assert abs(dy - pitch_mm) < 1e-10 "Y bin width $(dy) != $(pitch_mm) mm"

    #@info "Grid setup: $(nbins_x) × $(nbins_y) bins with $(pitch_mm) mm pitch"
    #@info "X range: [$(x_min), $(x_max)] mm"
    #@info "Y range: [$(y_min), $(y_max)] mm"

    return GridConfig(x_min, x_max, y_min, y_max, nbins_x, nbins_y, dx, dy)
end

"""
    diffuse_and_histogram!(hist_counts, x_hits, y_hits, electrons, grid::GridConfig, sigma_mm)

Core diffusion algorithm: generate diffused electrons and bin them into histogram.
Returns the number of electrons captured within the histogram bounds.
"""
function diffuse_and_histogram!(hist_counts, x_hits, y_hits, electrons, grid::GridConfig, sigma_mm)
    total_electrons = sum(electrons)
    electrons_counted = 0

    @inbounds for i in eachindex(x_hits)
        x_center = x_hits[i]
        y_center = y_hits[i]
        n_electrons = electrons[i]

        # Generate diffused positions for each electron (2D Gaussian)
        for _ in 1:n_electrons
            x_diffused = x_center + randn() * sigma_mm
            y_diffused = y_center + randn() * sigma_mm

            # Find bin indices
            ix = Int(floor((x_diffused - grid.x_min) / grid.dx)) + 1
            iy = Int(floor((y_diffused - grid.y_min) / grid.dy)) + 1

            # Check if within bounds
            if 1 <= ix <= grid.nbins_x && 1 <= iy <= grid.nbins_y
                hist_counts[ix, iy] += 1
                electrons_counted += 1
            end
        end
    end

    return electrons_counted, total_electrons
end

"""
    histogram_to_dataframe(hist_counts, grid::GridConfig) -> DataFrame

Convert 2D histogram to DataFrame with columns: x, y, electrons.
"""
function histogram_to_dataframe(hist_counts, grid::GridConfig)
    # Calculate bin centers
    x_edges = range(grid.x_min, grid.x_max, length=grid.nbins_x+1)
    y_edges = range(grid.y_min, grid.y_max, length=grid.nbins_y+1)
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    y_centers = (y_edges[1:end-1] .+ y_edges[2:end]) ./ 2

    # Build DataFrame (preallocate for efficiency)
    n_total = grid.nbins_x * grid.nbins_y
    result_x = Vector{Float64}(undef, n_total)
    result_y = Vector{Float64}(undef, n_total)
    result_intensity = Vector{Int}(undef, n_total)

    idx = 1
    for ix in 1:grid.nbins_x
        for iy in 1:grid.nbins_y
            result_x[idx] = x_centers[ix]
            result_y[idx] = y_centers[iy]
            result_intensity[idx] = hist_counts[ix, iy]
            idx += 1
        end
    end

    return DataFrame(
        x = result_x,
        y = result_y,
        electrons = result_intensity
    )
end

"""
    histogram_to_dataframe_ions(hist_counts, grid::GridConfig) -> DataFrame

Convert 2D histogram to DataFrame with columns: x, y, n_ions (for ion collection).
"""
function histogram_to_dataframe_ions(hist_counts, grid::GridConfig)
    # Calculate bin centers
    x_edges = range(grid.x_min, grid.x_max, length=grid.nbins_x+1)
    y_edges = range(grid.y_min, grid.y_max, length=grid.nbins_y+1)
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    y_centers = (y_edges[1:end-1] .+ y_edges[2:end]) ./ 2

    # Build DataFrame (preallocate for efficiency)
    n_total = grid.nbins_x * grid.nbins_y
    result_x = Vector{Float64}(undef, n_total)
    result_y = Vector{Float64}(undef, n_total)
    result_n_ions = Vector{Int}(undef, n_total)

    idx = 1
    for ix in 1:grid.nbins_x
        for iy in 1:grid.nbins_y
            result_x[idx] = x_centers[ix]
            result_y[idx] = y_centers[iy]
            result_n_ions[idx] = hist_counts[ix, iy]
            idx += 1
        end
    end

    return DataFrame(
        x = result_x,
        y = result_y,
        n_ions = result_n_ions
    )
end

"""
    diffuse_xy_image_mc(df::DataFrame; sigma_mm=0.6, nbins=100, nsigma=3.0) -> DataFrame

Diffuse electron hits with variable bin size.
Returns DataFrame with x, y coordinates and electron counts per bin.

# Arguments
- `df`: DataFrame with x, y positions and optionally electrons column
- `sigma_mm`: Diffusion sigma in mm (default: 0.6)
- `nbins`: Number of bins in each dimension (default: 100)
- `nsigma`: Number of sigmas for padding around data range (default: 3.0)
"""
function diffuse_xy_image_mc(df::DataFrame; sigma_mm=0.6, nbins=100, nsigma=3.0)
    # Extract data
    x_hits, y_hits, electrons = extract_hits_and_electrons(df)

    # Setup grid with variable bin size
    grid = setup_grid_variable_bins(x_hits, y_hits; nbins=nbins, nsigma=nsigma, sigma_mm=sigma_mm)

    # Initialize histogram
    hist_counts = zeros(Int, grid.nbins_x, grid.nbins_y)

    # Perform diffusion and binning
    electrons_counted, total_electrons = diffuse_and_histogram!(
        hist_counts, x_hits, y_hits, electrons, grid, sigma_mm
    )

    # Optional: report capture fraction
    # capture_fraction = electrons_counted / total_electrons
    # @info "Captured $(electrons_counted) out of $(total_electrons) electrons ($(round(capture_fraction*100, digits=1))%)"

    # Convert to DataFrame
    return histogram_to_dataframe(hist_counts, grid)
end

"""
    diffuse_xy_ions_mc(df::DataFrame; sigma_mm=0.6, pitch_mm=1.0, nsigma=3.0) -> (GridConfig, DataFrame)

Diffuse hits on a fixed grid with specified pitch (e.g., 1 mm for ion collection grid).
Bins are aligned to exact pitch boundaries to match physical detector geometry.

Returns grid configuration and DataFrame with columns: x, y, n_ions

# Arguments
- `df`: DataFrame with x, y positions and optionally electrons column
- `sigma_mm`: Diffusion sigma in mm (default: 0.6)
- `pitch_mm`: Grid pitch in mm (default: 1.0)
- `nsigma`: Number of sigmas for padding around data range (default: 3.0)
"""
function diffuse_xy_ions_mc(df::DataFrame; sigma_mm=0.6, pitch_mm=1.0, nsigma=3.0)
    # Extract data
    x_hits, y_hits, electrons = extract_hits_and_electrons(df)

    # Setup grid with fixed pitch
    grid = setup_grid_fixed_pitch(x_hits, y_hits; pitch_mm=pitch_mm, nsigma=nsigma, sigma_mm=sigma_mm)

    # Initialize histogram
    hist_counts = zeros(Int, grid.nbins_x, grid.nbins_y)

    # Perform diffusion and binning
    electrons_counted, total_electrons = diffuse_and_histogram!(
        hist_counts, x_hits, y_hits, electrons, grid, sigma_mm
    )

    # Report capture fraction
    #capture_fraction = electrons_counted / total_electrons
    #@info "Captured $(electrons_counted) out of $(total_electrons) electrons ($(round(capture_fraction*100, digits=1))%)"

    # Convert to DataFrame with n_ions column
    ions_df = histogram_to_dataframe_ions(hist_counts, grid)

    return grid, ions_df
end

"""
    plot_diffused_xy(diffused_df::DataFrame; title="Diffused XY Image", titlefontsize=10, colorbar=true, intensity_col=:electrons, zero_as_nan=false)

Plot diffused electron distribution as a heatmap.

# Arguments
- `diffused_df`: DataFrame with x, y, and intensity columns (output from diffuse_xy_* functions)
- `title`: Plot title
- `titlefontsize`: Font size for title
- `colorbar`: Show colorbar (default: true)
- `intensity_col`: Column name to use for intensity values (default: :electrons)
- `zero_as_nan`: If true, replace zeros with NaN for white background (default: false)
"""
function plot_diffused_xy(diffused_df::DataFrame; title="Diffused XY Image", titlefontsize=10,
                           colorbar=true, intensity_col=:electrons, zero_as_nan=false)

    # Get intensity matrix and unique coordinates
    x_unique, y_unique, intensity_matrix = get_intensity_matrix(diffused_df, intensity_col)
    return plot_imatrix(intensity_matrix, x_unique, y_unique;
                        title=title,
                        zero_as_nan=zero_as_nan,
                        titlefontsize=titlefontsize,
                        colorbar=colorbar)

    # # Replace zeros with NaN if requested (for white background)
    # if zero_as_nan
    #     intensity_matrix[intensity_matrix .== 0] .= NaN
    # end

    # # Determine plot extent
    # xlim = (minimum(x_unique), maximum(x_unique))
    # ylim = (minimum(y_unique), maximum(y_unique))

    # # Create heatmap
    # p = heatmap(x_unique, y_unique, intensity_matrix',
    #             xlabel="X (mm)",
    #             ylabel="Y (mm)",
    #             title=title,
    #             titlefontsize=titlefontsize,
    #             colorbar=colorbar,
    #             color=:viridis,
    #             clims=(0, maximum(intensity_matrix)),
    #             xlims=xlim,
    #             ylims=ylim,
    #             xrotation=45,
    #             xtickfontsize=8,
    #             ytickfontsize=8,
    #             bottom_margin=5Plots.mm)

    # return p
end


function plot_imatrix(intensity_matrix::Matrix{Float64}, x_unique::Vector{Float64}, y_unique::Vector{Float64}; 
                      title="intensity matrix", zero_as_nan=false,
                      titlefontsize=10, colorbar=true)

    
    # Replace zeros with NaN if requested (for white background)
    if zero_as_nan
        intensity_matrix[intensity_matrix .== 0] .= NaN
    end

    # Determine plot extent
    xlim = (minimum(x_unique), maximum(x_unique))
    ylim = (minimum(y_unique), maximum(y_unique))

    # Create heatmap
    p = heatmap(x_unique, y_unique, intensity_matrix',
                xlabel="X (mm)",
                ylabel="Y (mm)",
                title=title,
                titlefontsize=titlefontsize,
                colorbar=colorbar,
                color=:viridis,
                #clims=(0, maximum(intensity_matrix)),
                xlims=xlim,
                ylims=ylim,
                xrotation=45,
                xtickfontsize=8,
                ytickfontsize=8,
                bottom_margin=5Plots.mm)

    return p
end


"""
    compute_background_from_matrix(intensity_matrix::Matrix{<:Real}; σ::Real=10.0)

Compute background from intensity matrix using Gaussian filtering.
Similar to compute_background_from_stack but for 2D matrices.

# Arguments
- `intensity_matrix`: 2D matrix of intensity values
- `σ`: Standard deviation for Gaussian kernel (default: 10.0)

# Returns
- Smoothed background matrix
"""
function compute_background_from_matrix(intensity_matrix::Matrix{<:Real}; σ::Real=10.0)
    # Apply Gaussian filter to smooth the background
    background = ImageFiltering.imfilter(intensity_matrix, ImageFiltering.Kernel.gaussian(σ))
    return background
end

"""
    subtract_background_from_matrix(intensity_matrix::Matrix{<:Real}, background::Matrix{<:Real}) -> Matrix{Float64}

Subtract a 2D background from an intensity matrix, clamping negative values to zero.
Similar to subtract_background_from_stack but for 2D matrices.

# Arguments
- `intensity_matrix::Matrix{<:Real}`: 2D matrix of intensity values
- `background::Matrix{<:Real}`: 2D matrix representing the estimated background

# Returns
- `Matrix{Float64}`: A new 2D matrix where the background has been subtracted and values are clamped to ≥ 0
"""
subtract_background_from_matrix(intensity_matrix::Matrix{<:Real},
                                 background::Matrix{<:Real}) = max.(intensity_matrix .- background, 0.0)

"""
    rebin_matrix(intensity_matrix::Matrix{<:Real}, x_unique::Vector{Float64}, y_unique::Vector{Float64}, binning::Int) -> (Matrix{Float64}, Vector{Float64}, Vector{Float64})

Rebin a 2D intensity matrix and its coordinate arrays by an integer factor, summing values in each bin.

# Arguments
- `intensity_matrix::Matrix{<:Real}`: 2D matrix to rebin
- `x_unique::Vector{Float64}`: X-coordinates corresponding to matrix columns
- `y_unique::Vector{Float64}`: Y-coordinates corresponding to matrix rows
- `binning::Int`: Rebinning factor (1 = no rebinning, 2 = 2×2 binning, etc.)

# Returns
- `(Matrix{Float64}, Vector{Float64}, Vector{Float64})`: Rebinned matrix, x_coords, y_coords

# Notes
- If dimensions are not exact multiples of binning, the matrix and coordinates are trimmed
- Values in each bin are summed (appropriate for count data)
- Coordinate arrays are rebinned by taking the center of each bin

# Example
```julia
# Get intensity matrix
x_unique, y_unique, intensity_matrix = get_intensity_matrix(ions_df, :n_ions)

# No rebinning
matrix, x, y = rebin_matrix(intensity_matrix, x_unique, y_unique, 1)

# 2×2 binning
matrix_2x, x_2x, y_2x = rebin_matrix(intensity_matrix, x_unique, y_unique, 2)

# 5×5 binning
matrix_5x, x_5x, y_5x = rebin_matrix(intensity_matrix, x_unique, y_unique, 5)
```
"""
function rebin_matrix(intensity_matrix::Matrix{<:Real},
                      x_unique::Vector{Float64},
                      y_unique::Vector{Float64},
                      binning::Int)
    if binning == 1
        return Float64.(intensity_matrix), x_unique, y_unique
    end

    rows, cols = size(intensity_matrix)

    # Calculate new dimensions (trim to fit exact multiples)
    new_rows = rows ÷ binning
    new_cols = cols ÷ binning

    # Trim matrix to fit exact multiples of binning
    rows_trim = new_rows * binning
    cols_trim = new_cols * binning
    M_trim = intensity_matrix[1:rows_trim, 1:cols_trim]

    # Rebin coordinates by averaging bins
    x_trim = x_unique[1:rows_trim]
    y_trim = y_unique[1:cols_trim]

    # Compute new coordinate arrays (center of each bin)
    x_rebinned = [mean(x_trim[i:i+binning-1]) for i in 1:binning:rows_trim]
    y_rebinned = [mean(y_trim[j:j+binning-1]) for j in 1:binning:cols_trim]

    # Reshape into (binning, new_rows, binning, new_cols)
    M_reshaped = reshape(M_trim, binning, new_rows, binning, new_cols)

    # Reorder dimensions to group bins: (new_rows, new_cols, binning, binning)
    M_reordered = permutedims(M_reshaped, (2, 4, 1, 3))

    # Sum over each bin (dims 3 and 4)
    M_binned = sum(M_reordered, dims=(3, 4))

    # Remove singleton dimensions
    M_binned = dropdims(M_binned, dims=(3, 4))

    return Float64.(M_binned), x_rebinned, y_rebinned
end


function get_intensity_matrix(diffused_df::DataFrame, intensity_col::Symbol)
    # Get unique x and y coordinates
    x_unique = sort(unique(diffused_df.x))
    y_unique = sort(unique(diffused_df.y))

    # Determine grid dimensions
    nbins_x = length(x_unique)
    nbins_y = length(y_unique)

    # Create intensity matrix properly by mapping coordinates
    # Don't assume ordering - build the matrix correctly
    intensity_matrix = zeros(Float64, nbins_x, nbins_y)

    # Create lookup dictionaries for indices
    x_to_idx = Dict(x => i for (i, x) in enumerate(x_unique))
    y_to_idx = Dict(y => j for (j, y) in enumerate(y_unique))

    # Fill in the matrix with actual data
    for row in eachrow(diffused_df)
        i = x_to_idx[row.x]
        j = y_to_idx[row.y]
        intensity_matrix[i, j] = row[intensity_col]
    end

    return x_unique, y_unique, intensity_matrix
end
