"""
Analysis Functions for Image Processing and Peak Detection

This module contains functions for analyzing fluorescence microscopy data,
including peak detection and local maxima identification.
"""

using DataFrames
using Images
using SparseArrays

export detect_local_maxima

"""
    detect_local_maxima(frame::AbstractMatrix{<:Real}; 
                        threshold=0.0, dx=0, dy=0, 
                        window_size=(3,3)) -> DataFrame

Detect local maxima in a 2D image frame using a sliding window approach.

# Arguments
- `frame`: 2D image matrix
- `threshold`: Minimum intensity a peak must exceed (default: 0.0)
- `dx`: Margin to exclude on left and right edges (default: 0)
- `dy`: Margin to exclude on top and bottom edges (default: 0)
- `window_size`: Size of sliding window for maxima detection (default: (3,3))

# Returns
DataFrame with columns: i (row), j (column), intensity

# Notes
- Uses replicate padding to handle borders during window operations
- Excludes peaks within dx/dy pixels of frame edges
- Window size must be odd in both dimensions

# Examples
```julia
# Create a test image with peaks
frame = [1.0 2.0 1.0; 2.0 5.0 2.0; 1.0 2.0 1.0]
peaks = detect_local_maxima(frame; threshold=3.0)

# Use custom window size
peaks = detect_local_maxima(frame; threshold=2.0, window_size=(5,5))
```
"""
function detect_local_maxima(frame::AbstractMatrix{<:Real}; 
                            threshold=0.0, dx=0, dy=0, 
                            window_size=(3,3))
    
    # Validate window size
    if any(iseven.(window_size))
        throw(ArgumentError("Window size must be odd in both dimensions"))
    end
    
    # Calculate center index for the window
    center_idx = div(prod(window_size) + 1, 2)
    
    # Find local maxima using sliding window
    is_max = mapwindow(x -> x[center_idx] == maximum(x), frame, window_size; 
                       border=Pad(:replicate))
    
    # Apply threshold and find candidates
    candidates = findall(is_max .& (frame .> threshold))
    
    # Pre-filter candidates based on edge constraints
    height, width = size(frame)
    valid_candidates = filter(candidates) do I
        i, j = Tuple(I)
        return i > dy && j > dx && i <= height - dy && j <= width - dx
    end
    
    # Extract coordinates and intensities
    n_peaks = length(valid_candidates)
    i_vals = Vector{Int}(undef, n_peaks)
    j_vals = Vector{Int}(undef, n_peaks)
    intensities = Vector{Float64}(undef, n_peaks)
    
    for (idx, I) in enumerate(valid_candidates)
        i, j = Tuple(I)
        i_vals[idx] = i
        j_vals[idx] = j
        intensities[idx] = frame[I]
    end
    
    return DataFrame(i = i_vals, j = j_vals, intensity = intensities)
end