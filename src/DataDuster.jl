# DataDuster.jl
# Julia translation of MATLAB DataDuster.m core logic
# Removes NaN, Inf, and obvious artifacts from data vectors

module DataDuster

using Statistics

export clean_data

"""
    clean_data(data::Vector{<:Real}) -> Vector{Float64}

Cleans a numeric data vector by removing NaNs, Infs, and isolated spikes.
"""
function clean_data(data::Vector{<:Real})
    # Remove NaNs and Infs
    finite_data = filter(isfinite, data)

    # Optional: remove large outlier spikes (if desired)
    # Example: spike filtering using rolling window or z-score
    cleaned = copy(finite_data)
    N = length(cleaned)
    if N < 3
        return cleaned
    end

    # Simple spike suppression using a local median filter
    for i in 2:(N - 1)
        neighborhood = [cleaned[i - 1], cleaned[i + 1]]
        local_median = median(neighborhood)
        if abs(cleaned[i] - local_median) > 10 * std(neighborhood)
            cleaned[i] = local_median
        end
    end

    return cleaned
end

end # module
