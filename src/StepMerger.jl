# StepMerger.jl
# Translated from MATLAB StepMerger.m
# Merges short dwell steps from an existing step fit

module StepMerger

using Statistics

export merge_short_dwells

"""
    merge_short_dwells(data::Vector{<:Real}, fit::Vector{<:Real}; min_dwell::Int=2) -> Vector{Float64}

Merges short plateaus in `fit` that are shorter than `min_dwell`.
Returns a new, smoothed step fit.
"""
function merge_short_dwells(data::Vector{<:Real}, fit::Vector{<:Real}; min_dwell::Int=2)
    N = length(data)
    diffs = diff(fit)
    change_pts = findall(!iszero, diffs)
    boundaries = [0; change_pts; N]

    new_fit = similar(fit)
    i = 1
    while i < length(boundaries)
        start_idx = boundaries[i] + 1
        stop_idx = boundaries[i+1]
        dwell = stop_idx - start_idx + 1

        if dwell < min_dwell && i > 1 && i < length(boundaries) - 1
            # Merge with previous or next depending on similarity
            prev_level = mean(data[boundaries[i-1]+1:boundaries[i]])
            next_level = mean(data[boundaries[i+1]+1:boundaries[i+2]])
            curr_level = mean(data[start_idx:stop_idx])
            if abs(curr_level - prev_level) < abs(curr_level - next_level)
                new_fit[start_idx:stop_idx] .= prev_level
            else
                new_fit[start_idx:stop_idx] .= next_level
            end
        else
            # Keep current plateau
            level = mean(data[start_idx:stop_idx])
            new_fit[start_idx:stop_idx] .= level
        end

        i += 1
    end

    return new_fit
end

end # module