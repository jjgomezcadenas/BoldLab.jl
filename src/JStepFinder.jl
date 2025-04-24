module JStepFinder
export stepfindcore
using Statistics

"""
    stepfindcore(dataX::AbstractVector{T}; tresH=0.05, N_iter=0) where {T<:Real}

Performs an iterative step fitting process to find step locations in a 1D signal `dataX`.
Returns the best fit, counterfit, split log, final S-curve, and number of best steps.
"""
function stepfindcore(dataX::AbstractVector{T}; tresH=0.05, N_iter=0) where {T<:Real}
    N_iter = (N_iter == 0 || N_iter > length(dataX) ÷ 4) ? length(dataX) ÷ 4 : N_iter
    S_curve, splitlog = A121_split_until_ready(dataX, N_iter)
    best_shot, S_curve_fin = A122_eval_Scurve(S_curve, tresH)

    if best_shot > 0
        indices_bestfit, indices_counterfit = Splitlog2FinalIndices(splitlog, best_shot)
        FitX = Indices2Fit(dataX, indices_bestfit)
        cFitX = Indices2Fit(dataX, indices_counterfit)
    else
        FitX = zeros(length(dataX))
        cFitX = zeros(length(dataX))
        best_shot = 0
    end
    return FitX, cFitX, splitlog, S_curve_fin, best_shot
end


"""
    A121_split_until_ready(dataX, N_iter)

Performs iterative segmentation and step placement.
"""

function A121_split_until_ready(dataX::AbstractVector{T}, N_iter::Int) where {T<:Real}
    split_table = A1210_set_up_splitlogtable(dataX)
    splitlog = zeros(Int, N_iter, 3)
    FitX = fill(mean(dataX), length(dataX))
    cFitX = zeros(length(dataX))

    i_start, i_stop, i_next = Int.(split_table[2, 1:3])
    avl, avr = split_table[2, 4:5]
    cFitX[i_start+1:i_next+1] .= avl
    cFitX[i_next+2:i_stop+1] .= avr

    S_curve = zeros(N_iter)
    c = 1

    for _ in 1:N_iter
        seg_lengths = split_table[:, 2] .- split_table[:, 1] .+ 1
        rankings = split_table[:, 6]
        valid_idxs = findall(i -> seg_lengths[i] > 2 && rankings[i] > 0, 1:size(split_table, 1))
        isempty(valid_idxs) && break

        best_row_idx = valid_idxs[argmax(rankings[valid_idxs])]
        FitX = A1211_adapt_fit(FitX, split_table, best_row_idx)
        split_table, splitlog_entry = A1212_adapt_splitlog_table(dataX, split_table, best_row_idx)
        cFitX = A1213_adapt_counterfit(cFitX, dataX, split_table, best_row_idx)
        splitlog[c, :] = splitlog_entry
        S_curve[c] = mean((dataX .- cFitX).^2) / mean((dataX .- FitX).^2)
        c += 1
    end

    return S_curve, splitlog
end


"""
    A1210_set_up_splitlogtable(dataX)

Initializes the split table with dummy boundaries and the full trace.
"""

function A1210_set_up_splitlogtable(dataX::AbstractVector{T}) where {T<:Real}
    Na = length(dataX)
    i_nxt, avl, avr, rankit, _ = splitFast(dataX)
    new_row = [0, Na-1, i_nxt, avl, avr, rankit]
    return vcat([-1 -1 -1 0.0 0.0 0.0], new_row', [Na Na Na 0.0 0.0 0.0])
end


"""
    A122_eval_Scurve(S_curve, acceptance)

Returns the optimal step count based on the S-curve.
"""

function A122_eval_Scurve(S_curve::AbstractVector{T}, acceptance=0.15) where {T<:Real}
    S1 = max.(S_curve, 1.0)
    LS = length(S_curve)
    baseline = range(1, stop=S_curve[LS], length=LS)
    S_curve_fin = S1 .- baseline
    i_max = argmax(S_curve_fin)
    best_shot = S_curve_fin[i_max] > acceptance ? i_max : 0
    return best_shot, S_curve_fin
end


"""
    Splitlog2FinalIndices(splitlog, best_shot)

Extract best-fit and counter-fit indices from split log.
"""
function Splitlog2FinalIndices(splitlog, best_shot)
    indices_bestfit = sort(splitlog[1:best_shot, 1])
    all_counterfit = vec(splitlog[1:best_shot, 2:3])
    used_mask = in.(all_counterfit, Ref(indices_bestfit))
    indices_counterfit = sort(unique(all_counterfit[.!used_mask]))
    return indices_bestfit, indices_counterfit
end


"""
    Indices2Fit(dataX, indices)

Builds piecewise constant fit from data and indices.
"""

 function Indices2Fit(dataX::AbstractVector{T}, indices::Vector{Int}) where {T<:Real}
     FitX = similar(dataX)
    bounds = [0; indices; length(dataX)]
     for i in 1:length(bounds)-1
         lo = bounds[i]+1
         hi = bounds[i+1]
         FitX[lo:hi] .= mean(dataX[lo:hi])
     end
     return FitX
 end


 """
     splitFast(segment)
 
 Determines the best step fit in a 1D signal segment.
 Returns the best index to split, left and right means, ranking score, and error curve.
 """
 function splitFast(segment::AbstractVector{T}) where {T<:Real}
     Ns = length(segment)
     Nmin = 2
     invalidFit = false
 
     if Ns > 2
         var_q = ones(Ns)
         avl = segment[1]
         avr = sum(segment[2:end]) / (Ns - 1)
         ava = sum(segment) / Ns
 
         for ii in 2:Ns-1
             n_L = ii - 1
             n_R = Ns - ii + 1
 
             # update running averages
             avl = (avl * (n_L - 1) + segment[ii]) / n_L
             avr = (avr * (n_R + 1) - segment[ii]) / n_R
 
             delta_l = avl - ava
             delta_r = avr - ava
 
             # variance correction neutralized (varcor = 1)
             delta_q = delta_l^2 * n_L + delta_r^2 * n_R
             var_q[ii] = -delta_q
         end
 
         idx = argmin(var_q)
         if idx < Nmin || Ns - idx < Nmin
             invalidFit = true
         else
             avl_fin = mean(segment[1:idx])
             avr_fin = mean(segment[idx+1:end])
             rankit = (avr_fin - avl_fin)^2 * Ns
             errorcurve = var_q / Ns
         end
     else
         invalidFit = true
     end
 
     if invalidFit
         idx = 1
         avl_fin = segment[1]
         avr_fin = segment[1]
         rankit = 0.0
         errorcurve = zeros(length(segment))
     end
 
     return idx - 1, avl_fin, avr_fin, rankit, errorcurve
 end
 

"""
    A1211_adapt_fit(FitX, split_table, idx)

Adapt the fit trace using the split table at the given index.
"""

function A1211_adapt_fit(FitX::AbstractVector{T1}, 
                         split_table::AbstractMatrix{T2}, idx::Int) where {T1<:Real, T2<:Real}
    rank = split_table[idx, 6]
    if rank > 0
        leftplateau_istart = Int(split_table[idx, 1]) + 1
        leftplateau_istop = Int(split_table[idx, 3]) + 1
        rightplateau_istart = leftplateau_istop + 1
        rightplateau_istop = Int(split_table[idx, 2]) + 1
        leftplateau_level = split_table[idx, 4]
        rightplateau_level = split_table[idx, 5]
        FitX[leftplateau_istart:leftplateau_istop] .= leftplateau_level
        FitX[rightplateau_istart:rightplateau_istop] .= rightplateau_level
    end
    return FitX
end


"""
    A1212_adapt_splitlog_table(segm, split_table, idx)
 
Updates the split table by replacing one row with two new ones after performing a split.
Returns the updated split table and a log entry of the step split.
"""

function A1212_adapt_splitlog_table(segm::AbstractVector{T1}, 
                                    split_table::AbstractMatrix{T2}, idx::Int) where {T1<:Real, T2<:Real}
    best_row_entry = vec(split_table[idx, :])

    # LEFT plateau
    istart_1 = Int(best_row_entry[1])
    istop_1 = Int(best_row_entry[3])
    inxt_1, avl_1, avr_1, rankit_1, _ = splitFast(segm[istart_1+1:istop_1+1])
    new_row1 = [istart_1, istop_1, inxt_1 + istart_1, avl_1, avr_1, rankit_1]

    # RIGHT plateau
    istart_2 = Int(best_row_entry[3] + 1)
    istop_2 = Int(best_row_entry[2])
    inxt_2, avl_2, avr_2, rankit_2, _ = splitFast(segm[istart_2+1:istop_2+1])
    new_row2 = [istart_2, istop_2, inxt_2 + istart_2, avl_2, avr_2, rankit_2]

    splitlog_entry = [istop_1, inxt_1 + istart_1, inxt_2 + istart_2]

    block_before = split_table[1:idx-1, :]
    block_after = split_table[idx+1:end, :]
    updated_table = vcat(block_before, new_row1', new_row2', block_after)

    return updated_table, splitlog_entry
end

"""
A1213_adapt_counterfit(cFitX, dataX, split_table, idx)

Updates the counterfit trace for three plateaus using four adjacent split table entries.
"""
function A1213_adapt_counterfit(cFitX::AbstractVector{T1}, 
                                dataX::AbstractVector{T2}, 
                                split_table::AbstractMatrix{T3}, 
                                idx::Int) where {T1<:Real, T2<:Real, T3<:Real}
    
    # Ensure valid range
    if idx ≤ 1 || idx + 2 > size(split_table, 1)
        return cFitX
    end

    rankings = split_table[idx:idx+2, 6]
    ok2adapt = all(rankings .> 0)

    if ok2adapt
        left_start = Int(split_table[idx - 1, 3]) + 1
        left_end   = Int(split_table[idx, 3])
        center_start = left_end + 1
        center_end  = Int(split_table[idx + 1, 3])
        right_start = center_end + 1
        right_end   = Int(split_table[idx + 2, 3])

        if 1 ≤ left_start ≤ left_end ≤ length(dataX)
            left_level = mean(view(dataX, left_start:left_end))
            cFitX[left_start:left_end] .= left_level
        end

        if 1 ≤ center_start ≤ center_end ≤ length(dataX)
            center_level = mean(view(dataX, center_start:center_end))
            cFitX[center_start:center_end] .= center_level
        end

        if 1 ≤ right_start ≤ right_end ≤ length(dataX)
            right_level = mean(view(dataX, right_start:right_end))
            cFitX[right_start:right_end] .= right_level
        end
    end

    return cFitX
end

end # module
