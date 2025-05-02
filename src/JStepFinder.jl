module JStepFinder
export stepfindcore, auto_step_main
using Statistics
using LinearAlgebra
include("SimpleLogger.jl")
using .SimpleLogger
global_log_level = ERROR


"""
    auto_step_main(dataX; tresH=0.15, N_iter=10, tol=1e-3, max_iter=5)

Perform iterative step-finding on `dataX` until convergence.

# Arguments
- `dataX::AbstractVector{<:Real}`: The signal to analyze.
- `tresH::Real`: Step detection threshold.
- `N_iter::Int`: Maximum number of iterations inside `stepfindcore`.
- `tol::Real`: Relative tolerance for convergence.
- `max_iter::Int`: Maximum number of outer iterations.

# Returns
- `S_curve`: Last score curve from step detection.
- `best_shot::Int`: Best iteration index in last step detection.
- `Fit_final`: Final reconstructed fit to the data.
"""
function auto_step_main(dataX::AbstractVector{T}; 
                        tresH::Real=0.15, N_iter::Int=10, 
                        tol::Real=1e-3, max_iter::Int=5) where {T<:Real}

    #Fit_current = zeros(T, length(dataX))
    #last_best_shot = 0
    #last_S_curve = zeros(T, N_iter)  # Preallocate a dummy S_curve in case stepfindcore fails

    FitX, _, _, S_curve, best_shot = stepfindcore(dataX; tresH, N_iter)
    iter = -10
	cc = -10.0
    
    if best_shot == 0    
        return S_curve, best_shot, FitX, iter, cc
    end
    
    Fit_current = FitX
    last_best_shot = best_shot
    last_S_curve = S_curve

    debug(" best_shot = $(best_shot)")
    cc =-1.0
    for iter in 1:max_iter

        debug(" iter = $(iter)")
        residuX = dataX .- Fit_current

        debug(" R - D = $(sum(residuX))")
        Correction, _, _, S_curve, best_shot = stepfindcore(residuX; tresH, N_iter)

        if best_shot > 0
            debug("best_shot of Correction =$(best_shot)")
        
            debug(" C - D = $(sum(dataX - Correction))")
            Fit_updated = append_fitx(Correction, Fit_current, dataX; Nmin=10)
            #Fit_updated = append_fitx(residuX, Fit_current, dataX; Nmin=5)

            debug(" U - D = $(sum( dataX - Fit_updated))")
        else
            # No good step found, return current state
            debug("No good step found, return current state")
            return last_S_curve, last_best_shot, Fit_current, iter, cc
        end

        # Convergence check
        cc = norm(Fit_updated - Fit_current) / (norm(Fit_current) + eps())
        debug("CC = $(cc)")
        if cc < tol
            debug("achieved convergence L - I =$(sum(Fit_updated - Fit_current ))")
            # eps() prevents division by zero
            last_best_shot = best_shot
            last_S_curve = S_curve
            Fit_current = Fit_updated
            return last_S_curve, last_best_shot, Fit_current, iter, cc
        end

        # Prepare for next iteration
        last_best_shot = best_shot
        last_S_curve = S_curve
        Fit_current = Fit_updated
    end

    return last_S_curve, last_best_shot, Fit_current, iter, cc
end



"""
    stepfindcore(dataX::AbstractVector{T}; tresH=0.05, N_iter=0) where {T<:Real}

Performs an iterative step fitting process to find step locations in a 1D signal `dataX`.
Returns the best fit, counterfit, split log, final S-curve, and number of best steps.
"""


function stepfindcore(dataX::AbstractVector{T}; tresH::Real=0.15, N_iter::Int=10) where {T<:Real}
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


"""
    AppendFitX(newFitX::AbstractVector{T}, FitX::AbstractVector{T}, dataX::AbstractVector{T}) where {T<:Real}

Combine two fit traces (`FitX` and `newFitX`) into a single consistent stepwise signal.

# Behavior
- Merges the two fits by summing them.
- Detects step changes in the combined fit.
- If two step changes are closer together than a minimum number of points (`Nmin=2`), 
  re-fits that segment using `splitFast` to enforce a clean step.

# Arguments
- `newFitX`: New fit trace (1D vector).
- `FitX`: Existing fit trace (1D vector).
- `dataX`: Original data trace (needed to redo local fits when step locations are too close).

# Returns
- `combiFitX`: Combined and corrected fit trace.

"""
function AppendFitX(newFitX::AbstractVector{T}, FitX::AbstractVector{T}, dataX::AbstractVector{T};
    Nmin=2) where {T<:Real}
    # Combine the two fits
    combiFitX = FitX .+ newFitX
    Lx = length(combiFitX)

    # Find where steps happen (where the value changes)
    ixes0 = findall(diff(combiFitX) .!= 0)

    if !isempty(ixes0)
        # Pad with start and end to define segments
        ixes = vcat(-1, ixes0, Lx)

        # Find where steps are too close (less than Nmin points apart)
        whereblips = findall(diff(ixes) .< Nmin)

        # For each problematic region, re-fit using splitFast
        for ix in whereblips
            lo = ix - 1
            ixlo = ixes[lo]
            ixhi = ixes[lo + 3]

            segment = dataX[(ixlo+2):(ixhi)]  # Julia is 1-based indexing
            idx, avl, avr, rankit, errorcurve = splitFast(segment)

            # Apply the new local fit
            combiFitX[(ixlo+2):(ixlo+1+idx+1)] .= avl
            combiFitX[(ixlo+2+idx):(ixhi)] .= avr
        end
    end

    return combiFitX
end


"""
    append_fitx(newFitX::AbstractVector{T}, FitX::AbstractVector{T}, dataX::AbstractVector{T}; Nmin::Int=2) where {T<:Real}

Merge two fitted traces (`FitX` and `newFitX`) and remove steps that are closer together than `Nmin` points.  
For each close pair of steps, re-fit the region using `splitFast`.

# Arguments
- `newFitX`, `FitX`: Step-fitted traces of same length.
- `dataX`: Original signal used for refitting.
- `Nmin`: Minimum allowed spacing between consecutive steps (default: 2).

# Returns
- `combiFitX`: Combined step-fit trace, with merged and re-fitted regions where needed.
"""
function append_fitx(newFitX::AbstractVector{T}, 
                     FitX::AbstractVector{T}, dataX::AbstractVector{T}; Nmin::Int=2) where {T<:Real}
    # Combine fits
    combiFitX = newFitX .+ FitX
    Lx = length(combiFitX)

    # Find step changes (where signal changes)
    diffs = diff(combiFitX)
    ixes0 = findall(!≈(0.0), diffs) .+ 1  # step indices
    if isempty(ixes0)
        return combiFitX
    end

    # Pad start and end for indexing
    ixes = vcat(0, ixes0, Lx)

    # Identify segments that are too short
    seg_lengths = diff(ixes)
    close_idx = findall(<(Nmin), seg_lengths)

    for ci in close_idx
        # Re-fit over a region that spans the close segment and its neighbors
        lo = max(ci - 1, 1)
        hi = min(ci + 2, length(ixes) - 1)
        ixlo = ixes[lo] + 1
        ixhi = ixes[hi]

        # Extract segment from original data
        segment = dataX[ixlo:ixhi]
        idx, avl, avr, _, _ = splitFast(segment)

        # Replace that segment in the combined fit
        combiFitX[ixlo:ixlo+idx-1] .= avl
        combiFitX[ixlo+idx:end][1:ixhi - (ixlo + idx) + 1] .= avr
    end

    return combiFitX
end



end # module
