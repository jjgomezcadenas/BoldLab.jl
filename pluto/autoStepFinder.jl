### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 98eb0dbc-0da3-11f0-3773-ef1c04c57127
using Pkg; Pkg.activate(ENV["JBoldLab"])

# ╔═╡ 54acbf3b-d144-49cc-a668-bb88566d6f85
begin
	using Revise
	using PlutoUI
	using CSV
	using DataFrames
	using Plots 
	using Measures
	using Printf
	using InteractiveUtils
	using LsqFit
	using Statistics
	using LinearAlgebra
	using StatsBase
	using DelimitedFiles
	using Images, FileIO
	using BoldLab
	using Random
end

# ╔═╡ b571e2a9-5e93-4ac5-8af1-d8c66466336c
names(BoldLab)

# ╔═╡ e121bde7-2005-4c22-80b8-de679cb80edd
function A1211_adapt_fit(FitX::Vector, split_table, idx::Int)
    rank = split_table[idx, 6]
    if rank > 0
        istart = Int(round(split_table[idx, 1]))
        istop  = Int(round(split_table[idx, 3]))
        iend   = Int(round(split_table[idx, 2]))
        FitX[istart:istop] .= split_table[idx, 4]
        FitX[istop+1:iend] .= split_table[idx, 5]
    end
    return FitX
end

# ╔═╡ 58b3137b-7e1f-45d0-918e-93f10cb55dea
function A1213_adapt_counterfit(cFitX::Vector, dataX::Vector, split_table, idx::Int)
    rankings = split_table[idx:idx+1, 6]
    if all(rankings .> 0)
        i1 = Int(round(split_table[idx-1, 3])) + 1
        i2 = Int(round(split_table[idx, 3]))
        i3 = Int(round(split_table[idx+1, 3]))
        i4 = Int(round(split_table[idx+2, 3]))

        L = length(dataX)

        if i1 > 0 && i2 <= L && i1 <= i2
            cFitX[i1:i2] .= mean(dataX[i1:i2])
        end
        if i2+1 > 0 && i3 <= L && i2+1 <= i3
            cFitX[i2+1:i3] .= mean(dataX[i2+1:i3])
        end
        if i3+1 > 0 && i4 <= L && i3+1 <= i4
            cFitX[i3+1:i4] .= mean(dataX[i3+1:i4])
        end
    end
    return cFitX
end

# ╔═╡ 33484d1b-58db-479c-956e-98037ed7be7e
function indices_to_fit(dataX::Vector, indices::Vector{Int}; how::Symbol=:mean)
    FitX = similar(dataX)
    indices_ext = [-1; indices; length(dataX)]
    for i in 1:(length(indices_ext) - 1)
        lo = indices_ext[i] + 1
        hi = indices_ext[i+1]
        if lo >= 1 && hi <= length(dataX) && lo <= hi
            FitX[lo:hi] .= how == :mean ? mean(dataX[lo:hi]) : median(dataX[lo:hi])
        end
    end
    return FitX
end

# ╔═╡ 1ab697d4-de94-45bd-8127-6f94258a9ec7
function fit_to_steps(dataX::Vector, FitX::Vector)
    Lx = length(dataX)
    globalnoise = std(diff(dataX .- FitX)) / √2
    ixes0 = findall(diff(FitX) .!= 0)
    ixes = [-1; ixes0; Lx]
    Lix = length(ixes) - 2
    steptable = zeros(Lix, 8)

    for ii in 1:Lix
        ix_pre = ixes[ii]
        ix = ixes[ii + 1]
        ix_aft = ixes[ii + 2]

        # Skip if any index is out of bounds
        if ix < 1 || ix + 1 > length(FitX) || ix_pre + 1 < 1 || ix_aft > length(FitX)
            continue
        end

        lev_pre = FitX[ix]
        lev_aft = FitX[ix + 1]
        step = lev_aft - lev_pre
        dwell_pre = ix - ix_pre
        dwell_aft = ix_aft - ix
        error_pred = 2 * sqrt(globalnoise^2 / dwell_pre + globalnoise^2 / dwell_aft) / √2

        rms_pre = (ix_pre + 1 <= ix && ix <= length(dataX)) ? std(dataX[(ix_pre + 1):ix]) : 0.0
        rms_aft = (ix + 1 <= ix_aft && ix_aft <= length(dataX)) ? std(dataX[(ix + 1):(ix_aft)]) : 0.0
        error_meas = 2 * sqrt(rms_pre^2 / dwell_pre + rms_aft^2 / dwell_aft) / √2
        steptable[ii, :] = [ix, lev_pre, lev_aft, step, dwell_pre, dwell_aft, error_pred, error_meas]
    end

    return steptable
end

# ╔═╡ d1b6be43-ff7c-4721-8657-4b0578829ad3
begin # StepFindCore



end

# ╔═╡ da6d50b4-2d05-40bc-af0d-8e715112c6bd
begin # DataDuster
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
end

# ╔═╡ ae100887-9df6-41ce-8e7a-79032c8a3ddd
begin # StepMaker
	"""
    generate_steps(n_steps::Int; 
                   minstep::Float64=1.0, 
                   maxstep::Float64=4.0, 
                   flat::Bool=true, 
                   gauss_noise::Bool=true,
                   dwell::Int=100) -> Vector{Float64}

Generate a stepwise trace with `n_steps` transitions, with random or flat step heights,
optional Gaussian noise, and constant dwell time.
"""
function generate_steps(n_steps::Int; 
                        minstep::Float64=1.0, 
                        maxstep::Float64=4.0, 
                        flat::Bool=true, 
                        gauss_noise::Bool=true,
                        dwell::Int=100)
    total_length = n_steps * dwell
    trace = Float64[]
    current_level = 0.0

    for i in 1:n_steps
        step_size = flat ? minstep : (minstep + rand() * (maxstep - minstep))
        direction = rand(Bool) ? 1 : -1
        current_level += direction * step_size
        segment = fill(current_level, dwell)
        if gauss_noise
            segment .+= randn(dwell)
        end
        append!(trace, segment)
    end

    return trace
end

end # module


# ╔═╡ a58fcf00-96a7-48d2-af21-ec316662d542
begin #StepMerger
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

# ╔═╡ 5a5fd414-a681-4d7d-9ec6-93b206146ed5
function split_fast(segment::Vector{<:Real}; demo::Bool=false)
    N = length(segment)
    if N <= 2
        return 1, segment[1], segment[1], 0.0, zeros(N)
    end

    var_q = ones(N)
    avl = segment[1]
    avr = sum(segment[2:end]) / (N - 1)
    ava = mean(segment)

    for ii in 2:N-1
        n_L = ii - 1
        n_R = N - ii + 1
        avl = (avl * (n_L - 1) + segment[ii]) / n_L
        avr = (avr * (n_R + 1) - segment[ii]) / n_R
        delta_l = avl - ava
        delta_r = avr - ava
        delta_q = delta_l^2 * n_L + delta_r^2 * n_R
        var_q[ii] = -delta_q
    end

    idx = argmin(var_q)
    if idx < 2 || N - idx < 2
        return 1, segment[1], segment[1], 0.0, zeros(N)
    end

    avl_fin = mean(segment[1:idx])
    avr_fin = mean(segment[idx+1:end])
    rankit = (avr_fin - avl_fin)^2 * N
    errorcurve = var_q ./ N
    return idx, avl_fin, avr_fin, rankit, errorcurve
end

# ╔═╡ fefc626d-3728-4973-a015-d59952007754
function A1210_set_up_splitlogtable(dataX::Vector{<:Real})
    Na = length(dataX)
    i_nxt, avl, avr, rankit, _ = split_fast(dataX)

    if rankit == 0.0 || i_nxt <= 1 || i_nxt >= Na - 1
        @warn "split_fast returned invalid result — returning fallback split_table"
        # fallback with 1-row dummy 2D matrix
        return reshape([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1, :)
    end

    new_row1 = [1.0, Na, i_nxt, avl, avr, rankit]
    split_table = vcat(
        reshape([-1.0, -1.0, -1.0, 0.0, 0.0, 0.0], 1, :),
        reshape(new_row1, 1, :),
        reshape([Na+1, Na+1, Na+1, 0.0, 0.0, 0.0], 1, :)
    )
    return split_table
end

# ╔═╡ 1881e12a-bdd3-4afd-ac2f-b5a45eb3ce80
function A1212_adapt_splitlog_table(dataX::Vector, split_table, idx::Int, demo::Float64)
    entry = split_table[idx, :]

    istart_1 = Int(round(entry[1]))
    istop_1 = Int(round(entry[3]))

    try
        if (istart_1 > 0) && (istop_1 <= length(dataX)) && (istart_1 <= istop_1)
            seg1 = dataX[istart_1:istop_1]
            i1, avl1, avr1, r1, _ = split_fast(seg1; demo=demo > 0)
            new_row1 = [istart_1, Int(round(entry[2])), i1 + istart_1 - 1, avl1, avr1, r1]

            istart_2 = Int(round(entry[3])) + 1
            istop_2  = Int(round(entry[2]))
            seg2 = dataX[istart_2:istop_2]
            i2, avl2, avr2, r2, _ = split_fast(seg2; demo=demo > 0)
            new_row2 = [istart_2, istop_2, i2 + istart_2 - 1, avl2, avr2, r2]

            splitlog_entry = [istop_1, new_row1[3], new_row2[3]]
            split_table = vcat(split_table[1:idx-1, :], reshape(new_row1, 1, :), reshape(new_row2, 1, :), split_table[idx+1:end, :])
            return split_table, splitlog_entry
        else
            return split_table, [0, 0, 0]
        end
    catch e
        @warn "Error in A1212_adapt_splitlog_table: $e"
        return split_table, [0, 0, 0]
    end
end

# ╔═╡ 993170ad-c13f-401e-8f92-291033bfdfd6
function A121_split_until_ready(dataX::Vector{<:Real}, N_iter::Int, demo::Float64)
    split_table = A1210_set_up_splitlogtable(dataX)

    # ✅ Defensive check — must have at least 3 rows (padding + valid row)
    if size(split_table, 1) < 3
        @warn "split_table has fewer than 3 rows; skipping step fitting"
        return zeros(N_iter), zeros(Int, N_iter, 3)
    end

    row2 = split_table[2, :]  # second row: actual initial segment info

    if length(row2) < 6
        @warn "split_table[2] has fewer than 6 elements; invalid format"
        return zeros(N_iter), zeros(Int, N_iter, 3)
    end

    i_start = Int(row2[1])
    i_stop  = Int(row2[2])
    i_next  = Int(row2[3])
    avl     = row2[4]
    avr     = row2[5]

    splitlog = zeros(Int, N_iter, 3)
    FitX = fill(mean(dataX), length(dataX))
    cFitX = zeros(length(FitX))

    # ✅ Bounds-safe assignment
    if i_start <= i_next < length(dataX)
        cFitX[i_start:i_next] .= avl
    end
    if i_next + 1 <= i_stop <= length(dataX)
        cFitX[i_next + 1:i_stop] .= avr
    end

    S_curve = zeros(N_iter)
    c = 0

    for ii in 1:N_iter
        segment_lengths = split_table[:, 2] .- split_table[:, 1] .+ 1
        rankings = split_table[:, 6]
        ix_valids = findall(x -> x > 2, segment_lengths) ∩ findall(x -> x > 0, rankings)

        if !isempty(ix_valids)
            subsel_idx = argmax(rankings[ix_valids])
            best_row_idx = ix_valids[subsel_idx]
            FitX = A1211_adapt_fit(FitX, split_table, best_row_idx)
            split_table, splitlog_entry = A1212_adapt_splitlog_table(dataX, split_table, best_row_idx, demo)
            cFitX = A1213_adapt_counterfit(cFitX, dataX, split_table, best_row_idx)
            splitlog[c+1, :] = splitlog_entry
            S_curve[c+1] = mean((dataX .- cFitX).^2) / mean((dataX .- FitX).^2)
            c += 1
        end
    end

    return S_curve, splitlog
end

# ╔═╡ 7520a659-e175-43d0-a8e0-02c436d3d09e
function append_fitX(newFitX::Vector, FitX::Vector, dataX::Vector)
    combiFitX = newFitX .+ FitX
    ixes0 = findall(diff(combiFitX) .!= 0)

    if !isempty(ixes0)
        ixes = [-1; ixes0; length(combiFitX)]
        Nmin = 2
        too_close = findall(diff(ixes) .< Nmin)

        for ix in too_close
            lo = ix - 1
            if lo + 4 > length(ixes)
                continue
            end
            ixlo = ixes[lo + 1]
            ixhi = ixes[lo + 4]
            segment = dataX[(ixlo + 1):(ixhi)]
            idx, avl, avr, _, _ = split_fast(segment)
            combiFitX[ixlo + 1:ixlo + idx] .= avl
            combiFitX[ixlo + idx + 1:ixhi] .= avr
        end
    end

    return combiFitX
end

# ╔═╡ c1ec08dc-d7e8-4a2a-92c1-24cdb3a550f9
function A122_eval_Scurve(S_curve::Vector{<:Real}, acceptance::Real=0.15)
    S1 = copy(S_curve)
    S1[S_curve .< 1.0] .= 1.0
    LS = length(S1)
    baseline = range(1, S1[end], length=LS)
    S_curve_fin = S1 .- baseline
    i_max = argmax(S_curve_fin)
    best_shot = S_curve_fin[i_max] > acceptance ? i_max : 0
    return best_shot, S_curve_fin
end

# ╔═╡ 7671fd09-13fb-482c-8e59-e9c04f6d1cb1
function Splitlog2FinalIndices(splitlog::Matrix{<:Integer}, best_shot::Int)
    indices_bestfit = sort(splitlog[1:best_shot, 1])
    counter_all = vcat(splitlog[1:best_shot, 2], splitlog[1:best_shot, 3])
    is_used = zeros(Int, length(counter_all))
    for (i, val) in enumerate(indices_bestfit)
        is_used[findall(x -> x == val, counter_all)] .= 1
    end
    indices_counterfit = sort(counter_all[is_used .== 0])
    return indices_bestfit, indices_counterfit
end

# ╔═╡ 4ec5ffa1-88f0-4bbf-a0ef-04f09678823b
function stepfindcore(dataX::Vector{<:Real}; demo=0.0, tresH=0.15, N_iter=100)
    N = length(dataX)
    N_iter = N_iter == 0 || N_iter > div(N, 4) ? div(N, 4) : N_iter

    S_curve, splitlog = A121_split_until_ready(dataX, N_iter, demo)
    best_shot, S_curve_fin = A122_eval_Scurve(S_curve, tresH)

    if best_shot > 0
        indices_bestfit, indices_counterfit = Splitlog2FinalIndices(splitlog, best_shot)
        FitX = indices_to_fit(dataX, indices_bestfit, how=:mean)
        cFitX = indices_to_fit(dataX, indices_counterfit, how=:mean)
    else
        FitX = zeros(length(dataX))
        cFitX = zeros(length(dataX))
        best_shot = 0
    end

    return FitX, cFitX, splitlog, S_curve_fin, best_shot
end

# ╔═╡ f4442f06-e6ea-45a6-bfbc-611a1e0f9698
function auto_stepfinder(data::Vector{Float64}; demo::Float64=0.0, N_iter::Int=100, tresH::Float64=0.15)
    FitX = zeros(length(data))
    Fits = Matrix{Float64}(undef, 0, length(data))
    S_curves = Vector{Vector{Float64}}()
    best_shots = Int[]

    for pass in 1:3
        residu = data .- FitX
        newFitX, _, _, S_curve, best_shot = stepfindcore(residu; demo=demo, tresH=tresH, N_iter=N_iter)
        FitX = append_fitX(newFitX, FitX, data)

        if best_shot > 0
            Fits = vcat(Fits, copy(FitX)')
            push!(S_curves, S_curve)
            push!(best_shots, best_shot)
        end
    end

    steptable = fit_to_steps(data, FitX)
    return FitX, Fits, S_curves, best_shots, steptable
end

# ╔═╡ 1f125790-241d-49e2-94e2-37881cb0517b
md"""

## Analysis
"""

# ╔═╡ e91084b3-c77d-49af-8e81-5370935ab3ce
data_path = "/Users/jjgomezcadenas/JJDropbox/2025/BOLD/jacobkers-AutoStepFinder-73f0028/test_data/Test_Data_AutoStepfinder/testdata_regular_steps.txt"

# ╔═╡ 8e3e2670-b7f4-4256-a9ac-ecd5fb9b7d19
data = readdlm(data_path, Float64)  # load as Float64 vector

# ╔═╡ cf7b955c-9fdf-4587-89a3-611d0eb0bc40
dataX = data[:, 2] 

# ╔═╡ 21bd743b-62d0-46e5-a897-6c3440181453
length(dataX)

# ╔═╡ 2b77ec77-f566-452d-8825-c138f7de35b3
fX, fs, sc, bs, st = auto_stepfinder(dataX; demo=0.0)

# ╔═╡ 952eb4be-4adc-48f2-acce-c98cb24d7e0a
data2X = generate_steps(5)  # 5 steps, default settings

# ╔═╡ 03a3fcb7-e770-489b-9030-9a093f8bf490
FitX, Fits, S_curves, best_shots, steptable = auto_stepfinder(data2X; N_iter=20, tresH=0.25)

# ╔═╡ 1c22f9fc-95d2-4c7d-aebe-dfcd91890da9
begin
	plot(data2X, label="Raw Data", lw=1.5)
	plot!(FitX, label="Final Fit", lw=2)
end

# ╔═╡ d0fa27c3-d4bf-44f9-a8da-cbef2a7b4ccc
begin
	step_positions = steptable[:, 1]
	step_heights = steptable[:, 4]
	bar(step_positions, step_heights, label="Steps", xlabel="Index", ylabel="Step Size")
end

# ╔═╡ f879d0fc-20d7-402c-820d-c871f73af126
step_positions

# ╔═╡ 8496aae3-6d9c-46d2-8b9b-5382a8556fc3
step_heights

# ╔═╡ Cell order:
# ╠═98eb0dbc-0da3-11f0-3773-ef1c04c57127
# ╠═54acbf3b-d144-49cc-a668-bb88566d6f85
# ╠═b571e2a9-5e93-4ac5-8af1-d8c66466336c
# ╠═fefc626d-3728-4973-a015-d59952007754
# ╠═e121bde7-2005-4c22-80b8-de679cb80edd
# ╠═993170ad-c13f-401e-8f92-291033bfdfd6
# ╠═1881e12a-bdd3-4afd-ac2f-b5a45eb3ce80
# ╠═58b3137b-7e1f-45d0-918e-93f10cb55dea
# ╠═33484d1b-58db-479c-956e-98037ed7be7e
# ╠═1ab697d4-de94-45bd-8127-6f94258a9ec7
# ╠═7520a659-e175-43d0-a8e0-02c436d3d09e
# ╠═d1b6be43-ff7c-4721-8657-4b0578829ad3
# ╠═da6d50b4-2d05-40bc-af0d-8e715112c6bd
# ╠═ae100887-9df6-41ce-8e7a-79032c8a3ddd
# ╠═a58fcf00-96a7-48d2-af21-ec316662d542
# ╠═4ec5ffa1-88f0-4bbf-a0ef-04f09678823b
# ╠═f4442f06-e6ea-45a6-bfbc-611a1e0f9698
# ╠═5a5fd414-a681-4d7d-9ec6-93b206146ed5
# ╠═c1ec08dc-d7e8-4a2a-92c1-24cdb3a550f9
# ╠═7671fd09-13fb-482c-8e59-e9c04f6d1cb1
# ╠═1f125790-241d-49e2-94e2-37881cb0517b
# ╠═e91084b3-c77d-49af-8e81-5370935ab3ce
# ╠═8e3e2670-b7f4-4256-a9ac-ecd5fb9b7d19
# ╠═cf7b955c-9fdf-4587-89a3-611d0eb0bc40
# ╠═21bd743b-62d0-46e5-a897-6c3440181453
# ╠═2b77ec77-f566-452d-8825-c138f7de35b3
# ╠═952eb4be-4adc-48f2-acce-c98cb24d7e0a
# ╠═03a3fcb7-e770-489b-9030-9a093f8bf490
# ╠═1c22f9fc-95d2-4c7d-aebe-dfcd91890da9
# ╠═d0fa27c3-d4bf-44f9-a8da-cbef2a7b4ccc
# ╠═f879d0fc-20d7-402c-820d-c871f73af126
# ╠═8496aae3-6d9c-46d2-8b9b-5382a8556fc3
