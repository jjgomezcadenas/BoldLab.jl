import numpy as np


def stepfindcore(dataX=1,tresH=0.05, N_iter=0, demo=0):

    # run an iterative step fit deep into overfitting:
    if N_iter == 0 or N_iter > int(len(dataX) / 4):
        N_iter = int(len(dataX) / 4)
        
    S_curve, splitlog = A121_split_until_ready(dataX, N_iter, demo)

    # estimate best step number in this iteration range:
    best_shot, S_curve_fin = A122_eval_Scurve(S_curve, tresH)


    # re-build best fit from 'split log':
    if best_shot > 0:
        indices_bestfit, indices_counterfit = Splitlog2FinalIndices(splitlog, best_shot)
        FitX = Indices2Fit(dataX, indices_bestfit, how="mean")
        cFitX = Indices2Fit(dataX, indices_counterfit, how="mean")

    else:
        FitX = np.zeros_like(dataX)
        cFitX = np.zeros_like(dataX)
        best_shot = 0

    return FitX, cFitX, splitlog, S_curve_fin, best_shot


def A121_split_until_ready(dataX=1, N_iter=50, demo=0):
    """
    iterative splitting of ever-smaller data segments, adding one step location add a time.
    For each added step location, the new fit and new counterfit are calculated.
    %This fit pair -fit and counterfit- are then used to calculate the S-value of this fit pair
    For efficient recording of the iteration sequence, all basic information
    (step locations, step level averages left and right of the step etc.)
    is stored in an ever-expanding 'split_table'.
    The final information that is passed is the 'splitlog, a single-column list of the indices
     to which each new stepfit was added, in order of use.
    @author: jkerssemakers feb 23 2022"""

    # set up
    split_table = A1210_set_up_splitlogtable(dataX)
    splitlog = np.zeros((N_iter, 3), dtype=np.uint32)

    # build initial fit:
    FitX = np.full_like(dataX, np.mean(dataX))

    # build initial counterfit:
    cFitX = np.zeros_like(FitX)
    i_start = int(split_table[1, 0])
    i_stop = int(split_table[1, 1])
    i_next = int(split_table[1, 2])
    avl = split_table[1, 3]
    avr = split_table[1, 4]
    cFitX[i_start : i_next + 1] = avl
    cFitX[i_next + 1 : i_stop + 1] = avr

    # build initial S-curve
    S_curve = np.zeros(N_iter)
    c = 0
    for ii in range(1, N_iter + 1, 1):
        # for row entries of valid-to-fit segments:
        # check validity
        # properties from split table:
        segment_lengths = np.ndarray.flatten(split_table[:, 1] - split_table[:, 0]) + 1
        rankings = np.ndarray.flatten(split_table[:, 5])

        # valid indices (expand the condition):
        ix_valids = np.argwhere((segment_lengths > 2) & (rankings > 0))

        if ix_valids.size > 0:
            # if any left, find the one best to split next:
            # if not, stop
            subsel_idx = np.argmax(rankings[ix_valids])

            best_row_idx = int(ix_valids[subsel_idx])
            best_row_entry = np.ndarray.flatten(split_table[best_row_idx, :])

            # adapt the fit:
            FitX = A1211_adapt_fit(FitX, split_table, best_row_idx)

            # replace this row by two new ones; this includes two new stepfits:
            split_table, splitlog_entry = A1212_adapt_splitlog_table(
                dataX, split_table, best_row_idx, demo
            )
            # only after this, adapt the counterfit:
            cFitX = A1213_adapt_counterfit(cFitX, dataX, split_table, best_row_idx)

            # keep track of step placement order of counterfit
            # remove the entry equal to int(best_row_entry[2])
            # add the two new entries
            # keep track of step placement order of fit
            splitlog[c, 0:3] = splitlog_entry

            # adapt S-curve;
            # #includes correction for different degrees of freedom Fit & CFit
            corF = c / (c + 1)
            corF = 1
            S_curve[c] = (
                corF * np.mean((dataX - cFitX) ** 2) / np.mean((dataX - FitX) ** 2)
            )
            c = c + 1

    return S_curve, splitlog


def A1210_set_up_splitlogtable(dataX):
    """# set  up split log prior to iteration:
    first and last rows are dummy fillers
    ''new row'' contains entries corresponding to whole data trace"""
    Na = len(dataX)
    i_nxt, avl, avr, rankit, errorcurve = splitFast(dataX)
    new_row1 = [0, Na, i_nxt, avl, avr, rankit]
    split_table = np.array([[-1, -1, -1, 0, 0, 0], new_row1, [Na, Na, Na, 0, 0, 0]])
    return split_table


def A1211_adapt_fit(FitX, split_table, idx):
    """adapt the fit trace using split table"""
    # is it valid to split?
    rank = split_table[idx, 5]
    ok2adapt = (rank) > 0
    if ok2adapt:
        # indices:
        leftplateau_istart = int(split_table[idx, 0])
        leftplateau_istop = int(split_table[idx, 2])
        rightplateau_istart = int(leftplateau_istop) + 1
        rightplateau_istop = int(split_table[idx, 1])
        # plateau levels:
        leftplateau_level = split_table[idx, 3]
        rightplateau_level = split_table[idx, 4]
        FitX[leftplateau_istart:leftplateau_istop] = leftplateau_level
        FitX[rightplateau_istart:rightplateau_istop] = rightplateau_level
    return FitX


def A1213_adapt_counterfit(cFitX, dataX, split_table, idx):
    """adapt the counterfit trace for three new plateaus using four adjacent split table entries"""
    # check if middle two segments have valid entries:
    rankings = split_table[idx : idx + 2 :, 5]
    ok2adapt = all(rankings)
    if ok2adapt:
        # indices for three segments:
        leftplateau_istart = int(split_table[idx - 1, 2]) + 1
        leftplateau_istop = int(split_table[idx, 2])
        centerplateau_istart = leftplateau_istop + 1
        centerplateau_istop = int(split_table[idx + 1, 2])
        rightplateau_istart = centerplateau_istop + 1
        rightplateau_istop = int(split_table[idx + 2, 2])
        # plateau levels:
        leftplateau_level = np.mean(dataX[leftplateau_istart:leftplateau_istop])
        centerplateau_level = np.mean(dataX[centerplateau_istart:centerplateau_istop])
        rightplateau_level = np.mean(dataX[rightplateau_istart:rightplateau_istop])
        cFitX[leftplateau_istart:leftplateau_istop] = leftplateau_level
        cFitX[centerplateau_istart:centerplateau_istop] = centerplateau_level
        cFitX[rightplateau_istart:rightplateau_istop] = rightplateau_level

    return cFitX


def A1212_adapt_splitlog_table(segm=1, split_table=1, idx=1, demo=0):
    """Iterative update of a table describing a split sequence.
    This def adapts the splitlog table for a newly added step location:
    two new plateau-property rows replace the location of one old one
    Each row of this 'split table  describes a data segment:
    col1: 'start' index of segment
    col2: 'stop' index of segment
    if this segment would be split, it would be according to:
    col3: 'next' index:  split here
    col4: 'avl': left average after split
    col5: 'avr': right average after split
    col6: 'rankit' : the relative prominence step(see Split_fast)
    These last three columns are found by performing a 'proof-split' on this segment
    @author: jkerssemakers march 2 2022
    """
    # select the one row to adapt:
    best_row_entry = np.ndarray.flatten(split_table[idx, :])

    # construct new row 1: describing LEFT plateau:
    istart_1 = int(best_row_entry[0])
    istop_1 = int(best_row_entry[2])
    inxt_1, avl_1, avr_1, rankit_1, errorcurve_1 = splitFast(
        segm[istart_1 : istop_1 + 1], demo
    )
    new_row1 = [istart_1, istop_1, inxt_1 + istart_1, avl_1, avr_1, rankit_1]

    # construct new row 2: describing RIGHT plateau
    istart_2 = int(best_row_entry[2] + 1)
    istop_2 = int(best_row_entry[1])
    inxt_2, avl_2, avr_2, rankit_2, errorcurve_2 = splitFast(
        segm[istart_2 : istop_2 + 1], demo
    )
    new_row2 = [istart_2, istop_2, inxt_2 + istart_2, avl_2, avr_2, rankit_2]

    splitlog_entry = [istop_1, inxt_1 + istart_1, inxt_2 + istart_2]

    # replace fomer row entry by the two new row entries
    L_slt = len(split_table[:, 0])
    block_before = split_table[0:idx, :]
    block_after = split_table[idx + 1 : L_slt + 1]
    split_table = np.vstack(
        [
            block_before,
            np.array(new_row1),
            np.array(new_row2),
            block_after,
        ]
    )
    return split_table, splitlog_entry


def A122_eval_Scurve(S_curve, acceptance=0.15):
    """find best number of steps based on S-curve"""
    # shave off low part
    ix_low = np.ndarray.flatten(np.argwhere(S_curve < 1))
    S1 = S_curve.copy()
    S1[ix_low] = 1
    # remove baseline
    LS = len(S_curve)
    baseline = np.linspace(1, S_curve[LS - 1], LS)
    S_curve_fin = S1 - baseline
    # estimate best number of steps
    i_max = np.argmax(S_curve_fin)
    if S_curve_fin[i_max] > acceptance:
        best_shot = i_max + 1
    else:
        best_shot = 0
    return best_shot, S_curve_fin


def Splitlog2FinalIndices(splitlog, best_shot):
    """build indice where steps should be placed for the best fit and its associated counterfit"""
    indices_bestfit = sorted(splitlog[0:best_shot, 0])
    # collect ALL indices that were at least once part of a counter-fit
    indices_counterfit_all = np.ndarray.flatten(splitlog[0:best_shot, 1:3])
    # select all counterfit indices NOT existing in the best-fit collection;
    # these define the associated counter-fit
    is_it_used = np.zeros_like(indices_counterfit_all, dtype=int)
    LF = len(indices_bestfit)
    for ii in range(0, LF, 1):
        thisindex = indices_bestfit[ii]
        ix = np.argwhere(indices_counterfit_all == thisindex)
        is_it_used[ix] = 1

    indices_counterfit = sorted(
        np.ndarray.flatten(indices_counterfit_all[np.argwhere(is_it_used == 0)])
    )
    return indices_bestfit, indices_counterfit


def splitFast(segment, demo=0):
    """
    Determines the best step-fit in a one-dim array 'segment'.
    To save time, use of function 'mean' is avoided in the deep loop
    If the result is invalid, this is passed via 'rank' value=0
    Jacob Kers 2021
    """
    invalidFit = 0
    Nmin = 2
    Ns = len(segment)
    if Ns > 2:
        var_q = np.ones(Ns)
        # set up averaging buffers
        avl = segment[0]
        avr = np.sum(segment[1 : Ns + 1]) / (Ns - 1)
        ava = np.sum(segment) / Ns
        for ii in range(1, Ns - 1, 1):
            # current plateau lengths:
            n_L = ii
            n_R = Ns - ii
            # update left average:
            avl = (avl * (n_L - 1) + segment[ii]) / (n_L)
            # update right average:
            avr = (avr * (n_R + 1) - segment[ii]) / (n_R)
            # relative change in averages left and right
            delta_l = avl - ava
            delta_r = avr - ava
            # total variance correction (NEUTRALIZED):
            varcor_L = 1 + 0 * (n_L + 1) / (n_L)
            varcor_R = 1 + 0 * (n_R + 1) / (n_R)
            # total variance:
            delta_q = delta_l**2 * n_L * varcor_L + delta_r**2 * (n_R) * varcor_R
            var_q[ii] = -delta_q
        # wrapping up:
        idx = np.argmin(var_q)
        if (idx < Nmin - 1) | (Ns - idx < Nmin - 1):
            invalidFit = 1
        else:
            avl_fin = np.mean(segment[0:idx])
            avr_fin = np.mean(segment[idx + 1 : Ns + 1])
            rankit = (avr_fin - avl_fin) ** 2 * Ns
            errorcurve = var_q / Ns
            
    else:
        invalidFit = 1
    if invalidFit:
        idx = 0
        avl_fin = segment[0]
        avr_fin = segment[0]
        rankit = 0
        errorcurve = segment * 0
    return idx, avl_fin, avr_fin, rankit, errorcurve


def Indices2Fit(dataX, indices, how="mean"):
    """ "This function builds a step fit"""
    ixlo = 0
    LX = len(dataX)
    FitX = np.zeros_like(dataX)
    indices_ext = np.append(-1, indices)
    indices_ext = np.append(indices_ext, LX)
    Lix = len(indices_ext)
    for ii in range(0, Lix - 1, 1):
        ixlo = indices_ext[ii] + 1
        ixhi = indices_ext[ii + 1] + 1
        if ixhi >= ixlo:
            if how == "mean":
                FitX[ixlo:ixhi] = np.mean(dataX[ixlo:ixhi])
            elif how == "median":
                FitX[ixlo:ixhi] = np.median(dataX[ixlo:ixhi])
        else:
            FitX[ixlo:ixhi] = dataX[ixlo:ixhi]

    return FitX
