#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# cell calling / count cutoff / curve calculation util funcitons

import numpy as np
import pandas as pd
import scipy.interpolate
import matplotlib.pyplot as plt


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    pass
else:
    pass


### ------------------------------ Count curve cutoff stuff ---------------------------


def bc_tscp_count_df(tas_df, cell_bc=None, sort_index=False):
    """Get sorted high-to-low transcript count per cell (barcode) dataframe

    tas_df = transcript assignment dataframe
    cell_bc = optional cells to use; Dataframe with index = 'bc_wells' or list of bc
    sort_index = flag to sort by index rather than descending count

    Return dataframe
    """
    tscp_per_cell = tas_df.groupby("bc_wells").size().sort_values(ascending=False)
    bcc_df = pd.DataFrame(tscp_per_cell, columns=["tscp_count"])
    bcc_df.index.name = "bc_wells"

    # Given cell barcodes? Could already be dataframe or just list
    cellbc_df = None
    if cell_bc is not None:
        if isinstance(cell_bc, pd.DataFrame):
            cellbc_df = cell_bc
        else:
            # Make sure cell_bc is a list (Pandas doesn't allow a set)
            cellbc_df = pd.DataFrame(index=list(cell_bc))
            cellbc_df.index.name = "bc_wells"

    # Make dataframe with given-cell rows if given cell barcodes
    if cellbc_df is not None:
        bcc_df = cellbc_df.join(bcc_df)
        bcc_df["tscp_count"].fillna(0, inplace=True)
        bcc_df["tscp_count"] = bcc_df["tscp_count"].astype(int)

    # Sort
    # By (barcode) index
    if sort_index:
        bcc_df = bcc_df.sort_index()
    # By descending count, then assending bc (which is index)
    else:
        bcc_df = bcc_df.sort_values(["tscp_count", "bc_wells"], ascending=[False, True])

    return bcc_df


def load_frag_count_df(fname, bc1_seq=True, verb=True):
    """Load ATAC fragments file

    bc1_seq = flag to include round 1 barcode sequence in dataframe

    Return dataframe with (fragment) counts and rank
    """
    cols = "chrom,st,en,bc,count".split(",")
    frags_df = pd.read_csv(fname, sep="\t", header=None, names=cols)
    if verb:
        print(f"Loaded {fname}; {frags_df.shape}")

    # Counts of lines for each barcode (not 'count' of reads)
    counts = frags_df.groupby("bc").size().sort_values(ascending=False)
    frags = list(counts)
    ranks = range(1, len(counts) + 1)
    df = pd.DataFrame({"rank": ranks, "counts": frags}, index=counts.index)

    # Split out round 1 barcode?
    if bc1_seq:
        df["bc1_seq"] = df.index.str[:8]

    if verb:
        print(f"Count df with {len(df)} ranks")
    return df


def loglog_count_df(df, rcol="rank", ccol="counts", copy=False):
    """Log10 transform count and rank in dataframe

    rcol = column with ranks
    ccol = column with counts
    copy = flag to copy dataframe; else modify given

    Return dataframe
    """
    if copy:
        df = df.copy()
    df[ccol] = np.log10(df[ccol])
    df[rcol] = np.log10(df[rcol] + 1)
    return df


### ------------------------------- General 'curve' array utils ---------------------------


def count_vals_in_range(counts, min_v=0, max_v=0):
    """Get number of counts in given range

    counts = array or series (e.g. counts per cell, as rank curve)
    min_v = possible min value (e.g. to qualify cells)
    max_v = possible max value (e.g. to qualify cells)

    Return int
    """
    # Get to numpy array and then bound with any cutoffs
    counts = np.array(counts)
    if min_v and max_v:
        counts = (counts >= min_v) & (counts <= max_v)
    elif min_v:
        counts = counts >= min_v
    elif max_v:
        counts = counts >= max_v
    # Number qualified = sum of bool
    return int(counts.sum())


def mean_win_smooth(vals, win_frac):
    """Smooth values by rolling mean of window +/- from each point

    vals = array or series
    win_frac = data fraction for +/- window size

    Return list same length as input
    """
    vals = np.array(vals)
    # Bound window fraction
    win_frac = min(1, max(0, win_frac))
    win_side = round(len(vals) * win_frac)
    new_vals = []
    for i in range(len(vals)):
        first = max(0, i - win_side)
        last = min(i + win_side + 1, len(vals))
        v = sum(vals[first:last]) / (last - first)
        new_vals.append(v)
    return np.array(new_vals)


def first_derivative(vals, win_frac):
    """Get slope as delta vals +/- window for each point

    vals = array or series
    win_frac = data fraction for +/- window size

    Return list same length as input
    """
    vals = np.array(vals)
    # Bound window fraction
    win_frac = min(1, max(0, win_frac))
    win_side = round(len(vals) * win_frac)
    new_vals = []
    for i in range(len(vals)):
        first = max(0, i - win_side)
        last = min(i + win_side + 1, len(vals) - 1)
        v = vals[last] - vals[first]
        new_vals.append(v)
    return np.array(new_vals)


def index_past_val(vals, target, ascending=False, last=False):
    """Scan array for (first) value "past" target

    vals = array or series (i.e. sorted counts representing rank curve)
    target = target value to find (on curve)
    ascending = flag to scan "upward" past target value, vs "downward" below target
    last = flag to set value to the last value even if target not found

    Return index; -1 if not found
    """
    i1 = -1
    for i in range(len(vals)):
        if ascending:
            if vals[i] > target:
                i1 = i
                break
        else:
            if vals[i] < target:
                i1 = i
                break
    # If last flag and not found, set to last index
    if last and (i1 < 0):
        i1 = len(vals) - 1
    return i1


def indexes_with_val(vals, target, slop=0):
    """Scan array for (first) values within +/- slop of target

    vals = array or series (i.e. sorted points representing ranked count curve)
    target = target value
    slop = margin above / below target to qualify

    Return tuple of array indexes, first,last with target value
    """
    min_t = target - slop
    max_t = target + slop
    i1 = i2 = -1
    for i in range(len(vals)):
        # In target range
        if min_t <= vals[i] <= max_t:
            # first index
            if i1 == -1:
                i1 = i
        else:
            # Have legit first index, then done
            if i1 != -1:
                i2 = i
                break
    # If legit first but not last, set last
    if (i1 >= 0) and (i2 < 0):
        i2 = len(vals)
    return (i1, i2)


def get_auc_index(vals, frac=0.5):
    """Get the index corresponding to area under curve of given fraction

    vals = array or series (i.e. sorted points representing ranked count curve)

    Return int
    """
    frac = max(0, min(frac, 1))
    total = vals.sum()
    index = len(vals) - 1
    r_sum = 0
    for i, v in enumerate(vals):
        r_sum += v
        if (r_sum / total) >= frac:
            index = i
            break
    return index


def curve_norm_local_slope(vals, target, frac=0.01):
    """Get negative normalized slope at target on count curve

    vals = sampled and scaled count "curve" y component
    target = target value to find (on curve)
    frac = fraction of curve +/-win that is evaluated for slope

    Return tuple (slope, dy, dx)
    """
    slope = dx = dy = None
    # Find first index with fewer transcripts than target
    i1 = index_past_val(vals, target, ascending=False)
    # Need at least two points and the curve can't be flat
    if (i1 > 0) and (vals[0] > vals[-1]):
        # First and last indexes with value before "past" index
        i1, i2 = indexes_with_val(vals, vals[i1 - 1])

        # Delta up/down indexes
        ijump = round(len(vals) * frac / 2)
        i1 = max(0, i1 - ijump)
        i2 = min(len(vals) - 1, i2 + ijump)
        # Sample window with at least two sample points
        if i1 < i2:
            # Global slope components
            g_dy = vals[0] - vals[-1]
            g_dx = len(vals)
            # Window slope components, normalized by global deltas
            dy = (vals[i2] - vals[i1]) / g_dy
            dx = (i2 - i1) / g_dx
            slope = dy / dx * -1
    return (slope, dy, dx)


### ------------------------- Cutoff algorithm stuff ---------------------------


def umb_parab_tran_func(n_steps, xfrac=0.5, yfrac=0.25, plot=False):
    """Get umbrella-like parabolic transfer function centered in X, Y bounded 0-1

    xfrac = fraction of X scaled below 1;
        If xfrac=0, func(x) = 1; If xfrac = 1, func = (x-centered) 1-x^2
    yfrac = fraction of Y scaled at edges; If yfrac=0, func(x) = 1

    return numpy array len n_steps
    """
    # Number of steps; Count if given list and make sure at least 2
    if not isinstance(n_steps, int):
        n_steps = len(n_steps)
    n_steps = max(2, n_steps)
    # Bound fractions
    xfrac = min(1, max(0, xfrac))
    yfrac = min(1, max(0, yfrac))

    # Get X and first Y (= x-center parabola)
    x = np.array([i / (n_steps - 1) for i in range(n_steps)])
    y = (x * 2 - 1) ** 2
    # Shift down to keep fraction of X center flat (this part = no change)
    y_shift = (1 - xfrac) ** 2
    y2 = np.clip(y - y_shift, 0, 1)
    # Scale Y to 0-1
    if y_shift < 1:
        y2p = y2 / (1 - y_shift)
    else:
        y2p = y2
    # Invert and scale to Y fraction
    yfin = 1 - y2p * yfrac

    if plot:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111)
        ax.plot(x, y, color=(0.5, 0.8, 0.4))
        ax.plot(x, y2, color=(0.7, 0.7, 0.3))
        ax.plot(x, y2p, color=(0.9, 0.4, 0.2))
        ax.plot(x, yfin, color=(1, 0, 0), marker="o")
        ax.vlines(xfrac / 2, 0, 1, linestyle="dotted")
        ax.vlines(1 - xfrac / 2, 0, 1, linestyle="dotted")
        ax.hlines(1 - yfrac, 0, 1, linestyle="dotted")

    return yfin


def calc_count_cutoff(counts, recipe, plot=None, verb=True, return_recipe=True):
    """Calculate count cutoff

    counts = array or series of counts (i.e. rank curve); raw counts, not log scale
    recipe = dict of params (below)
    plot = flag or dict of params (below)
    verb = flag to report some story if problems encountered
    return_recipe = flag to include given recipe parameters with results

    recipe = {
        'min_cnt': 30,
        'use_cnt': 0,
        'min_cells': 10,
        'max_cells': 30000,
        'min_cell_auc': 0.4,
        'use_cells': 0,
        'cv_steps': 1000,
        'cv_smooth_win': 0.02,
        'cv_slope_win': 0.02,
        'edge_pen_x': 0.5,
        'edge_pen_y': 0.25,
        'cnt_scale_fac': 1,
    }

    If plot is given, plot will be generated.
        If plot is a dict, these params are used.
        Useful to set limits to 'zoom in' to regions of curve.
        Also can set plot title and an output filename (to save figure)
        plot = {
            'xlim': (2, 3),
            'ylim': (1.5, 4),
            'title': 'some title',
            'filename': 'filename.png',
        }

    Return tuple (status, dict)
    """
    # Make sure array with decreasing rank of counts
    counts = sorted(np.array(counts), reverse=True)

    # Return recipe along with results?
    if return_recipe:
        results = {k: v for k, v in recipe.items()}
    else:
        results = {}
    results["problem"] = ""

    # Hardcoded plotting stuff; If plot given (e.g. True), make sure dict
    win_off = 0.2
    if plot:
        if not isinstance(plot, dict):
            plot = {"plot": True}

    # Given specific count or cell numbers to use; Adjust bounds as needed
    if recipe["use_cnt"] > 0:
        recipe["min_cnt"] = min(recipe["min_cnt"], recipe["use_cnt"])
    elif recipe["use_cells"] > 0:
        recipe["use_cells"] = min(recipe["use_cells"], len(counts) - 1)
        recipe["min_cells"] = min(recipe["min_cells"], recipe["use_cells"])
        recipe["max_cells"] = max(recipe["max_cells"], recipe["use_cells"])

    # Unpack / prepare parameters; Bounds and targets to log
    min_cnt = np.log10(max(1, recipe["min_cnt"]))
    min_cells = np.log10(max(1, recipe["min_cells"]))
    max_cells = np.log10(recipe["max_cells"])
    min_cell_auc = recipe["min_cell_auc"]
    scale_fac = recipe["cnt_scale_fac"]
    cv_steps = recipe["cv_steps"]

    # Total barcodes (cells) and max X in log space
    nbc_tot = len(counts)

    # Create 'curve' in log log space (X = rank, Y = count)
    x = np.log10(np.arange(1, nbc_tot + 1))
    y = np.log10(counts)
    f = scipy.interpolate.interp1d(x, y, kind="linear")

    # Sampling positions for 'curve' and interpolated Y values
    x_curve = np.linspace(0, np.log10(nbc_tot), num=cv_steps)
    y_curve = f(x_curve)
    y_smooth = mean_win_smooth(y_curve, recipe["cv_smooth_win"])
    # Slope; Multiply by -1 to flip over
    y_slope = first_derivative(y_smooth, recipe["cv_slope_win"]) * -1

    ct_past = auc_past = 0
    # Where does curve cross min_cnt?; Set index past cells with at least min_cnt
    if y_smooth[0] > min_cnt:
        # If there are at least enough at first, all may qualify, so last=True
        ct_past = index_past_val(y_smooth, min_cnt, ascending=False, last=True)
    # Min fraction of tscp in cells?
    if min_cell_auc:
        auc_past = get_auc_index(y_smooth, frac=min_cell_auc)

    if plot:
        if "figsize" in plot:
            figsize = plot["figsize"]
        else:
            figsize = (12, 8)
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        ax.plot(x_curve, y_smooth, color=(0.9, 0.2, 0.1))
        ax.plot(x_curve, y_slope, color=(0.9, 0.2, 0.1))
        # No smoothing slope; Just for plotting
        y_slope_raw = first_derivative(y_curve, recipe["cv_slope_win"]) * -1
        ax.plot(x_curve, y_curve, color=(0.9, 0.6, 0.6))
        ax.plot(x_curve, y_slope_raw, color=(0.9, 0.6, 0.6))
        # Diagonal
        ax.plot(
            [0, x_curve[-1]],
            [y_curve[0], y_curve[-1]],
            color=(0.9, 0.6, 0.6),
            linestyle="dotted",
        )

        # Min and max cells
        ax.vlines(
            min_cells, 0, max(y_smooth), linestyle="dotted", color=(0.3, 0.3, 0.9)
        )
        ax.vlines(
            max_cells, 0, max(y_smooth), linestyle="dotted", color=(0.3, 0.3, 0.9)
        )
        ax.hlines(min_cnt, 0, x_curve[-1], linestyle="dotted", color=(0.6, 0.3, 0.9))
        if ct_past > 0:
            ax.vlines(
                x_curve[ct_past],
                0,
                max(y_smooth),
                linestyle="dotted",
                color=(0.6, 0.3, 0.9),
            )
        if auc_past > 0:
            ax.vlines(
                x_curve[:auc_past],
                0,
                y_smooth[:auc_past],
                alpha=0.05,
                color=(0, 0.8, 0),
            )

    # Get X-curve indexes past min and max cells; For max, can include all, so last=True
    cu_past = index_past_val(x_curve, min_cells, ascending=True)
    cu_past = max(cu_past, auc_past)
    cd_past = index_past_val(x_curve, max_cells, ascending=True, last=True)

    # Window on curve for threshold search
    # Start:end shifted toward more cells at this point; Can limit later
    win_st = max(0, cu_past)
    win_end = min(ct_past, cd_past) + 1
    # Window still must fit in overall curve
    win_end = min(win_end, len(x_curve) - 1)
    # Make sure window of at least two points
    win_span = win_end - win_st
    if win_span < 2:
        story = f"calc_count_cutoff: Bad sample window: start {win_st} end {win_end}"
        story += f"input counts: {len(counts)}"
        story += f"input recipe: {recipe}"
        if verb:
            print(story)
        results["problem"] = story
        return False, results

    # Edge penalty transfer function to flatten window ends
    # First arg just used for length
    n_steps = len(x_curve[win_st:win_end])
    y_tran = umb_parab_tran_func(
        n_steps, xfrac=recipe["edge_pen_x"], yfrac=recipe["edge_pen_y"]
    )
    y_tran = mean_win_smooth(y_tran, recipe["cv_smooth_win"])
    y_tslope = y_slope[win_st:win_end] * y_tran

    # Use or calculate cutoff
    if recipe["use_cnt"] > 0:
        # Given count thresh to use
        thresh_raw = threshold = recipe["use_cnt"]
    elif recipe["use_cells"] > 0:
        # Given cell number; use count at that (1-based) rank cell
        thresh_raw = threshold = counts[recipe["use_cells"] - 1]
    else:
        # Find cutoff de novo
        # Get index for peak in (possibly) tranformed slope curve
        max_y = max_i = 0
        for i in range(len(y_tslope)):
            if y_tslope[i] > max_y:
                max_y = y_tslope[i]
                max_i = i
        # Thresh = index in smoothed curve
        i_raw = max_i + win_st
        thresh_raw = y_smooth[i_raw]
        # Final count threshold = log10 to normal, then scaled
        thresh_raw = 10**thresh_raw
        threshold = thresh_raw * scale_fac

    # Final index into curve
    thresh_final = np.log10(threshold) + 1 / cv_steps
    i_fin = index_past_val(y_smooth, thresh_final, ascending=False, last=True)

    results["threshold"] = threshold
    results["thresh_raw"] = thresh_raw
    results["nbc_tot"] = nbc_tot
    results["thwin_min"] = 10 ** x_curve[win_st]
    results["thwin_max"] = 10 ** x_curve[win_end]
    results["thpoint_frac_cell"] = i_fin / cv_steps
    results["thpoint_frac_cnt"] = thresh_final / max(y_curve)
    # Cutoff point slope, different size windows
    slope, dy, dx = curve_norm_local_slope(y_curve, thresh_final, frac=0.01)
    results["thpoint_f01_slope"] = slope
    results["thpoint_f01_dy"] = dy
    results["thpoint_f01_dx"] = dx
    slope, dy, dx = curve_norm_local_slope(y_curve, thresh_final, frac=0.03)
    results["thpoint_f03_slope"] = slope
    results["thpoint_f03_dy"] = dy
    results["thpoint_f03_dx"] = dx
    slope, dy, dx = curve_norm_local_slope(y_curve, thresh_final, frac=0.1)
    results["thpoint_f10_slope"] = slope
    results["thpoint_f10_dy"] = dy
    results["thpoint_f10_dx"] = dx
    # Cell counts
    results["ncell_thresh"] = count_vals_in_range(counts, min_v=results["threshold"])
    results["ncell_raw"] = count_vals_in_range(counts, min_v=results["thresh_raw"])

    if plot:
        # Plot transfer function and scaled slope; with vertical offset
        ax.plot(
            x_curve[win_st:win_end],
            y_slope[win_st:win_end] + win_off * 1,
            color=(0.0, 0.7, 0.9),
            linestyle="dashed",
        )
        ax.plot(x_curve[win_st:win_end], y_tslope + win_off * 2, color=(0.0, 0.7, 0.9))
        ax.plot(x_curve[win_st:win_end], y_tran, color="purple", linestyle="dashed")

        story = f"Raw {round(results['thresh_raw'])} = {results['ncell_raw']} cells"
        # print(story)
        thresh = np.log10(results["thresh_raw"])
        ax.hlines(
            thresh,
            0,
            max(x_curve),
            linestyles="dashed",
            color=(0.0, 0.7, 0.9),
            label=story,
        )
        cells = np.log10(results["ncell_thresh"])
        ax.vlines(cells, 0, max(y_curve), linestyles="dashed", color=(0.0, 0.7, 0.9))

        story = f"Final {round(results['threshold'])} = {results['ncell_thresh']} cells"
        # print(story)
        thresh = np.log10(results["threshold"])
        ax.hlines(
            thresh,
            0,
            max(x_curve),
            linestyles="dashed",
            color=(0.1, 0.1, 0.9),
            label=story,
        )
        ax.vlines(cells, 0, max(y_curve), linestyles="dashed", color=(0.1, 0.1, 0.9))

        # If given explicit values, over color
        if recipe["use_cnt"] > 0:
            story = f"Given count {round(recipe['use_cnt'])}"
            ax.hlines(
                thresh,
                0,
                max(x_curve),
                linewidth=6,
                color=(0.8, 0.9, 0, 0.5),
                label=story,
            )
        if recipe["use_cells"] > 0:
            story = f"Given cells {round(recipe['use_cells'])}"
            ax.vlines(
                cells,
                0,
                max(y_curve),
                linewidth=6,
                color=(0.8, 0.9, 0, 0.5),
                label=story,
            )

        ax.legend()
        if "xlim" in plot:
            plt.xlim(plot["xlim"])
        if "ylim" in plot:
            plt.ylim(plot["ylim"])

        ax.set_xlabel("log10( Rank )", fontsize=15)
        ax.set_ylabel("log10( Count )", fontsize=15)
        if "title" in plot:
            ax.set_title(plot["title"], fontsize=20)

        if "filename" in plot:
            fig.savefig(plot["filename"])
            print("New file:", plot["filename"])

        plt.show()

    return True, results
