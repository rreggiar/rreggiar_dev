#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import re
import numpy as np
import pandas as pd

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import pb_const as const
    import utils
else:
    from splitpipe import pb_const as const
    from splitpipe import utils


def get_stat_dict_and_df(stat_dict, stat_df, col="value", to_dict=False, to_df=False):
    """Get both dict and dataframe from given input (either one)

    stat_dict = dict (with stats) or None
    stat_df = dataframe (with stats in col) or None

    Return tuple dict, df
    """
    new_dict = new_df = None
    # If neither, both
    if (not to_dict) and (not to_df):
        to_dict = to_df = True

    if isinstance(stat_df, pd.DataFrame) and to_dict:
        new_dict = stat_df[col].to_dict()
    if isinstance(stat_dict, dict) and to_df:
        new_df = pd.DataFrame.from_dict(stat_dict, orient="index", columns=[col])

    return new_dict, new_df


# ----------------- analysis summary keys ----------------------------------------

ASUM_KEY_DEFS = r"""
    # Key regex in order; More specific should be first
    # Species prefixes are genome names. These may have dash so '[\w-]' not just '\w' at start

    ^sample_well_count$
    ^(targeted_)?number_of_cells$       # Total cells; May be prefixed with targeted
    ^[\w-]*number_of_cells$             # Cells per species; Prefix genome name can have dash
    ^[\w-]*median_tscp_per_cell$
    ^[\w-]*median_tscp_at50$
    ^[\w-]*median_genes_per_cell$
    ^mean_reads_per_cell$
    ^number_of_reads$
    ^number_of_tscp$

    ^expected_number_guides$
    ^detected_number_guides$
    ^frac_cells_with_\d+_crispr_guide$     # Focal barcoding genes in cells
    ^frac_cells_with_\d+plus_crispr_guide$
    ^num_cells_with_\d+_crispr_guide$
    ^num_cells_with_\d+plus_crispr_guide$

    ^target_number$                     # Target enrichment stats
    ^targeted_m?read_fraction$
    ^targeted_tscp_fraction$
    ^targeted_gene_fraction$

    ^[\w-]*number_of_tscp               # Per species; Prefix can have dash
    ^(targeted_)?sequencing_saturation$
    ^valid_barcode_fraction$
    ^\w*_Q30$                           # Prefix is simple keyword
    ^transcriptome_map_fraction$
    ^tso_fraction_in_read1$
    ^[\w-]*fraction_m?reads_in_cells    # Per species; Prefix can have dash
    ^[\w-]*fraction_tscp_in_cells
    ^[\w-]*fraction_exonic
    ^cell_tscp_cutoff$
    ^cell_tscp_f01_slope$
    ^mapto_[\w-]*                 # Suffix is whatever seq is names; Can have dash

    ^tcr_
    ^guide_tscp_cutoff$           #CRISPR tscp cutoff
    ^guide_read_cutoff$           #CRISPR read cutoff
"""

ASUM_KEY_DEFS_CRISPR = r"""
    # Key regex in order; More specific should be first
    # Species prefixes are genome names. These may have dash so '[\w-]' not just '\w' at start

    ^sample_well_count$
    ^(targeted_)?number_of_cells$       # Total cells; May be prefixed with targeted
    ^(targeted_)?sequencing_saturation$
    ^[\w-]*number_of_cells$             # Cells per species; Prefix genome name can have dash
    ^[\w-]*median_tscp_per_cell$
    ^[\w-]*median_tscp_at50$
    ^[\w-]*median_genes_per_cell$
    ^mean_reads_per_cell$
    ^number_of_reads$
    ^number_of_tscp$

    ^expected_number_guides$
    ^detected_number_guides$
    #^frac_cells_with_\d+_focal_gene$     # Focal barcoding genes in cells
    #^frac_cells_with_\d+plus_focal_gene$
    #^num_cells_with_\d+_focal_gene$
    #^num_cells_with_\d+plus_focal_gene$

    ^transcriptome_map_fraction$
    ^valid_barcode_fraction$
    ^\w*_Q30$                           # Prefix is simple keyword
    ^cell_tscp_cutoff$
    ^guide_tscp_cutoff$                 #AK trying to change cell_tscp_cutoff for crispr
"""


def get_asum_out_keys(in_keys, is_crispr=False, targeted=False):
    """Get ordered subset of given keys for analysis summary files

    in_keys = list of keys to screen

    Return tuple(list of ordered keys, set of discarded keys)
    """
    key_defs = ASUM_KEY_DEFS
    if is_crispr:
        key_defs = ASUM_KEY_DEFS_CRISPR

    # Collect in-order list of keys, values
    key_list = []
    for line in key_defs.split("\n"):
        parts = line.split("#")[0].strip().split()
        if len(parts) < 1:
            continue

        regex = parts[0]
        # Collect all matching keys first, so can sort
        sm_keys = []
        for k in in_keys:
            # Skip if already have
            if k in key_list:
                continue
            if re.match(regex, k):
                sm_keys.append(k)

        for k in sorted(sm_keys):
            key_list.append(k)

    # Second pass for specific discards
    keep_list = []
    for k in key_list:
        # Ignore
        if "raw_" in k:
            continue
        # Targeted; Ignore non-targeted versions
        if targeted:
            if "fraction_exonic" in k:
                if "targeted" not in k:
                    continue
            if "sequencing_saturation" in k:
                if "targeted" not in k:
                    continue
        keep_list.append(k)

    discards = set(in_keys) - set(keep_list)

    return keep_list, discards


def get_ana_sum_dict_and_df(stat_dict, stat_df, col="value"):
    """Get both dict and dataframe from given input (either one)

    Return tuple dict, df
    """
    # Get both dict and dataframe (for limited and all value outputs)
    if isinstance(stat_df, pd.DataFrame):
        stat_dict = stat_df[col].to_dict()
    elif isinstance(stat_dict, dict):
        stat_df = pd.DataFrame.from_dict(stat_dict, orient="index", columns=[col])
    else:
        raise ValueError(
            f"Need either dict or dataframe; Got {type(stat_dict)} {type(stat_df)}"
        )
    return stat_dict, stat_df


def get_clean_ana_sum_df(
    spipe, stat_data, col="value", is_crispr=False, targeted=False
):
    """Clean up analysis summary dataframe from (full) stat input

    stat_data = dict or dataframe with full stats
    is_crispr = flag to format for crispr / focal

    Return dataframe
    """
    stat_dict, stat_df = get_ana_sum_dict_and_df(stat_data, stat_data, col=col)
    # Always have number_of_cells
    if "number_of_cells" not in stat_df.index:
        stat_df.loc["number_of_cells"] = 0
    stat_dict, stat_df = get_ana_sum_dict_and_df(stat_dict, stat_df, col=col)

    keep_list, discards = get_asum_out_keys(
        stat_dict, is_crispr=is_crispr, targeted=targeted
    )
    keep_df = stat_df.loc[keep_list]

    new_index = []
    for key in keep_df.index:
        if is_crispr:
            key = key.replace("focal_gene", "guide")
            key = key.replace("transcriptome", "guide")

        key = key.replace("_mreads_", "_reads_")
        key = key.replace("_mread_", "_read_")
        key = key.replace("_focal_gene_", "_guide_")

        new_index.append(key)

    keep_df.index = new_index
    return keep_df


def save_analysis_summary(spipe, samp, stat_data, col="value", round_to=4):
    """Save all stats and analysis summary files for sample

    stat_data = stat dictionary or dataframe; Should have all stats (not just summary)
    col = dataframe column to use / create
    round_to = limit on precision

    Return nothing
    """
    targeted = spipe.is_targeted()

    # Full set; all values saved
    stat_dict, stat_df = get_ana_sum_dict_and_df(stat_data, stat_data, col=col)
    stat_df.index.name = "statistic"
    spipe.write_df(stat_df, "SFR_ALLSTATS", samp=samp)

    # Get dataframe with subset / reformatted keys for summary
    stat_df = get_clean_ana_sum_df(
        spipe, stat_data, col=col, is_crispr=False, targeted=targeted
    )

    # Collect values, clean up
    k_list = []
    v_list = []
    for k in stat_df.index:
        v = stat_df.loc[k, col]
        # Limit precision if number
        ok, v = utils.str_to_num(v, as_float=True)
        if ok:
            v = round(v, round_to)
        k_list.append(k)
        v_list.append(v)

    # Back to dataframe for saving
    for k, v in zip(k_list, v_list):
        stat_df.loc[k, col] = v

    # Save only (ordered) key subset of rows
    stat_df = stat_df.loc[k_list, :]
    stat_df.index.name = "statistic"
    spipe.write_df(stat_df, "SFR_ASUM_CSV", samp=samp, sort_index=False)


# ----------------- Rules for combing (DGE) stats

STAT_SUM_KEY_DEFS = r"""
    # Regular expression to match stats that should be sum
    #   Default is to combine via weighted mean
    # Third column is optional condition; 'combine' or 'single'
    # Order matters; Most specific patterns should be first

    # Keys starting with 'read_' are from preprocessing or align
    #   These only get summed for different sublibs in combine mode
    ^reads_.*$              sum     combine
    # Raw numbers only get summed for different sublibs
    ^raw_.*$                sum     combine

    # Sample well count mean in combine mode = mean, single = sum
    ^sample_well_count$     mean    combine
    ^sample_well_count$     sum     single

    # Number of cells always sum
    ^.*number_of_cells$     sum
    ^.*number_of.*$         sum
    ^.*num_cells_with.*$    sum
    ^.*count.*$             sum
    ^bc_edit_dist.*$        sum
    ^sample_tscp_fraction$  sum
"""


def get_stat_sum_mean_keys(in_keys, combine=True, def_op="mean"):
    """Get subsets of given keys that get combined via sum or mean

    in_keys = list of keys to screen

    Return tuple(list of sum keys, set of mean keys)
    """
    # Collect regular expressions for explicit sum or mean collection
    reg_list = []
    for line in STAT_SUM_KEY_DEFS.split("\n"):
        parts = line.split("#")[0].strip().split()
        if len(parts) < 2:
            continue
        logic = "" if len(parts) < 3 else parts[2]
        reg_list.append([parts[0], parts[1], logic])

    # Save keys for either sum or mean processing
    sum_list = []
    mean_list = []
    for k in in_keys:
        use_op = def_op
        # Screen reg expression collection to find any match
        for regex, operator, logic in reg_list:
            if re.match(regex, k):
                # Qualified logic for combine or single?
                # Only set operator if combine / not-combine is appropriate
                if logic == "combine":
                    if combine:
                        use_op = operator
                elif logic == "single":
                    if not combine:
                        use_op = operator
                # Otherwise just use operator
                else:
                    use_op = operator
                break

        if use_op == "sum":
            sum_list.append(k)
        elif use_op == "mean":
            mean_list.append(k)

    return sum_list, mean_list


def aggregate_asum_csv(spipe):
    """Aggregate analysis_summary csv files for samples together

    Returns nothing
    """
    spipe.report_run_story("Aggregating sample analysis_summary csv files")

    df_lis = []
    col_lis = []
    for samp in spipe.get_samples():
        # Make sure file was actually written
        fpath = spipe.filepath("SFR_ASUM_CSV", samp)
        if utils.check_infile(fpath):
            df = spipe.read_csv("SFR_ASUM_CSV", samp)
            df_lis.append(df.iloc[:, 0])
            col_lis.append(samp.get_name())

    if df_lis:
        # Concat column-wise and write
        ag_df = pd.concat(df_lis, axis=1)
        ag_df.columns = col_lis

        # Clean up extraneous
        keep_rows = []
        for row in ag_df.index:
            if "_mreads_" in row:
                continue
            keep_rows.append(row)
        ag_df = ag_df.loc[keep_rows, :]

        ag_df.index.name = "statistic"
        spipe.write_df(ag_df, "SF_AGG_ASUM_CSV", sort_index=False)


def updated_stats_dict(spipe, samp, stat_dict, targeted=False):
    """Update (initially combined) stats by recalculating some

    stat_df = (merged) all_stats dataframe
    targeted = flag for targeted keywords or not

    Return dict of updated values
    """
    # Extra part of keywords if targeted
    xkey = "targeted_" if targeted else ""

    # Saturation from saved counts (per sample)
    fkey = "SFR_SAMP_MRTG_CT"
    rtg_df = spipe.check_read_csv(fkey, samp=samp)
    if rtg_df is None:
        story = f"No MRGT count file for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)
        return stat_dict

    rtg_dict = rtg_df["count"].to_dict()
    stat_dict[f"{xkey}sequencing_saturation"] = (
        1 - rtg_dict[f"{xkey}total_tscp_count"] / rtg_dict[f"{xkey}total_mread_count"]
    )

    # Read estimate and sample fraction; Need raw counts in stats_dict
    num_tscp = stat_dict["number_of_tscp"]
    num_targ_tscp = (
        stat_dict["targeted_number_of_tscp"]
        if "targeted_number_of_tscp" in stat_dict
        else 0
    )

    stats = get_samp_well_and_tscp_info(spipe, samp, stat_dict, num_tscp, num_targ_tscp)
    stat_dict.update(stats)
    stats = get_samp_reads_info(spipe, samp, stat_dict)
    stat_dict.update(stats)

    dge_fkeys = const.out_filepaths(spipe, samp, targeted=targeted, as_key=True)
    cell_df = spipe.check_read_csv(dge_fkeys["cell"], samp=samp)
    if cell_df is None:
        story = f"No cell file for {samp.get_name()} ({dge_fkeys['cell']})"
        spipe.report_warning(story)
        return stat_dict

    # Calc stats for cell data; Add saturation to new dict too
    stats = get_cell_stat_info(spipe, samp, stat_dict, cell_df, targeted=targeted)
    stat_dict.update(stats)

    # Explicitly calculate this for combine mode
    stat_dict["mean_reads_per_cell"] = (
        stat_dict["number_of_reads"] / stat_dict[f"{xkey}number_of_cells"]
    )

    # Targeted gene counts
    if targeted:
        eg_df = spipe.read_csv("SFR_EXPRESS_GENE", samp=samp)
        new_dict = get_targeted_stat_info(spipe, stat_dict, eg_df)
        stat_dict.update(new_dict)

    # CRISPR guide counts; Assigned guide files must already have been merged
    if spipe.is_focal_crispr():
        new_dict = get_focal_crispr_stat_info(spipe, samp)
        stat_dict.update(new_dict)

    return stat_dict


#### ---------- Stat calc and collect functions; may be single or combine


def get_samp_well_and_tscp_info(spipe, samp, stat_dict, num_tscp, num_targ_tscp=0):
    """Add read count info to stats

    stat_dict = dict with raw numbers of reads and tscp

    Return dict
    """
    new_dict = {}
    new_dict["sample_well_count"] = samp.num_wells(all_sublibs=True)

    # Transcript / targeted tscp counts
    new_dict["number_of_tscp"] = num_tscp
    if num_targ_tscp:
        new_dict["targeted_number_of_tscp"] = num_targ_tscp

    # Fraction of transcripts for this sample
    new_dict["sample_tscp_fraction"] = (
        new_dict["number_of_tscp"] / stat_dict["raw_number_of_tscp"]
    )

    return new_dict


def get_samp_reads_info(spipe, samp, stat_dict):
    """Get estimated number of reads for sample

    Return dict
    """
    new_dict = {}
    vbc_frac = spipe.get_samp_vbc_count(samp, as_frac=True)
    new_dict["number_of_reads"] = stat_dict["raw_number_of_reads"] * vbc_frac
    return new_dict


#### ---------- Stat calc and collect functions; may be single or combine


def get_dge_tscp_cut_info(spipe, tscp_cut_info):
    """Add tscp cell cutoff info to stats

    tscp_cut_info = dict with tscp cutoff stuff

    Return dict
    """
    new_dict = {}

    # AK - Work around for crispr
    if spipe.is_focal_crispr():
        new_dict["guide_tscp_cutoff"] = tscp_cut_info["cell_tscp_cutoff"]
        new_dict["guide_read_cutoff"] = tscp_cut_info["cell_read_cutoff"]
    else:
        new_dict["cell_tscp_cutoff"] = tscp_cut_info["cell_tscp_cutoff"]
    # May not have extra cutoff stuff; e.g. if given barcodes
    try:
        new_dict["cell_tscp_f01_slope"] = tscp_cut_info["thpoint_f01_slope"]
    except Exception:
        pass

    # Make sure all values are float
    new_dict = {k: float(v) for k, v in new_dict.items()}
    return new_dict


def get_tas_stat_info(spipe, tas_df, targeted=False):
    """Collect transcript assignment based stats

    tas_df is tscp assignment df; Should be masked for sample

    Return dict
    """
    new_dict = {}
    # Extra part of keyword if targeted
    xkey = "targeted_" if targeted else ""

    # Saturation = 1 - (tscps / reads)
    new_dict[f"{xkey}sequencing_saturation"] = 1.0 - len(tas_df) / tas_df["count"].sum()

    # Per species; Regardless if only one
    for s in spipe.get_genome_list():
        # Species subset and number of transcripts
        ts_df = tas_df[tas_df["genome"] == s]
        new_dict[f"{s}_{xkey}number_of_tscp"] = len(ts_df)
        new_dict[f"{s}_{xkey}fraction_exonic"] = ts_df.exonic.mean()

    # Make sure all values are float
    new_dict = {k: float(v) for k, v in new_dict.items()}
    return new_dict


def get_mrtg_count_info(spipe, tas_df, bc_indexes, targeted=False):
    """Get mapped read, tscp, gene counts from transcript-assignment

    tas_df = transcript assignment df; Should be masked for sample
    bc_indexs = indices in cells

    Return dict of stats
    """
    new_dict = {}
    bc_indexes = set(bc_indexes)

    # Extra part of keyword if targeted
    xkey = "targeted_" if targeted else ""
    t_pref = xkey + "total_"
    c_pref = xkey + "cell_"

    # Cell subset to get read counts from
    cell_tas_df = tas_df[tas_df["bc_wells"].isin(bc_indexes)]

    # Counts of everything; 'tcount' = Total count
    c_dict = get_one_mrtg_counts(tas_df, pref=t_pref)
    new_dict.update(c_dict)
    # Counts in cells
    c_dict = get_one_mrtg_counts(cell_tas_df, pref=c_pref)
    new_dict.update(c_dict)

    # If multiple species, add counts for those explicitly
    species = spipe.get_genome_list()
    if len(species) > 1:
        for s in species:
            # Species, everything (i.e. total count)
            ts_df = tas_df[tas_df["genome"] == s].reindex()
            c_dict = get_one_mrtg_counts(ts_df, pref=f"{s}_{t_pref}")
            new_dict.update(c_dict)

            # Species, in cells
            ts_df = cell_tas_df[cell_tas_df["genome"] == s].reindex()
            c_dict = get_one_mrtg_counts(ts_df, pref=f"{s}_{c_pref}")
            new_dict.update(c_dict)

    return new_dict


def get_one_mrtg_counts(tas_df, pref=""):
    """Get counts for mapped read, tscp, gene from given tscp assign df

    tas_df = transcript assignment df; Should be masked (i.e. per sample, targeted)

    Return dict
    """
    new_dict = {}
    new_dict[pref + "mread_count"] = tas_df["count"].sum()
    new_dict[pref + "tscp_count"] = len(tas_df)
    new_dict[pref + "gene_count"] = len(tas_df["gene"].unique())
    return new_dict


def get_cell_stat_info(spipe, samp, stat_dict, cell_df, targeted=False):
    """Cell dataframe based counts / stats

    stat_dict = existing stats value (number_of_reads should be set for sample)
    cell_df = cell dataframe for current sample

    Return dict
    """
    new_dict = {}

    # Possible target-enrichment part of key names
    xkey = "targeted_" if targeted else ""

    # Load count totals for denoms; These keys have no 'targeted'
    rtg_df = spipe.read_csv("SFR_SAMP_MRTG_CT", samp=samp)
    rtg_dict = rtg_df["count"].to_dict()

    # Assumes number of reads has already been set for this sample
    new_dict[f"{xkey}number_of_cells"] = len(cell_df)
    new_dict["mean_reads_per_cell"] = (
        stat_dict["number_of_reads"] / new_dict[f"{xkey}number_of_cells"]
    )

    species = spipe.get_genome_list()
    # Number of multiplet cells
    if len(species) > 1:
        n_mult = len(cell_df[cell_df["species"] == "multiplet"])
        new_dict[f"multiplet_{xkey}number_of_cells"] = n_mult

    # Per species
    for s in species:
        # Species subset and number of transcripts
        sp_cell_df = cell_df[cell_df["species"] == s]
        new_dict[f"{s}_{xkey}number_of_cells"] = len(sp_cell_df)

        # Note; Med of all tscp counts, not species-specific
        v = np.median(sp_cell_df[f"{xkey}tscp_count"].values)
        new_dict[f"{s}_{xkey}median_tscp_per_cell"] = v

        v = (
            new_dict[f"{s}_{xkey}median_tscp_per_cell"]
            * 0.5
            / stat_dict[f"{xkey}sequencing_saturation"]
        )
        new_dict[f"{s}_{xkey}median_tscp_at50"] = v

        # Note; Med of all gene counts, not species-specific
        v = np.median(sp_cell_df[f"{xkey}gene_count"].values)
        new_dict[f"{s}_{xkey}median_genes_per_cell"] = v

        # If one genome, no genome-specific cols (so key is general)
        key = (
            f"{s}_{xkey}total_mread_count"
            if len(species) > 1
            else f"{xkey}total_mread_count"
        )
        num_reads = rtg_dict[key]
        num_reads_in_cells = sp_cell_df[f"{xkey}mread_count"].sum()
        new_dict[f"{s}_{xkey}fraction_reads_in_cells"] = num_reads_in_cells / num_reads

        key = (
            f"{s}_{xkey}total_tscp_count"
            if len(species) > 1
            else f"{xkey}total_tscp_count"
        )
        num_tscp = rtg_dict[key]
        num_tscp_in_cells = sp_cell_df[f"{xkey}tscp_count"].sum()
        new_dict[f"{s}_{xkey}fraction_tscp_in_cells"] = num_tscp_in_cells / num_tscp

    # Make sure all values are float
    new_dict = {k: float(v) for k, v in new_dict.items()}
    return new_dict


def get_read_samp_stat_info(spipe, samp, stat_dict):
    """Add read sampled info to stats

    TSO from previously collected lengths
    map_targs from previously collected file

    Return dict
    """
    new_dict = {}

    # TSO fraction
    new_dict["tso_fraction_in_read1"] = (
        stat_dict["reads_tso_trim"] / stat_dict["reads_valid_barcode"]
    )

    # Map targets
    if spipe.num_map_targs() > 0:
        mapt_df = spipe.read_csv("PF_MAP_TARG_CT", samp=samp)
        # Sample subset
        well_list = samp.get_well_list()
        mapt_df = mapt_df[mapt_df["well"].isin(well_list)]

        # Don't break into _for _rev; Just use combo = _any
        out_cols = [c for c in mapt_df.columns if c.endswith("_count_any")]
        # Primer type list; i.e. R,T
        pt_list = sorted(mapt_df["stype"].unique())
        for col in out_cols:
            # Column base name
            bname = col.replace("_count_any", "")

            denom = mapt_df["count"].sum()
            frac = 0 if (denom < 1) else mapt_df[col].sum() / denom
            key = f"mapto_{bname}_all"
            # print(key, frac)
            new_dict[key] = round(frac, 4)
            # If more than primer type, process per type
            if len(pt_list) > 1:
                for ptype in pt_list:
                    mapt2_df = mapt_df[mapt_df["stype"] == ptype]
                    denom2 = mapt2_df["count"].sum()
                    frac = 0 if (denom2 < 1) else mapt2_df[col].sum() / denom2
                    key = f"mapto_{bname}_{ptype}"
                    # print(key, ptype, frac)
                    new_dict[key] = round(frac, 4)

    # Make sure all values are float
    new_dict = {k: float(v) for k, v in new_dict.items()}
    return new_dict


def get_stat_fract_info(spipe, stat_dict):
    """Calculate misc fractional stats

    Return dict
    """
    new_dict = {}

    # Raw number of reads, as valid barcode not scaled by sample
    new_dict["valid_barcode_fraction"] = (
        stat_dict["reads_valid_barcode"] / stat_dict["raw_number_of_reads"]
    )
    new_dict["transcriptome_map_fraction"] = (
        stat_dict["reads_map_transcriptome"] / stat_dict["reads_valid_barcode"]
    )

    # Make sure all values are float
    new_dict = {k: float(v) for k, v in new_dict.items()}
    return new_dict


def get_targeted_stat_info(spipe, stat_dict, eg_df):
    """Add target enrichment stats

    dg_df = expressed gene dataframe

    Return dict
    """
    new_dict = {}

    # Number of targets
    targen_info = spipe.get_targen_info()
    new_dict["target_number"] = targen_info["target_number"]

    # Nontarget = total - targeted transcript count
    new_dict["nontargeted_number_of_tscp"] = (
        stat_dict["total_tscp_count"] - stat_dict["targeted_total_tscp_count"]
    )
    new_dict["nontargeted_number_of_mreads"] = (
        stat_dict["total_mread_count"] - stat_dict["targeted_total_mread_count"]
    )

    new_dict["targeted_tscp_fraction"] = (
        stat_dict["targeted_total_tscp_count"] / stat_dict["total_tscp_count"]
    )
    new_dict["targeted_mread_fraction"] = (
        stat_dict["targeted_total_mread_count"] / stat_dict["total_mread_count"]
    )

    # Genes with min cells
    min_cells = spipe.get_par_val("dge_gene_cell_targeted_min", as_int=True)
    # Counts with min and more
    cut_list = [min_cells, min_cells * 10]
    counts = utils.get_df_cutoff_counts(eg_df, "cell_count", cutoff=cut_list)
    new_dict["targeted_gene_fraction"] = counts[0] / new_dict["target_number"]
    new_dict["targeted_gene10_fraction"] = counts[1] / new_dict["target_number"]

    # Make sure all values are float
    new_dict = {k: float(v) for k, v in new_dict.items()}
    return new_dict


def get_focal_crispr_stat_info(spipe, samp, mtx=None, max_count=1):
    """Get focal barcoding / crispr stats

    mtx = (sparse) matrix cell x gene (guide); Counts < thresh should be zero
    max_count = max count to tally values for (start at zero, and include over)

    Return dict
    """
    # Get matrix if not given
    if mtx is None:
        f_dict = const.out_filepaths(spipe, samp)
        mtx = utils.read_mtx(f_dict["mtx"])

    mtx_counts = mtx.astype(bool).sum(axis=1)
    new_dict = {}
    n = 0
    while n <= max_count:
        key = f"num_cells_with_{n}_crispr_guide"
        # Result of this is 2D matrix
        new_dict[key] = mtx_counts[mtx_counts == n].shape[1]
        n += 1
    # Last one is over count
    key = f"num_cells_with_{n}plus_crispr_guide"
    new_dict[key] = mtx_counts[mtx_counts >= n].shape[1]

    # Now fracitons
    f_dict = {}
    n_cells = mtx.shape[0]
    for k, v in new_dict.items():
        new_k = k.replace("num_", "frac_")
        f_dict[k] = v
        f_dict[new_k] = v / n_cells

    # Number of guides with any counts
    # gene_info = spipe.get_gene_info()
    # gene_names = gene_info["all_genes"]["gene_name"]
    f_dict["expected_number_guides"] = spipe.num_genes()

    # guide cell count file
    fkey = "SFR_FB_CELL_CT" if spipe.is_focal_crispr(focal=True) else "SFR_CRSP_CELL_CT"
    df = spipe.read_csv(fkey, samp=samp)
    # Key may be 'combined' or 'count'
    if "count" in df.columns:
        c_col = "count"
    elif "combined" in df.columns:
        c_col = "combined"
    else:
        raise ValueError(f"Could not find count in guide count file cols: {df.columns}")
    f_dict["detected_number_guides"] = len(df[df[c_col] > 0])

    return f_dict
