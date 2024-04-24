#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import sys
import traceback
import pandas as pd
import numpy as np
import scipy

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
    import ccutils
    import combine as comb
    import crispr
    import pb_const as const
    import pb_stats
    import utils
else:
    from splitpipe import bcutils
    from splitpipe import ccutils
    from splitpipe import combine as comb
    from splitpipe import crispr
    from splitpipe import pb_const as const
    from splitpipe import pb_stats
    from splitpipe import utils


# ---------------------------------------------------------------------------
def have_dge_ins(spipe, verb=True):
    """Check if DGE input files/data exist

    Return tuple (status, list of files)
    """
    ok = False

    check_files = []
    # Transcript assignment file
    check_files.append(spipe.filepath("PF_TSCP_ASSIGN", None))
    # Stat files if not combine mode
    check_files.append(spipe.filepath("PF_STAT_SEQ", None))
    check_files.append(spipe.filepath("PF_STAT_PIPE", None))
    # If list is good, we have inputs
    if utils.check_infile(check_files, verb=verb):
        ok = True

    # Genome working files; Add name to "file" list for possible later print
    if ok:
        if spipe.is_focal_crispr(map_direct=True):
            check_files.append("spipe.get_guide_info")
            if not spipe.get_guide_info():
                ok = False
        else:
            check_files.append("spipe.get_gene_info")
            if not spipe.get_gene_info():
                ok = False

    # Min number of transcripts; not combine mode
    if ok and not spipe.get_mode("comb"):
        dge_min = spipe.get_par_val("dge_min_start_tscp", as_int=True)
        stat_dict = spipe.get_seq_pipe_stats_dict(verb=False)
        n_tscp = int(stat_dict["number_of_tscp"])
        if n_tscp < dge_min:
            ok = False
            story = f"Number of transcripts {n_tscp} too few for DGE (min {dge_min})"
            check_files = [story]
    return (ok, check_files)


def have_dge_outs(spipe, verb=False, samp=None):
    """Check if DGE output files exist

    Return tuple (status, list of files)
    """
    ok = False
    check_files = []
    # File path keys depend on use case
    dge_fkeys = const.out_filepaths(spipe, use_case=True, samp=None, as_key=True)
    # DGE matrix files
    if samp:
        check_files.append(spipe.filepath(dge_fkeys["mtx"], samp))
    else:
        for samp in spipe.get_samples():
            check_files.append(spipe.filepath(dge_fkeys["mtx"], samp))

    # If list is good, we have outs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return (ok, check_files)


def run_dge(spipe):
    """Run DGE step

    Return status
    """
    spipe.report_proc_step("DGE")
    # Input and output status and list
    i_ok, i_list = have_dge_ins(spipe)
    o_ok, o_list = have_dge_outs(spipe)

    # Run only if don't have output or fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for dge: {bad_list}"
            spipe.set_problem(story)
        else:
            # spipe.report_run_story("Processing global (not per-sample) data")
            # Per-sample processing
            generate_all_dge_outputs(spipe)

        # If all good, aggregate summary files
        if spipe.no_problems():
            pb_stats.aggregate_asum_csv(spipe)
    else:
        spipe.report_run_story("Skipping DGE and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")

    ok = spipe.no_problems()
    if ok:
        dge_clean_up(spipe)

    spipe.report_proc_step("DGE", status=ok)
    return ok


def generate_all_dge_outputs(spipe):
    """Function to loop through every sample and process

    Returns nothing; If any problems, these are reported / remembered
    """
    spipe.report_run_story(f"Processing data for {spipe.num_samples()} samples")
    # Each sample
    for i, samp in enumerate(spipe.get_samples()):
        print()
        story = f"DGE processing sample ({i+1}) {samp.get_name()}"
        spipe.report_run_story(story)
        spipe.report_run_mem(story=f"DGE samp {i+1}")

        # Catch possible problems
        try:
            # Processing different for combine vs normal
            if spipe.get_mode("comb"):
                ok = comb.combine_single_dge_output(spipe, samp)
            else:
                ok = generate_single_dge_output(spipe, samp)
        except Exception as e:
            ok = False
            # Dump traceback, stdout and log
            traceback.print_exc(file=sys.stdout)
            sys.stdout.flush()
            traceback.print_exc(file=spipe.get_log_fh())
            story = f"Failed single_dge_output: {e}"
            spipe.set_problem(story)

        # Problem? report / possibly bail
        if not ok:
            # Problem report, possibly bail
            story = f"DGE problem with sample {samp.get_name()}"
            spipe.set_problem(story)
            samp.set_status(False, "DGE")
            if spipe.keep_going():
                spipe.report_run_story("Keep-going True... So, keep on going!")
                spipe.clear_problems()
                print()
            else:
                break


# ------------------ Per sample processing ----------------------------------


def generate_single_dge_output(spipe, samp):
    """Make DGE outputs for given single sample

    Return status
    """
    if samp.is_meta():
        spipe.report_run_story(f"Meta sample {samp.get_name()}; combining DGE parts")
        ok = gen_meta_sample_dge_output(spipe, samp)
    else:
        spipe.report_run_story(f"Simple sample {samp.get_name()}; calculating DGE")
        ok = gen_simp_sample_dge_output(spipe, samp)
    return ok


def gen_simp_sample_dge_output(spipe, samp):
    """Make DGE outputs for given simple (non-meta) sample

    Return status
    """
    targeted = spipe.is_targeted()
    focal = spipe.is_focal_crispr()
    targ_tas_df = None
    num_targ_tscp = 0

    # Get starting (global) stats collection; Sequencing and pipeline
    stat_dict = spipe.get_seq_pipe_stats_dict()

    # Transcript assignment dataframe, for sample only
    tas_df = spipe.get_tscp_assign_df(samp=samp, focal=focal)

    # Get non filtered num_tscp if focal
    if focal:
        # tas_df_fb_uf = spipe.get_tscp_assign_df(samp=samp)
        tas_df_fb_uf = tas_df

        num_tscp = len(tas_df_fb_uf)
    else:
        num_tscp = len(tas_df)

    # Also targeted only tscp collection?
    if targeted:
        targ_tas_df = spipe.get_tscp_assign_df(samp=samp, targeted=True, focal=focal)
        num_targ_tscp = len(targ_tas_df)

    # DGE outputs; Some data written to file but also returns dict
    dge_var_dict = make_use_case_dge_outputs(spipe, samp, tas_df, targ_tas_df)

    # Sequentialy build up stat fields
    stats = pb_stats.get_samp_well_and_tscp_info(
        spipe, samp, stat_dict, num_tscp, num_targ_tscp
    )
    stat_dict.update(stats)
    stats = pb_stats.get_samp_reads_info(spipe, samp, stat_dict)
    stat_dict.update(stats)

    # Cutoff calculation stuff
    if dge_var_dict:
        stats = pb_stats.get_dge_tscp_cut_info(spipe, dge_var_dict["tscp_cut_info"])
        stat_dict.update(stats)

    # Transcript count stuff (saturation...)
    if focal:
        stats = pb_stats.get_tas_stat_info(spipe, tas_df_fb_uf)
    else:
        stats = pb_stats.get_tas_stat_info(spipe, tas_df)
    if targeted:
        t_stats = pb_stats.get_tas_stat_info(spipe, targ_tas_df, targeted=True)
        stats.update(t_stats)
    stat_dict.update(stats)

    # Counts of reads, tscp, genes; Collect and save to file
    if dge_var_dict:
        bc_indexes = set(dge_var_dict["bc_indexes"])
        if focal:
            stats = pb_stats.get_mrtg_count_info(spipe, tas_df_fb_uf, bc_indexes)
        else:
            stats = pb_stats.get_mrtg_count_info(spipe, tas_df, bc_indexes)
        if targeted:
            t_stats = pb_stats.get_mrtg_count_info(
                spipe, targ_tas_df, bc_indexes, targeted=True
            )
            stats.update(t_stats)
        stat_dict.update(stats)
        write_tas_mrtg_counts(spipe, samp, stats)

    # Info from cell data
    if dge_var_dict:
        stats = pb_stats.get_cell_stat_info(
            spipe, samp, stat_dict, dge_var_dict["cell_df"], targeted=targeted
        )
        stat_dict.update(stats)

    # Info from read sample; TSO, map targets
    stats = pb_stats.get_read_samp_stat_info(spipe, samp, stat_dict)
    stat_dict.update(stats)

    # Add targeted fields from already-collected data
    if targeted:
        stats = pb_stats.get_targeted_stat_info(
            spipe, stat_dict, dge_var_dict["expressed_genes"]
        )
        stat_dict.update(stats)

    # Add focal barcoding fields
    if focal:
        # Threshold to mask too-low tscp elements
        thresh = dge_var_dict["tscp_cut_info"]["cell_tscp_cutoff"]
        # This modifies DGE matrix (less than thresh removed); Don't need later
        mtx = utils.mtx_zero_less_than(dge_var_dict["dge_matrix"], thresh, copy=False)
        stats = pb_stats.get_focal_crispr_stat_info(spipe, samp, mtx)
        stat_dict.update(stats)

    # Add fractions from various already-collected data
    stats = pb_stats.get_stat_fract_info(spipe, stat_dict)
    stat_dict.update(stats)

    # Save analysis summary stats csv file
    pb_stats.save_analysis_summary(spipe, samp, stat_dict)

    # Read and tscp subsampling (saturation plots); Focal barcode has no prefix
    if dge_var_dict:
        spipe.report_run_story2("Getting read subsample stats")
        cell_df = dge_var_dict["cell_df"]
        ss_tas_df = targ_tas_df if targeted else tas_df
        gene_ss_df, tscp_ss_df = proc_subsample_counts(
            spipe, ss_tas_df, cell_df, stat_dict
        )
        # Note: Cannot write as int as may have NA (two species, different numbers)
        spipe.write_df(gene_ss_df, "SFR_SS_GENE_CT", samp=samp)
        spipe.write_df(tscp_ss_df, "SFR_SS_TSCP_CT", samp=samp)

    return True


def make_use_case_dge_outputs(spipe, samp, tas_df, targ_tas_df):
    """Make DGE outputs for given sample based on different use cases

    tas_df = tscp assignment df for sample
    targ_tas_df = targeted tscp assignment df for sample; None if not targeted

    Outputs written to file

    Return dict with DGE outputs
    """
    # Use case dictates processing and outputs
    #   normal           filtered, unfiltered
    #   target_only      targeted_filtered, targeted_unfiltered
    #   target_parent    targeted_filtered, targeted_unfiltered
    #   focal_bc         focal / crispr filtered, unfiltered
    targeted = True if targ_tas_df is not None else False

    # Different use cases yield different DGE outputs
    # Normal case
    if spipe.get_use_case("normal"):
        # No filter, no target genes, and no existing barcode constraints
        spipe.report_run_story2("Getting unfiltered DGE matrix")
        dge_var_dict = make_dge_parts(spipe, samp, tas_df)
        if dge_var_dict:
            write_dge_files(spipe, samp, dge_var_dict, dir_key="DIR_S_DGE_N_UF")

        spipe.report_run_story2("Getting filtered DGE matrix")
        dge_var_dict = make_dge_parts(spipe, samp, tas_df, filtered=True)
        if dge_var_dict:
            write_dge_files(spipe, samp, dge_var_dict, dir_key="DIR_S_DGE_N_F")

        # Sample (report subdir) count files; per species, per well
        if dge_var_dict:
            write_samp_rep_count_files(spipe, samp, tas_df, dge_var_dict)

    # Targeted, no parent
    elif spipe.get_use_case("target_only"):
        # Only get unfiltered for enrichment stats; Don't save anything
        spipe.report_run_story2("Getting unfiltered DGE matrix")
        dge_var_dict = make_dge_parts(spipe, samp, tas_df)
        # Save unfiltered read,tscp,gene counts for enrichment fractions
        if dge_var_dict:
            cols = "mread_count,tscp_count,gene_count".split(",")
            uf_bc_rtg = dge_var_dict["cell_df"][cols]

        spipe.report_run_story2("Getting targeted, unfiltered DGE matrix")
        dge_var_dict = make_dge_parts(spipe, samp, targ_tas_df, targeted=True)
        if dge_var_dict:
            write_dge_files(
                spipe, samp, dge_var_dict, dir_key="DIR_S_DGE_T_UF", targeted=True
            )

        spipe.report_run_story2("Getting targeted, filtered  DGE matrix")
        dge_var_dict = make_dge_parts(
            spipe, samp, targ_tas_df, targeted=True, filtered=True
        )
        if dge_var_dict:
            write_dge_files(
                spipe, samp, dge_var_dict, dir_key="DIR_S_DGE_T_F", targeted=True
            )

        # Sample report count files; Targeted tscp only
        if dge_var_dict:
            write_samp_rep_count_files(
                spipe, samp, targ_tas_df, dge_var_dict, targeted=True
            )
            write_samp_targeted_fracs(spipe, samp, dge_var_dict["cell_df"], uf_bc_rtg)

    # Targeted with parent
    elif spipe.get_use_case("target_parent"):
        # Only get unfiltered for enrichment stats; Don't save anything
        spipe.report_run_story2("Getting unfiltered DGE matrix")
        dge_var_dict = make_dge_parts(spipe, samp, tas_df)
        # Save unfiltered read,tscp,gene counts for enrichment fractions
        if dge_var_dict:
            cols = "mread_count,tscp_count,gene_count".split(",")
            uf_bc_rtg = dge_var_dict["cell_df"][cols]

        spipe.report_run_story2("Getting target_parent, unfiltered DGE matrix")
        dge_var_dict = make_dge_parts(spipe, samp, targ_tas_df, targeted=True)
        if dge_var_dict:
            write_dge_files(
                spipe, samp, dge_var_dict, dir_key="DIR_S_DGE_T_UF", targeted=True
            )

        spipe.report_run_story2("Getting target_parent, filtered  DGE matrix")
        dge_var_dict = make_dge_parts(
            spipe, samp, targ_tas_df, targeted=True, filtered=True
        )
        if dge_var_dict:
            write_dge_files(
                spipe, samp, dge_var_dict, dir_key="DIR_S_DGE_T_F", targeted=True
            )

        # Sample report count files; Targeted tscp only
        if dge_var_dict:
            write_samp_rep_count_files(
                spipe, samp, targ_tas_df, dge_var_dict, targeted=True
            )
            write_samp_targeted_fracs(spipe, samp, dge_var_dict["cell_df"], uf_bc_rtg)

        spipe.report_run_story2("Getting focal / CRISPR data matrix")
    # Focal / CRISPR case
    elif spipe.is_focal_crispr():
        # Process as non-filtered (so all barcodes from parent)
        #   and non-targeted (so no "target_" keys, etc)
        # Starting point transcripts may or may not be targeted
        if targeted:
            use_tas_df = targ_tas_df
        else:
            use_tas_df = tas_df

        thresh = 1
        # thresh = spipe.get_par_val("dge_tscp_read_fb_min", as_int=True)

        dge_var_dict = make_dge_parts(spipe, samp, use_tas_df, targeted=targeted)
        if dge_var_dict:
            # save cutoff info
            dge_var_dict["tscp_cut_info"] = get_crispr_cut_info(spipe)
            # Unfiltered DGE without applying threshold to matrix
            fkey = (
                "DIR_S_DGE_FB_UF"
                if spipe.is_focal_crispr(focal=True)
                else "DIR_S_DGE_CRSP_UF"
            )
            write_dge_files(spipe, samp, dge_var_dict, dir_key=fkey, targeted=targeted)

        # Sample report count files; Focal guide mapping files
        if dge_var_dict:
            # Utilize automatic crispr thresholding to filter df
            # If given manual thresholds, use those and no automated thresholding
            given_cell_bc = spipe.get_cell_bc(samp)
            given_cell_bc = list(given_cell_bc)
            use_tas_df = use_tas_df[use_tas_df["bc_wells"].isin(given_cell_bc)]
            min_count_read = spipe.get_par_val("dge_tscp_read_fb_min", as_int=True)
            min_count_tscp = spipe.get_par_val("dge_gene_tscp_fb_min", as_int=True)
            if (min_count_read > 1) or (min_count_tscp > 1):
                spipe.report_run_story2("Manual thresholds provided for filtering")
                story = f"Filtering using {min_count_read} read cutoff and {min_count_tscp} tscp cutoff"
                spipe.report_run_story2(story)

                use_tas_df_f, read_cutoff, tscp_cutoff = crispr.crispr_get_cutoff(
                    use_tas_df,
                    manual=True,
                    manual_read_cutoff=min_count_read,
                    manual_tscp_cutoff=min_count_tscp,
                )

                # Remake mtx and objects using filtered data
                dge_var_dict = make_dge_parts(
                    spipe, samp, use_tas_df_f, targeted=targeted
                )
                dge_var_dict["tscp_cut_info"] = get_crispr_cut_info(
                    spipe
                )  # Refill the required data
            # If no manual thresholds given (default), perform automated thresholding
            else:
                spipe.report_run_story2(
                    "Estimating read and transcript cutoffs for filtering guides"
                )
                use_tas_df_f, read_cutoff, tscp_cutoff = crispr.crispr_get_cutoff(
                    use_tas_df
                )

                # Remake parts using filtered data
                dge_var_dict = make_dge_parts(
                    spipe, samp, use_tas_df_f, targeted=targeted
                )

                # Update crispr thresholds
                crispr_info = {
                    "algorithm": "Automatic thresholds",
                    "cell_mm_cutoff_raw": spipe.get_par_val(
                        "dge_gene_tscp_fb_min", as_int=True
                    ),
                    "cell_tscp_cutoff": tscp_cutoff,
                    "cell_read_cutoff": read_cutoff,
                    "num_in_tscp": 0,
                    "num_out_tscp": 0,
                }
                dge_var_dict["tscp_cut_info"] = crispr_info

            # Zero out below-threshold matrix elements
            thresh = dge_var_dict["tscp_cut_info"]["cell_tscp_cutoff"]
            spipe.report_run_story2(f"Masking (focal) DGE elements < {thresh} to zero")
            mtx = dge_var_dict["dge_matrix"]
            story = utils.mtx_stat_str(mtx, min_list=[1, 3, thresh])
            print("FFF1", story)
            mtx = utils.mtx_zero_less_than(mtx, thresh, copy=False)
            story = utils.mtx_stat_str(mtx, min_list=[1, 3, thresh])
            print("FFF2", story)
            dge_var_dict["dge_matrix"] = mtx

            # Filtered DGE with thresholded matrix
            fkey = (
                "DIR_S_DGE_FB_F"
                if spipe.is_focal_crispr(focal=True)
                else "DIR_S_DGE_CRSP_F"
            )
            write_dge_files(spipe, samp, dge_var_dict, dir_key=fkey, targeted=targeted)

            write_samp_rep_count_files(
                spipe, samp, use_tas_df, dge_var_dict, targeted=targeted
            )
            crispr.write_guide_cell_map_files(spipe, samp, dge_var_dict, thresh)

    else:
        use_case = spipe.get_use_case()
        story = f"generate_single_dge_output unexpected use_case {use_case}"
        raise ValueError(story)

    return dge_var_dict


def make_dge_parts(spipe, samp, tas_df, filtered=False, targeted=False):
    """Generate DGE matrix and associated objects

    tas_df = transcript assignment df; Should be masked for sample, target enrichment
    filtered = flag to generate filtered matrix; i.e. calc cell cutoff
    targeted = flag to use only targeted subset of genes (if exists)

    Return dict with matrix, expressed genes, bc_indexes (cells), cell_df, ... more
    """
    # Have cell barcodes already (i.e. from parent)
    if spipe.have_cell_bc():
        given_bc = True
        cellbc_df = spipe.get_cell_bc(samp)
    else:
        given_bc = False
        cellbc_df = None

    # Story to report
    filt_part = "Filtered" if filtered else "Unfiltered"
    targ_part = "Genes targeted" if targeted else "Genes all"
    cbc_part = "Barcodes given" if given_bc else "Barcodes TBD"
    description = f"{samp.get_name()} {filt_part}, {targ_part}, {cbc_part}"
    spipe.report_run_story(f">> Getting DGE; {description}")

    tscp_cut_info = {}
    tscp_thresh = 0
    # tscp per cell; Barcodes are dataframe index
    tscp_counts_df = ccutils.bc_tscp_count_df(
        tas_df, cell_bc=cellbc_df, sort_index=True
    )
    spipe.report_run_story2(f"Raw transcripts per cell {tscp_counts_df.shape}")

    # Transcript threshold and cells (bc indexes) to keep
    # Given cell barcodes already? (This may be None)
    if given_bc:
        bc_indexes = tscp_counts_df.index.values
        spipe.report_run_story2(
            f"Cells given and filter, so use cells {len(bc_indexes)}"
        )
        tscp_cut_info["cell_call"] = "given_barcodes"
        # Threshold is still used for filtered case when given bc
        if filtered:
            tscp_thresh = spipe.get_par_val("dge_tscp_cut_targeted_min", as_int=True)

    # Use tscp counts to define cells
    else:
        spipe.report_run_story2("Getting tscp cutoff to call cells")
        # Get barcode and gene thresholds. Depends filtering or not
        if filtered:
            ok, tscp_cut_info = spipe.get_tscp_cutoff(
                tscp_counts_df, samp, targeted=targeted
            )
            if not ok:
                story = f"No filtered DGE for sample {samp.get_name()} (cutoff failed)"
                spipe.report_warning("Failed to get cutoff; NO DGE")
                samp.set_status(False, "DGE")
                return None
            tscp_thresh = tscp_cut_info["cell_tscp_cutoff"]
        else:
            if targeted:
                tscp_thresh = spipe.get_par_val(
                    "dge_tscp_cut_targeted_umin", as_int=True
                )
            else:
                tscp_thresh = spipe.get_par_val("dge_tscp_cut_umin", as_int=True)

    tscp_cut_info["cell_tscp_cutoff"] = tscp_thresh
    # Filter by tscp cutoff
    tscp_counts_df = tscp_counts_df[tscp_counts_df["tscp_count"] >= tscp_thresh]
    bc_indexes = tscp_counts_df.index.values
    spipe.report_run_story2(
        f"Cells filtered by tscp cutoff {tscp_thresh}: cells {len(bc_indexes)}"
    )

    # Special-case specific inputs ...
    xkey = ""
    species = spipe.get_genome_list()
    if spipe.is_focal_crispr(map_direct=True):
        guide_info = spipe.get_guide_info()
        gene_id_to_genome = guide_info["gene_id_to_genome"]
        all_genes = guide_info["all_genes"]
        gene_dict = guide_info["gene_id_to_idx"]
    elif targeted:
        gene_info = spipe.get_gene_info()
        gene_id_to_genome = gene_info["gene_id_to_genome"]
        # Target gene info
        targen_info = spipe.get_targen_info()
        all_genes = targen_info["target_genes"]
        gene_dict = targen_info["gene_id_to_idx"]
        # xkey for extra part of targeted data keys
        xkey = "targeted_"
    else:
        gene_info = spipe.get_gene_info()
        gene_id_to_genome = gene_info["gene_id_to_genome"]
        all_genes = gene_info["all_genes"]
        gene_dict = gene_info["gene_id_to_idx"]

    ### Build DGE matrix
    # Collect barcode-gene pairs for matrix
    cell_dict = dict(zip(bc_indexes, range(len(bc_indexes))))
    rows, cols, vals = [], [], []
    for bci, g in zip(tas_df["bc_wells"], tas_df["gene"]):
        try:
            cell_dict[bci]
        except Exception:
            pass
        else:
            rows.append(cell_dict[bci])
            cols.append(gene_dict[g])
            vals.append(1)
    rows.append(len(cell_dict) - 1)
    cols.append(len(gene_dict) - 1)
    vals.append(0)
    # int32 goes to 2 billion ... should be enough
    dge_matrix = scipy.sparse.csr_matrix((vals, (rows, cols)), dtype=np.int32)

    # Filter cells to any threshold
    if tscp_thresh > 0:
        # Subset of bc_indexes (cells) with number of transcripts passing cutoff
        # Sum tscp counts across genes (col = sum(1))
        cell_mask = np.array(dge_matrix.sum(1)).flatten() >= tscp_thresh
        dge_matrix = dge_matrix[cell_mask, :]
        bc_indexes = tscp_counts_df.index.values[cell_mask]
    spipe.report_run_story2(f"dge_matrix {dge_matrix.shape}")

    # Table of expressed genes
    expressed_genes = get_expressed_genes_df(dge_matrix, all_genes)
    story = f"expressed_genes {description}, genes={len(expressed_genes)}"
    spipe.report_run_story2(story)

    ####
    # Collect (separated by species) counts of reads, tscps, genes
    # Cell subset to get read counts from
    tas2_df = tas_df[tas_df["bc_wells"].isin(set(bc_indexes))]

    species_gene_inds = {}
    species_read_counts = {}
    species_tscp_counts = {}
    species_gene_counts = {}
    for s in species:
        species_gene_inds[s] = np.where(
            all_genes["gene_id"].apply(lambda s: gene_id_to_genome[s]) == s
        )[0]

        species_read_counts[s] = (
            tas2_df[tas2_df["genome"] == s].groupby("bc_wells")["count"].sum()
        )
        species_tscp_counts[s] = (
            tas2_df[tas2_df["genome"] == s].groupby("bc_wells")["count"].size()
        )
        species_gene_counts[s] = pd.Series(
            index=bc_indexes,
            data=np.array((dge_matrix[:, species_gene_inds[s]] > 0).sum(1)).flatten(),
        )

    # Create count dataframes; save
    species_read_counts = (
        pd.DataFrame.from_dict(species_read_counts).fillna(0).astype(int)
    )
    species_read_counts = species_read_counts.reindex(bc_indexes)
    species_read_counts.index.name = "bc_wells"

    species_tscp_counts = (
        pd.DataFrame.from_dict(species_tscp_counts).fillna(0).astype(int)
    )
    species_tscp_counts = species_tscp_counts.reindex(bc_indexes)
    species_tscp_counts.index.name = "bc_wells"

    species_gene_counts = pd.DataFrame.from_dict(species_gene_counts)
    species_gene_counts = species_gene_counts.reindex(bc_indexes)
    species_gene_counts.index.name = "bc_wells"

    # Species assignment
    # Init all to 'multiplet', reset subset barcode (cell) has more than purity thresh
    species_assignments = pd.Series(["multiplet" for i in range(len(bc_indexes))])
    sp_purity = spipe.get_par_val("dge_species_purity_thresh", as_float=True)
    for s in species:
        species_assignments.loc[
            np.where((species_tscp_counts[s] / species_tscp_counts.sum(1)) > sp_purity)
        ] = s

    ####
    # Cell dataframe
    cell_df = pd.DataFrame(index=bc_indexes)
    cell_df.index.name = "bc_wells"
    cell_df["sample"] = samp.get_name()
    cell_df["species"] = species_assignments.values
    # Count column names depend on targeted or not
    cell_df[f"{xkey}gene_count"] = np.array((dge_matrix > 0).sum(1)).flatten()
    cell_df[f"{xkey}tscp_count"] = np.array(dge_matrix.sum(1)).flatten()
    # Add mapped read count col (Uses bc as index for df alignment)
    cell_df[f"{xkey}mread_count"] = species_read_counts.sum(1)

    # If multiple species, add counts for those explicitly
    species = spipe.get_genome_list()
    if len(species) > 1:
        for s in species:
            gene_col = f"{s}_{xkey}gene_count"
            cell_df[gene_col] = species_gene_counts[s]
            tscp_col = f"{s}_{xkey}tscp_count"
            cell_df[tscp_col] = species_tscp_counts[s]
            read_col = f"{s}_{xkey}mread_count"
            cell_df[read_col] = species_read_counts[s]

    # Add extra, seaparated well index and label cols
    bc_info = spipe.get_bc_info()
    cell_df = bcutils.sep_df_bci_cols(
        cell_df, bc_info=bc_info, bc_wind=True, bc_well=True
    )

    # If all-well sample, (re)set sample name to other samples
    if samp.is_allwell() and not spipe.get_mode("comb"):
        # If multi-plate, two barcodes
        wsamp_dict = spipe.get_wind_to_sample_dict()
        if samp.is_multi_plate():
            pass
        else:
            cell_df["sample"] = cell_df["bc1_wind"].map(wsamp_dict)

    # In case any zero count 'cells' (e.g. via given list), reset species
    cell_df.loc[cell_df[f"{xkey}tscp_count"] < 1, "species"] = "NA"
    # Make sure all values set; Counts as int, sample as something
    count_cols = [c for c in cell_df.columns if c.endswith("_count")]
    cell_df[count_cols] = cell_df[count_cols].fillna(0).astype(int)
    aw_name = spipe.get_allwell_sample(as_name=True)
    cell_df["sample"] = cell_df["sample"].fillna(aw_name)

    ####
    # Package to return
    new_dict = {
        "dge_matrix": dge_matrix,
        "bc_indexes": bc_indexes,
        "cell_df": cell_df,
        "expressed_genes": expressed_genes.sort_values("gene_id"),
        "tscp_cut_info": tscp_cut_info,
        "species_read_counts": species_read_counts,
        "species_tscp_counts": species_tscp_counts,
        "species_assignments": species_assignments,
    }

    story = (
        f"<< DGE matrix {description} {dge_matrix.shape}, bc_indexes {bc_indexes.shape}"
    )
    spipe.report_run_story2(story)
    return new_dict


def get_expressed_genes_df(dge_mtx, all_genes):
    """Get expressed gene dataframe

    dge_mtx = (sparse) DGE matrix
    all_genes = dataframe with all starting genes

    Return dataframe
    """
    # Subset of genes with at least min bc_indexes; No masking of matrix
    # Cell count = sum non-zero counts across rows (cells), transpose (col in dge, row in all_genes)
    eg_df = all_genes.copy()
    eg_ncell = (dge_mtx != 0).sum(0).T.astype(int)
    eg_ntscp = dge_mtx.sum(0).T.astype(int)
    eg_df["cell_count"] = pd.DataFrame(eg_ncell, columns=["cell_count"])
    eg_df["tscp_count"] = pd.DataFrame(eg_ntscp, columns=["tscp_count"])
    # Have to have at least one cell
    eg_df = eg_df[eg_df["cell_count"] > 0]
    return eg_df


def write_dge_files(spipe, samp, dge_var_dict, dir_key, targeted=False):
    """Write matrix, gene, cell files for given sample

    dge_var_dict = dict with dge stuff; barcode indexes, matrix
    dir_key = DGE output directory key
    targeted = flag to output only targeted genes (vs all)

    Return nothing
    """
    dge_matrix = dge_var_dict["dge_matrix"]
    cell_df = dge_var_dict["cell_df"]

    # Only if filtered (normal or targeted)
    if (dir_key == "DIR_S_DGE_N_F") or (dir_key == "DIR_S_DGE_T_F"):
        exp_genes = dge_var_dict["expressed_genes"]
    else:
        exp_genes = None

    # Filenames depend on dge dir
    file_keys = const.out_filepaths(
        spipe, samp, dir_key, as_key=True, targeted=targeted
    )

    # Sparse matrix
    story = utils.mtx_stat_str(dge_matrix)
    spipe.report_run_story2(story)
    mtx_fname = spipe.filepath(file_keys["mtx"], samp)
    spipe.write_mtx(dge_matrix, mtx_fname)

    # Copy all / target genes file to DGE output dir
    if spipe.get_par_val("dge_copy_dge_genes_file", as_bool=True):
        # cp_src = spipe.filepath(file_keys['src_gene'], samp)
        fkey = "PF_TARGET_GENES" if targeted else "PF_ALL_GENES"
        cp_src = spipe.filepath(fkey, samp)
        cp_dst = spipe.filepath(file_keys["gene"], samp)
        story = utils.copy_file(cp_src, cp_dst)
        spipe.report_run_story2(story)
    else:
        spipe.report_run_story2(
            f"Not copying {cp_src} to DGE dir for {samp.get_name()}"
        )

    # out_cols = [c for c in cell_df.columns if not c.endswith('_wind')]
    out_cols = cell_df.columns
    spipe.write_df(cell_df[out_cols], file_keys["cell"], samp=samp, index=True)

    # Save expressed genes file?
    if spipe.get_par_val("dge_save_express_gene_file", as_bool=True):
        if exp_genes is not None:
            spipe.write_df(exp_genes, "SFR_EXPRESS_GENE", samp=samp, index=False)


def write_samp_targeted_fracs(spipe, samp, cell_df, uf_bc_rtg):
    """Write target enrichment fractions for sample cells

    cell_df = cell dataframe
    uf_bc_rtg = unfiltered barcode dataframe with read, tscp, gene count cols

    Returns nothing
    """
    assert (
        spipe.is_targeted()
    ), f"write_samp_targeted_fracs only for targeted case ({spipe.get_use_case()})"

    # Per cell output; Unfiltered bc may have different
    # Copy in case setting NA to values
    # df = cell_df.join(uf_bc_rtg).copy()
    df = cell_df.join(uf_bc_rtg)

    # col_names = "targeted_mread_count,targeted_tscp_count,targeted_gene_count".split(',')
    # df[col_names].fillna(0, inplace=True)

    df["mread_fraction"] = df["targeted_mread_count"] / df["mread_count"]
    df["tscp_fraction"] = df["targeted_tscp_count"] / df["tscp_count"]
    df["gene_fraction"] = df["targeted_gene_count"] / df["gene_count"]
    # Limit precision of values
    col_names = "mread_fraction,tscp_fraction,gene_fraction".split(",")
    for col in col_names:
        df[col] = df[col].map(lambda x: "%2.4f" % x)
    spipe.write_df(df[col_names], "SFR_ENRICHMENT", samp=samp, index=True)


def tscp_or_cell_count_by_well(
    spipe, samp, cell_df, bc_round, do_cell=False, targeted=False
):
    """Get counts of cells / median tscp per well

    cell_df = cell dataframe
    bc_round = barcode round; 1, 2, 3
    do_cell = flag to get cell counts; Default = tscp
    targeted = flag to denote if targeted

    Return dataframe of counts / medians per well
    """
    # Column with wells for this round
    well_col = f"bc{bc_round}_well"
    # Count col may be targeted
    count_col = "targeted_tscp_count" if targeted else "tscp_count"

    well_list = []
    count_list = []
    for well in list(cell_df[well_col].unique()):
        # cells with this well
        c_df = cell_df[cell_df[well_col] == well]
        # number of cells or median of tscp
        if do_cell:
            n = len(c_df)
        else:
            n = c_df[count_col].median()
        well_list.append(well)
        count_list.append(n)

    df = pd.DataFrame({"well": well_list, "count": count_list})
    df.set_index("well", inplace=True)
    df["count"] = df["count"].astype(int)
    return df


def write_samp_rep_count_files(spipe, samp, tas_df, dge_var_dict, targeted=False):
    """Write out various (sample report subdir) count files

    Transcript cutoff info and counts
    Per species read and tscp counts

    Returns nothing
    """
    # Transcript cutoff calc info
    tscp_cut_info = dge_var_dict["tscp_cut_info"]
    if tscp_cut_info:
        # Round tscp_info values
        for k, v in tscp_cut_info.items():
            if not isinstance(v, str):
                v = round(v, 4)
            tscp_cut_info[k] = v
        df = pd.DataFrame.from_dict(tscp_cut_info, orient="index", columns=["value"])
        df.index.name = "statistic"
        spipe.write_df(df, "SFR_TSCP_CUTOFF", samp=samp, index=True, dump=True)

    # transcript counts per barcode
    tscp_counts = ccutils.bc_tscp_count_df(tas_df)
    tscp_counts = tscp_counts.rename(columns={"tscp_count": "count"})["count"]
    spipe.write_df(tscp_counts.astype(int), "SFR_TSCP_CT", samp=samp, index=True)

    # Per species counts
    species_read_counts = dge_var_dict["species_read_counts"].fillna(0)
    species_tscp_counts = dge_var_dict["species_tscp_counts"].fillna(0)
    spipe.write_df(
        species_read_counts.astype(int), "SFR_SPEC_READ_CT", samp=samp, index=True
    )
    spipe.write_df(
        species_tscp_counts.astype(int), "SFR_SPEC_TSCP_CT", samp=samp, index=True
    )

    # Per well counts
    cell_df = dge_var_dict["cell_df"]
    write_samp_per_well_count_files(spipe, samp, cell_df, targeted=targeted)


def write_samp_per_well_count_files(spipe, samp, cell_df, targeted=False):
    """Write out per-well counts of cells / medians of transcripts for sample

    Return nothing
    """
    for c in range(2):
        if c == 0:
            do_cell = False
            ok_pref = "SFR_TSCP"
        else:
            do_cell = True
            ok_pref = "SFR_CELL"
        # Three rounds
        for r in range(3):
            df = tscp_or_cell_count_by_well(
                spipe, samp, cell_df, r + 1, do_cell, targeted=targeted
            )
            okey = ok_pref + f"_R{r+1}W"
            spipe.write_df(df, okey, samp=samp, index=True)


def write_tas_mrtg_counts(spipe, samp, cnt_dict):
    """Write mapped read, tscp, gene counts to file

    Returns nothing
    """
    # Dataframe and save
    new_df = pd.DataFrame.from_dict(cnt_dict, orient="index", columns=["count"])
    new_df.index.name = "statistic"
    spipe.write_df(new_df, "SFR_SAMP_MRTG_CT", samp=samp)


def get_crispr_cut_info(spipe):
    """Get info ready to report

    Return dict
    """
    new_dict = {
        "algorithm": "fixed thresholds",
        "cell_mm_cutoff_raw": spipe.get_par_val("dge_gene_tscp_fb_min", as_int=True),
        "cell_tscp_cutoff": int(spipe.get_par_val("dge_gene_tscp_fb_min", as_int=True)),
        "cell_read_cutoff": int(spipe.get_par_val("dge_tscp_read_fb_min", as_int=True)),
        "num_in_tscp": 0,
        "num_out_tscp": 0,
    }
    return new_dict


# ------------------- Subsampling -------------------------------------------


def get_subsamp_x_list(max_x, step, n_step):
    """Get sampling points, zero to max (plot on X axis)

    Return list of x sample positions (for read subsampling)
    """
    nx = []
    s = 0
    x = step
    while x < max_x:
        nx.append(x)
        x += step
        s += 1
        if n_step and (s >= n_step):
            s = 0
            step *= 2
    if (x - step) < max_x:
        nx.append(max_x)
    return nx


def proc_subsample_counts(spipe, tas_df, cell_df, stat_dict):
    """Get subsampling count stats for reads

    tas_df = transcript assignment df; Should be masked for sample, etc

    Return tuple of dataframes (gene_ss_df, tscp_ss_df)
    """
    species_assignments = cell_df["species"]
    num_reads = stat_dict["number_of_reads"]
    # X-axis sampling parameters
    ssx_step = spipe.get_par_val("dge_ssc_xstep_size", as_int=True)
    ssx_n_step = spipe.get_par_val("dge_ssc_xstep_ndup", as_int=True)

    # Subsample reads
    species_read_proportions = (
        tas_df.groupby("genome").size() / tas_df.groupby("genome").size().sum()
    )
    # Dicts to collect series per species
    gene_counts_subsampled_dic = {}
    tscp_counts_subsampled_dic = {}

    species = species_assignments.unique()
    for s in species:
        if s == "multiplet":
            continue

        if s == "NA":
            continue

        n_cells = len(cell_df[cell_df["species"] == s])
        if n_cells < 1:
            continue

        # Species-normalized reads per cell = species (max) depth
        seq_depth = species_read_proportions[s] * num_reads / n_cells
        subsample_depths = np.array(get_subsamp_x_list(seq_depth, ssx_step, ssx_n_step))

        s_tas_df = tas_df[tas_df["genome"] == s]
        s_bc_indexes = cell_df[cell_df["species"] == s].index.values

        gene_counts, tscp_counts = subsamp_tas_gene_tscp(
            spipe, s_tas_df, s_bc_indexes, subsample_depths, what=s
        )

        gene_counts_subsampled_dic[s] = pd.Series(gene_counts).fillna(0)
        tscp_counts_subsampled_dic[s] = pd.Series(tscp_counts).fillna(0)

    # Dataframes from one or more species series; Output 'subsample' index
    gene_ss_df = pd.DataFrame.from_dict(gene_counts_subsampled_dic)
    gene_ss_df.index.name = "subsample"
    tscp_ss_df = pd.DataFrame.from_dict(tscp_counts_subsampled_dic)
    tscp_ss_df.index.name = "subsample"

    return gene_ss_df, tscp_ss_df


def subsamp_tas_gene_tscp(spipe, tas_df, bc_indexes, subsample_depths, what=""):
    """Subsample transcript assignment data to get tscp and gene count arrays

    tas_df = tscp assignment data; Already filtered for species, sample, etc
    bc_indexes = array of cell indexes
    subsample_depths = array of depths to sample
    what = value / flag for reporting

    Return two count dicts; gene, tscp
    """
    if what:
        spipe.report_run_story2(
            f"Subsampling {len(subsample_depths)} depths to {round(subsample_depths[-1])} for {what}"
        )

    ss_gene_counts = {0: 0}
    ss_tscp_counts = {0: 0}
    for i, ss_depth in enumerate(subsample_depths):
        if what:
            spipe.report_run_story3(f"Subsample depth[{i}] = {ss_depth}")

        # The last element has the full depth
        ss_fraction = ss_depth / subsample_depths[-1]
        # This is a way to subsample counts; x counts > 0 gives subset of tscp
        sub_sampled_counts = np.random.binomial(tas_df["count"], ss_fraction)

        ss_gene_counts[ss_depth] = (
            tas_df[sub_sampled_counts > 0]
            .groupby("bc_wells")
            .gene.apply(lambda x: len(np.unique(x)))
            .reindex(bc_indexes)
            .median()
        )

        ss_tscp_counts[ss_depth] = (
            tas_df[sub_sampled_counts > 0]
            .groupby("bc_wells")
            .polyN.size()
            .reindex(bc_indexes)
            .median()
        )

    return ss_gene_counts, ss_tscp_counts


# ------------- Combine meta sample parts -----------------------------------


def gen_meta_sample_dge_output(spipe, samp):
    """Make DGE outputs for given single meta sample

    Return status
    """
    assert samp.is_meta(), f"Sample is not meta {samp.get_name()}"

    # Target enrichment
    targeted = spipe.is_targeted()

    # List of DGE dir keys
    dir_key_list = const.use_case_out_dir_key(spipe, as_list=True)
    for dir_key in dir_key_list:
        spipe.report_run_story2(f"Combining DGE files; {dir_key}")
        if not comb.merge_multi_dge_files(spipe, samp, dir_key, combine=False):
            spipe.set_problem(
                f"Failed to combine DGE files for sample {samp.get_name()}"
            )
            return False

    # Expressed genes files?
    if spipe.get_par_val("dge_save_express_gene_file", as_bool=True):
        comb.merge_multi_exgene_files(spipe, samp, combine=False)

    # Analysis_summary and all_stats files; Some fields sum, some mean...
    #   Weights are for weighted average (number of reads) for mean aggregations
    stat_df, weights = comb.merge_multi_stats_df(spipe, samp, combine=False)
    asum_keys, _ = pb_stats.get_asum_out_keys(stat_df.index)

    # Dict with values from new merged column
    merge_col = stat_df.columns[0]
    stat_dict = stat_df.loc[:, merge_col].to_dict()

    # per-sample report files
    # Per barcode ('cell') counts
    comb.merge_multi_bc_counts(spipe, samp, combine=False)
    # Counts per well, etc
    comb.merge_multi_count_files(spipe, samp, weights, combine=False)
    # Focal / CRISPR guide counts
    if spipe.is_focal_crispr():
        comb.merge_multi_focal_files(spipe, samp, combine=False)

    # Subsample de novo from parts
    tas_df = spipe.get_tscp_assign_df(samp=samp, targeted=targeted)
    # New cell metadata for this sample
    fname_dict = const.out_filepaths(spipe, samp, targeted=targeted)
    cell_df = utils.read_csv(fname_dict["cell"]).set_index("bc_wells")
    # Get subsample data then save
    gene_ss_df, tscp_ss_df = proc_subsample_counts(spipe, tas_df, cell_df, stat_dict)
    # Note: Cannot write as int as may have NA (two species, different numbers)
    spipe.write_df(gene_ss_df, "SFR_SS_GENE_CT", samp=samp)
    spipe.write_df(tscp_ss_df, "SFR_SS_TSCP_CT", samp=samp)

    # Some stats get recalculated de novo
    spipe.report_run_story2("Recalculating some merged stats")
    up_dict = pb_stats.updated_stats_dict(spipe, samp, stat_dict, targeted=targeted)
    for k, v in up_dict.items():
        stat_df.loc[k, merge_col] = v

    pb_stats.save_analysis_summary(spipe, samp, stat_df, col="combined")

    return True


def dge_clean_up(spipe):
    """Clean up DGE step things"""
    # Don't need tscp assignment loaded any more
    spipe.free_tscp_assign()
