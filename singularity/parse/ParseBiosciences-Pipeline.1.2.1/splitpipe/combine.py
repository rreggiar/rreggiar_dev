#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import re
import pandas as pd
import numpy as np


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
    import crispr
    import dge
    import pb_const as const
    import pb_stats
    import tcr
    import utils
else:
    from splitpipe import bcutils
    from splitpipe import crispr
    from splitpipe import dge
    from splitpipe import pb_const as const
    from splitpipe import pb_stats
    from splitpipe import tcr
    from splitpipe import utils


# ---------------------------------------------------------------------------
#
def have_combine_ins(spipe, verb=True):
    """Check if combine mode inputs exist"""
    # Return tuple (status, list of files)
    return True, []


def have_combine_outs(spipe, verb=True):
    return True, []


def run_combine(spipe):
    """Run combine step

    Return status
    """
    assert spipe.get_mode("comb"), "run_combine is for comb mode"
    spipe.report_proc_step("combine")


    # TCR handles combine directly
    if spipe.is_tcr():
        ok = tcr.run_tcr(spipe)
    else:
        spipe.report_run_story("Processing global (not per-sample) data")
        comb_global_files(spipe)
        # Per-sample processing
        dge.generate_all_dge_outputs(spipe)
        # If all good, aggregate summary files
        if spipe.no_problems():
            pb_stats.aggregate_asum_csv(spipe)
        ok = spipe.no_problems()
        if ok:
            combine_clean_up(spipe)

    spipe.report_proc_step("combine", status=ok)
    return ok


# -------------- Combine mode, one time (not per sample) --------------------
#
def comb_global_files(spipe):
    """Combine global (not per sample) files for sublibraies

    Returns nothing; Errors may set problems that are stored internally
    """
    # For CRISPR (no genome), copy all_gene file from sublib
    if spipe.is_focal_crispr():
        fpath = spipe.filepath("PF_ALL_GENES", None)
        need_cp = bool(not utils.check_infile(fpath, verb=False))
        if need_cp:
            # Source from first sublib; May not exist yet, so from_par=True
            top_dir = spipe.get_sublibs(from_par=True)[0]
            in_path = spipe.filepath("PF_ALL_GENES", None, top_dir=top_dir)
            story = utils.copy_file(in_path, fpath)
            if not story:
                spipe.set_problem(f"Failed to copy gene data {in_path} to {fpath}")
                return

    # Transcript assignment files; Optional, but not TCR
    if spipe.get_par_val("comb_new_tscp", as_bool=True) and (not spipe.is_tcr()):
        # Output filename
        o_tas_fname = spipe.filepath("PF_TSCP_ASSIGN", None)
        tot_tscp = 0
        tot_bc = 0
        spipe.report_run_story2(
            f"Combining from tscp assign files from {spipe.num_sublibs()} sublibs"
        )
        for i, sublib in enumerate(spipe.get_sublibs()):
            s_top_dir = sublib.get_path()
            tas_df = spipe.get_tscp_assign_df(top_dir=s_top_dir, sub_ind=i, fresh=True)
            tot_tscp += len(tas_df)
            # First one, write, else append
            if i == 0:
                story = utils.write_df(
                    tas_df, o_tas_fname, index=False, compression="gzip", verb=True
                )
            else:
                story = utils.write_df(
                    tas_df,
                    o_tas_fname,
                    header=False,
                    index=False,
                    compression="gzip",
                    verb=True,
                    mode="a",
                )
            spipe.report_run_story2(story)
        spipe.report_run_story2(f"Combined total {tot_tscp} tscp in {tot_bc} barcodes")
    else:
        spipe.report_run_story2("Not combining tscp assign files (comb_new_tscp False)")

    # Valid barcode counts
    df, n = merge_multi_csvs_df(spipe, "PF_VAL_BC_DATA", None, axis=1, combine=True)
    df = aggregate_col_df(df, "count", mode="sum")
    spipe.write_df(df.astype(int), "PF_VAL_BC_DATA")
    # Multi plate valid bc counts?
    if spipe.num_sample_plates() > 1:
        df, n = merge_multi_csvs_df(
            spipe, "PF_VAL_2BC_DATA", None, axis=1, combine=True
        )
        df = aggregate_col_df(df, "count", mode="sum")
        spipe.write_df(df.astype(int), "PF_VAL_2BC_DATA")

    # Stat csv files; Pipeline stats
    df, n = merge_multi_csvs_df(spipe, "PF_STAT_PIPE", None, axis=1, combine=True)
    # Weights and lists of values that get sum or weighted mean
    w_cols = [x for x in df.columns if x != "value"]
    weights = list(df.loc["number_of_reads"][w_cols])
    sum_lis, mean_lis = pb_stats.get_stat_sum_mean_keys(df.index, combine=True)
    df = aggregate_col_df(df, "value", mode="sum", row_lis=sum_lis, na=0)
    df = aggregate_col_df(df, "value", mode="mean", weights=weights, row_lis=mean_lis)
    # Dump all lines to log; n_rows=0
    spipe.write_df(df, "PF_STAT_PIPE", n_rows=0)

    # Sequence stats
    df, n = merge_multi_csvs_df(spipe, "PF_STAT_SEQ", None, axis=1, combine=True)
    sum_lis, mean_lis = pb_stats.get_stat_sum_mean_keys(df.index, combine=True)
    df = aggregate_col_df(df, "value", mode="sum", row_lis=sum_lis, na=0)
    df = aggregate_col_df(df, "value", mode="mean", weights=weights, row_lis=mean_lis)
    spipe.write_df(df, "PF_STAT_SEQ")


# -------------- Combine mode, per each sample ------------------------------
#
def combine_single_dge_output(spipe, samp):
    """Combine sublibrary data for one given sample

    Return status
    """
    assert spipe.get_mode("comb"), "combine_single_dge_output for comb mode"

    # Target enrichment
    targeted = spipe.is_targeted()
    is_tcr = spipe.is_tcr()

    # Data files; DGE only; TCR not here
    if not is_tcr:
        # List of DGE dir keys; First is primary (e.g. filtered for analysis, report)
        dir_key_list = const.use_case_out_dir_key(spipe, as_list=True)
        for i, dir_key in enumerate(dir_key_list):
            # Non-primary combine is conditional
            if (i > 0) and (not spipe.get_par_val("comb_unfilt_dge", as_bool=True)):
                spipe.report_run_story2(f"Not combining DGE files; {dir_key}")
                continue

            spipe.report_run_story2(f"Combining DGE files; {dir_key}")
            if not merge_multi_dge_files(spipe, samp, dir_key, combine=True):
                spipe.set_problem(
                    f"Failed to combine DGE files for sample {samp.get_name()}"
                )
                return False
        # Expressed genes files?
        if spipe.get_par_val("dge_save_express_gene_file", as_bool=True):
            merge_multi_exgene_files(spipe, samp, combine=True)

    # Combine sample specific csv files (i.e. report dir)
    spipe.report_run_story2("Merging combine-mode csv files")

    # Analysis_summary files; Some fields sum, some mean...
    #   Weights are for weighted average (number of reads) for mean aggregations
    stat_df, weights = merge_multi_stats_df(spipe, samp, combine=True)
    asum_keys, _ = pb_stats.get_asum_out_keys(stat_df.index)

    # Dict with values from new merged column
    merge_col = stat_df.columns[0]
    stat_dict = stat_df.loc[:, merge_col].to_dict()

    # TCR metrics or non-TCR per-sample report files
    if is_tcr:
        # This only updates TCR files; Values in all-stats not fixed here
        merge_multi_tcr_metrics(spipe, samp)

    else:
        # Per barcode ('cell') counts
        merge_multi_bc_counts(spipe, samp, combine=True)
        # Subsample stuff
        merge_multi_ss_csvs(spipe, samp, combine=True)
        # Counts per well, etc
        merge_multi_count_files(spipe, samp, weights, combine=True)
        # Focal / CRISPR specific
        if spipe.is_focal_crispr():
            merge_multi_focal_files(spipe, samp)

        # Some stats get recalculated de novo
        spipe.report_run_story2("Recalculating some merged stats")
        up_dict = pb_stats.updated_stats_dict(spipe, samp, stat_dict, targeted=targeted)
        for k, v in up_dict.items():
            stat_df.loc[k, "combined"] = v

    pb_stats.save_analysis_summary(spipe, samp, stat_df, col="combined")

    return True


def merge_multi_dge_files(spipe, samp, dir_key, combine=True):
    """Combine multiple DGE files for sublibraries or meta sample parts

    combine = flag to combine sublibraries -vs- meta sample parts

    Return status
    """
    targeted = spipe.is_targeted()

    # DGE file keys and output filenames
    dge_fkeys = const.out_filepaths(
        spipe, samp, dir_key=dir_key, as_key=True, targeted=targeted
    )
    print(">>> merge_multi_dge_files", dir_key, dge_fkeys["mtx"])
    mtx_ofname = spipe.filepath(dge_fkeys["mtx"], samp)
    cell_ofname = spipe.filepath(dge_fkeys["cell"], samp)
    gene_ofname = spipe.filepath(dge_fkeys["gene"], samp)
    # Tallies
    total_cells = total_data = n_gene = 0

    # Combining sublibrary data
    if combine:
        part_list = spipe.get_sublibs()
        spipe.report_run_story2(f"DGE files from potentially {len(part_list)} sublibs")
    # Combining meta sample (non-meta sample) part data
    else:
        part_list = samp.get_meta_parts(unwrap=True)
        spipe.report_run_story2(f"DGE files from potentially {len(part_list)} parts")

    # Each sublib or sample part
    n_samps = 0
    n_good = 0
    for i, part in enumerate(part_list):
        if combine:
            # If combine sublib doesn't have current sample, ignore
            if samp.num_wells(sublib=i) < 1:
                spipe.report_run_story2(
                    f"No sample {samp.get_name()} in sublib {part.get_name()}; Ignoring"
                )
                continue
            # top_dir = path for sublib, sample = given sample
            s_top_dir = part.get_path()
            samp_part = samp
        # Meta sample; top_dir = default, sample = sample part
        else:
            s_top_dir = ""
            samp_part = part
        n_samps += 1

        # Check if output exists
        cell_fpath = spipe.filepath(dge_fkeys["cell"], samp_part, top_dir=s_top_dir)
        if not utils.check_infile(cell_fpath, verb=False, toxic=False):
            story = f" ({s_top_dir})" if s_top_dir else ""
            spipe.report_warning(
                f"No cells (DGE files) for {samp_part.get_name()}{story}"
            )
            continue
        n_good += 1

        # Cell data
        cell_df = spipe.read_csv(dge_fkeys["cell"], samp=samp_part, top_dir=s_top_dir)
        # Add sublib prex to cell barcodes if combine; Barcode is dataframe index.
        if combine:
            cell_df = bcutils.add_suf_to_df_bci(cell_df, i + 1, index=True)
        cell_df.index.name = "bc_wells"

        # Load matrix data
        mtx_fname = spipe.filepath(dge_fkeys["mtx"], samp_part, top_dir=s_top_dir)
        head_lines, data_lines, story = read_mtx_file_lines(mtx_fname)
        if not data_lines:
            spipe.set_problem(story)
            spipe.set_problem(f"Failed to load DGE mtx from {mtx_fname}")
            return False
        spipe.report_run_story2(story)

        # Remove any zero-value lines; e.g. at end of files
        data_lines = mtx_strip_zero_data_lines(data_lines)

        # Get n_cell and n_gene from mtx header
        # i.e. %Rows=cells (468), Cols=genes (112625)
        n_cell, n_gene = head_row_cell_col_gene(head_lines)
        n_data = len(data_lines)
        assert n_cell == len(
            cell_df
        ), f"Cell counts differ; mtx {n_cell}, cell {len(cell_df)}"

        # First sublib / part is special
        if n_good == 1:
            # Write matrix header once here; Will update numbers after loop
            write_mtx_header(mtx_ofname, n_cell, n_gene, n_data)
            # Write cells new, with header
            utils.write_df(cell_df, cell_ofname, header=True, index=True)
        else:
            # Increment cell (row) indexes in matrix data lines
            new_data_lines = []
            for line in data_lines:
                r, c, v = line.split()
                nline = f"{int(r)+total_cells} {c} {v}"
                new_data_lines.append(nline)
            data_lines = new_data_lines
            # Append cells, no header
            utils.write_df(cell_df, cell_ofname, header=False, index=True, mode="a")

        # Append matrix data
        with open(mtx_ofname, mode="a") as OUTFILE:
            print("\n".join(data_lines), file=OUTFILE)

        total_cells += n_cell
        total_data += n_data
        spipe.report_run_story2(f"Processed cells {len(cell_df)}, {total_cells} total")

    # Only save anything if any good inputs
    if n_good:
        # Finish matrix; Add final capping data line, then update header
        with open(mtx_ofname, mode="a") as OUTFILE:
            print(f"{total_cells} {n_gene} {0}", file=OUTFILE)
            total_data += 1

        write_mtx_header(mtx_ofname, total_cells, n_gene, total_data, update=True)
        spipe.report_run_story(
            f"Combined DGE matrix {total_data} elements saved to {mtx_ofname}"
        )
        # Copy genes file (all_genes or target_genes) to DGE output dir
        if spipe.get_par_val("dge_copy_dge_genes_file", as_bool=True):
            src = const.geno_filepaths(spipe, samp, targeted=targeted)["genes"]
            story = utils.copy_file(src, gene_ofname)
            spipe.report_run_story2(story)
        else:
            spipe.report_run_story2(
                f"Not copying {str} to DGE dir for {samp.get_name()}"
            )

    else:
        story = f"No DGE files combined for {samp_part.get_name()} ({dir_key})"
        spipe.report_warning(story)

    return True


def mtx_strip_zero_data_lines(data_lines):
    """Strip lines with zero value

    Return list of matrix data lines
    """
    new_lines = []
    for line in data_lines:
        r, c, v = line.split()
        if int(v) > 0:
            new_lines.append(line)
    return new_lines


def read_mtx_file_lines(fname, header="%", hdata=1):
    """Read from mtx (sparse matrix) file

    header = header line starting char
    hdata = how many data lines to include in header

    Return tuple of lists and str (header, data, story)
    """
    head_lines = []
    data_lines = []
    n = 0
    with open(fname) as ifile:
        while True:
            line = ifile.readline()
            if not line:
                break
            if line.startswith(header):
                head_lines.append(line.strip())
            else:
                n += 1
                if n <= hdata:
                    head_lines.append(line.strip())
                else:
                    data_lines.append(line.strip())
        # Check that header row and col comment
        if head_row_cell_col_gene(head_lines):
            story = f"Loaded {len(data_lines)} data from {fname}"
        else:
            story = "No 'rows=cells' 'cols=genes' listed in header"
            head_lines = data_lines = []

    return head_lines, data_lines, story


# Match something like:
#   %Rows=cells (468), Cols=genes (112625)
HEAD_R_C_REGEX = re.compile(r"%Rows=cells \((\d+)\), Cols=genes \((\d+)\)\s*")


#
def head_row_cell_col_gene(head):
    """Look for row=cell and col=gene annotation in header line(s)

    Return tuple(max_cell, max_gene)
    """
    # If list, check each element
    if isinstance(head, list):
        for line in head:
            tup = head_row_cell_col_gene(line)
            if tup:
                return tup
    else:
        # Should be string
        assert isinstance(head, str)
        mat = HEAD_R_C_REGEX.match(head)
        if mat:
            return (int(mat[1]), int(mat[2]))
    return None


# Fixed width lines so can replace / update later
MTX_HEAD_LINEWIDTH = 80


#
def write_mtx_header(fname, max_cell, max_gene, n_data, update=False):
    """Write (sparse matrix) mtx file header block

    Returns nothing
    """
    # Open file to update (i.e. overwrite) or write
    mode = "r+" if update else "w"
    with open(fname, mode=mode) as OUTFILE:
        # Header lines
        head1 = "%%MatrixMarket matrix coordinate integer general"
        head2 = f"%Rows=cells ({max_cell}), Cols=genes ({max_gene})"
        head3 = f"{max_cell} {max_gene} {n_data}"

        # Pad to fixed length and write
        print(utils.pad_line(head1, plen=MTX_HEAD_LINEWIDTH), file=OUTFILE)
        print(utils.pad_line(head2, plen=MTX_HEAD_LINEWIDTH), file=OUTFILE)
        print(utils.pad_line(head3, plen=MTX_HEAD_LINEWIDTH), file=OUTFILE)


def merge_multi_exgene_files(spipe, samp, combine=True):
    """Combine expressed gene files (if / when they exist)

    combine = flag to combine sublibrary data vs meta sample parts

    Return status
    """
    df_list = []
    num_df = 0
    id_names = {}
    id_genomes = {}

    # Combine sublibrary data
    if combine:
        # Load tables, collecting gene name and genome info across all
        for slib in samp.get_sublibs():
            # Load expressed gene table; If it exists
            s_top_dir = slib.get_path()
            df = spipe.check_read_csv("SFR_EXPRESS_GENE", samp=samp, top_dir=s_top_dir)
            if isinstance(df, pd.DataFrame):
                id_names.update(df["gene_name"].to_dict())
                id_genomes.update(df["genome"].to_dict())
                num_df += 1
            # Save tuple sublib and dataframe (or nothing)
            # (Save regardless so index in sublib list is correct)
            df_list.append((slib, df))
    # Combine meta sample (non-meta sample) part data
    else:
        for samp_part in samp.get_meta_parts(unwrap=True):
            df = spipe.check_read_csv("SFR_EXPRESS_GENE", samp=samp_part)
            if isinstance(df, pd.DataFrame):
                id_names.update(df["gene_name"].to_dict())
                id_genomes.update(df["genome"].to_dict())
                num_df += 1
            df_list.append((samp_part, df))

    # If nothing collected, all done
    if num_df < 1:
        return True

    # New dataframe with all combined gene info
    g_data = {"gene_name": id_names, "genome": id_genomes}
    new_df = pd.DataFrame.from_dict(g_data, orient="columns")
    # Add counts for each sublib / part; Keep list of cols to combine
    c_cols = []
    for i, (part, df) in enumerate(df_list):
        if isinstance(df, pd.DataFrame):
            if combine:
                name = f"s{i+1}"
            else:
                name = part.get_name()
            # Name col as sublib number (not name)
            new_col = f"cell_count__{name}"
            new_df[new_col] = df["cell_count"]
            c_cols.append(new_col)

    # Count = sum of parts; Also set count cols to int
    new_df["cell_count"] = new_df[c_cols].sum(axis=1)
    int_dict = {"cell_count": int}
    for c in c_cols:
        int_dict[c] = int
    new_df = new_df.fillna(0).astype(int_dict)

    # Put count up front in cols and sort index
    o_cols = ["gene_name", "genome", "cell_count"]
    o_cols += [c for c in c_cols if c != "cell_count"]
    new_df = new_df[o_cols].sort_index()
    new_df.index.name = "gene_id"
    spipe.write_df(new_df, "SFR_EXPRESS_GENE", samp=samp, index=True)

    return True


def merge_multi_csvs_df(
    spipe,
    fkey,
    samp,
    axis=0,
    add_n=True,
    addn_col=None,
    suf_col=True,
    name_col=False,
    combine=True,
):
    """Load and merge csv files indicated by fkey for all sublibraries into new dataframe

    fkey = file keyword
    samp = sample (or None)
    axis = how to combine; 0 = rows, 1 = cols
    add_n = flag to append sublibrary count (i.e. n) to index or column in axis=0=row mode
    addn_col = which column to add n to; If None and add_n, target index
    suf_col = flag to add sublib or metasample name as col suffix in axis=1=col mode
    name_col = flag to use sublib or metasample name for col name in axis=1=col mode
    combine = flag to use sublibraries else meta sample parts

    Return tuple (dataframe, number merged)
    """
    # Combining sublibrary data
    if combine:
        part_list = spipe.get_sublibs()
    # Combine meta sample (non-meta sample) part data
    else:
        part_list = samp.get_meta_parts(unwrap=True)
        # No addition of sublib number
        add_n = False

    df_lis = []
    # Collect explicit number for each part
    num_lis = []
    n = 0
    # Each sublib / meta sample part
    for i, part in enumerate(part_list):
        n += 1
        # If combine and and sublib doesn't have current sample, ignore
        if combine:
            # samp may be None (e.g. global files are not sample dependent)
            if samp and samp.num_wells(sublib=i) < 1:
                spipe.report_run_story2(
                    f"No sample {samp.get_name()} in {part.get_name()}; ignoring"
                )
                continue
            # top_dir = sublib and sample is sample itself
            top_dir = part.get_path()
            samp_part = samp
        # Sample part; default top_dir = current
        else:
            top_dir = ""
            samp_part = part

        # May not have?
        df = spipe.check_read_csv(fkey, samp=samp_part, top_dir=top_dir)
        if df is None:
            csv_fpath = spipe.filepath(fkey, samp_part, top_dir=top_dir, mkdir=False)
            name = part.get_name() if part else "<none>?"
            story = f"Combine for {name} missing file: {csv_fpath}"
            spipe.report_run_story2(story)
            continue

        # If merging col-wise, may need to clean up / modify cols
        if axis == 1:
            # Meta samples have extra columns; Only want first one, named as sample
            # Have to check have samp as this could be None (e.g. global files)
            if samp and samp.is_meta():
                df = pd.DataFrame(df.iloc[:, 0])
                df.columns = [samp.get_name()]

            # add suffix or replace with part name
            if suf_col or name_col:
                new_cols = []
                for col in df.columns:
                    if name_col:
                        new_cols.append(part.get_name())
                    else:
                        new_cols.append(f"{col}__{part.get_name()}")
                df.columns = new_cols
        df_lis.append(df)
        num_lis.append(n)

    # Merge what's collected
    if df_lis:
        if (axis == 0) and add_n:
            # Add (sublib) number to indexes
            for i, df in enumerate(df_lis):
                num = num_lis[i]
                if addn_col in df.columns:
                    df[addn_col] = df[addn_col].astype(str) + f"__s{num}"
                else:
                    df.index = df.index.astype(str) + f"__s{num}"
        df = pd.concat(df_lis, axis=axis)
    else:
        df = None

    # return df, df_lis
    return df, len(df_lis)


def merge_multi_stats_df(spipe, samp, combine=False):
    """Merge analysis_summary and allstats csv files for sublibraries / meta sample parts

    all_stats = flag to merge all_stats file, else analysis_summary
    combine = flag for combine mode

    Returns tuple (merged_df, list of weights per part)
    """
    # Combine columns for each (part for sample)
    df, n = merge_multi_csvs_df(spipe, "SFR_ALLSTATS", samp, axis=1, combine=combine)
    # Set any missing cell counts to zero
    df = replace_df_nan_vals(df, row_key="number_of_cells", na=0, substring=True)
    # Weights as number of reads for each part
    weights = list(df.loc["number_of_reads"])
    # List of cols that should be merged via sum or via (weighted) mean
    sum_lis, mean_lis = pb_stats.get_stat_sum_mean_keys(df.index, combine=combine)
    # New 'value' column gets sum, weighted average agregate from diff row subsets
    df = aggregate_col_df(df, "value", mode="sum", row_lis=sum_lis)
    df = aggregate_col_df(df, "value", mode="mean", weights=weights, row_lis=mean_lis)
    df = clean_merged_df_col_names(df)
    return df, weights


def clean_merged_df_col_names(
    df, first_col="combined", kill_str="value__", as_float=True
):
    """Clean up column names in merged dataframe

    Return dataframe
    """
    # Clean up column names; Replace value__<part> with <part>
    new_cols = []
    for col in df.columns:
        new_cols.append(col.replace(kill_str, ""))
    # Repace first col label (e.g. 'value' with 'combined')
    new_cols[0] = first_col
    df.columns = new_cols
    # Cast first col
    if as_float:
        df[first_col] = df[first_col].astype(float)
    return df


def merge_multi_bc_counts(spipe, samp, combine=True):
    """Merge per-barcode count csv files all sublibraries or meta sample parts

    combine = flag to use sublibraries else meta sample parts

    output written to files

    Returns nothing
    """
    # Combining sublibs = add to bc index, meta sample = no change
    add_n = True if combine else False

    # Total counts
    fkey = "SFR_TSCP_CT"
    df, n = merge_multi_csvs_df(spipe, fkey, samp, axis=0, add_n=add_n, combine=combine)
    if n:
        spipe.write_df(df["count"].astype(int), fkey, samp=samp)
    else:
        story = f"No tscp count files combined for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)

    # Per species counts
    fkey = "SFR_SPEC_TSCP_CT"
    df, n = merge_multi_csvs_df(spipe, fkey, samp, axis=0, add_n=add_n, combine=combine)
    if n:
        spipe.write_df(df.astype(int), fkey, samp=samp)
    else:
        story = f"No species tscp count files combined for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)

    fkey = "SFR_SPEC_READ_CT"
    df, n = merge_multi_csvs_df(spipe, fkey, samp, axis=0, add_n=add_n, combine=combine)
    if n:
        spipe.write_df(df.astype(int), fkey, samp=samp)
    else:
        story = f"No species read count files combined for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)

    # Enrichment if targeted and not focal_bc
    if spipe.is_targeted() and (not spipe.is_focal_crispr()):
        fkey = "SFR_ENRICHMENT"
        df, n = merge_multi_csvs_df(
            spipe, fkey, samp, axis=0, add_n=add_n, combine=combine
        )
        if n:
            spipe.write_df(df, fkey, samp=samp)
        else:
            story = f"No enrichment files combined for {samp.get_name()} ({fkey})"
            spipe.report_warning(story)


def merge_multi_ss_csvs(spipe, samp, combine=True):
    """Merge subsample csv files for all sublibraries, write output csv

    output written to files

    Returns nothing
    """
    # Simply combine by appending columns; No aggregate col(s) added
    #   Result can have many 'na' values for subsample values not in sublibs
    fkey = "SFR_SS_GENE_CT"
    df, n = merge_multi_csvs_df(spipe, fkey, samp, axis=1, combine=True)
    if n:
        spipe.write_df(df, fkey, samp=samp)
    else:
        story = f"No gene subsample files combined for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)

    fkey = "SFR_SS_TSCP_CT"
    df, n = merge_multi_csvs_df(spipe, fkey, samp, axis=1, combine=True)
    if n:
        spipe.write_df(df, fkey, samp=samp)
    else:
        story = f"No tscp subsample files combined for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)


def merge_multi_count_files(spipe, samp, weights, combine=True):
    """
    Merge multiple count files for sublibrary or meta sample parts

    weights = weights for each part (for weighted mean)
    combine = flag to use sublibraries else meta sample parts

    output written to files

    Returns nothing
    """
    # cutoff info = mean
    fkey = "SFR_TSCP_CUTOFF"
    df, n = merge_multi_csvs_df(spipe, fkey, samp, axis=1, combine=combine)
    if n:
        df.fillna(0, inplace=True)
        df = aggregate_col_df(df, "value", mode="mean", weights=weights)
        spipe.write_df(df, fkey, samp=samp)
    else:
        story = f"No tscp cutoff files combined for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)

    # MapReadTscpGene counts = sum
    # Load and col-combine all sublibary csvs; Aggreate new sum column
    fkey = "SFR_SAMP_MRTG_CT"
    df, n = merge_multi_csvs_df(spipe, fkey, samp, axis=1, combine=combine)
    if n:
        df.fillna(0, inplace=True)
        df = aggregate_col_df(df, "count", mode="sum")
        spipe.write_df(df.astype(int), fkey, samp=samp)
    else:
        story = f"No mrtg count files combined for {samp.get_name()} ({fkey})"
        spipe.report_warning(story)

    # For per-well tscp and cell counts, load merged cell_df then process
    dge_fkeys = const.out_filepaths(spipe, samp, as_key=True)
    cell_df = spipe.check_read_csv(dge_fkeys["cell"], samp=samp)
    if cell_df is not None:
        targeted = spipe.is_targeted()
        dge.write_samp_per_well_count_files(spipe, samp, cell_df, targeted=targeted)
    else:
        story = f"No cell files combined for {samp.get_name()} ({dge_fkeys['cell']})"
        spipe.report_warning(story)


def merge_multi_tcr_metrics(spipe, samp):
    """Merge TCR metric files for sample

    Returns nothing
    """
    # Always have unfiltered; Only filtered if parent
    fkey_list = ["SF_TCR_UF_STAT"]
    if spipe.get_use_case("tcr_parent"):
        fkey_list.append("SF_TCR_F_STAT")

    for fkey in fkey_list:
        # Get existing combine output and merged collection from sublibs
        tcr_df = spipe.read_csv(fkey, samp=samp)
        parts_df, num = merge_multi_csvs_df(spipe, fkey, samp, combine=True, axis=1)
        # Merge existing output with sublibs, clean names (e.g. col 1 = 'combined')
        df = pd.concat([tcr_df, parts_df], axis=1)
        df = clean_merged_df_col_names(df)

        spipe.write_df(df, fkey, samp=samp)


def merge_multi_focal_files(spipe, samp, combine=True):
    """Merge focal / CRSIPR files for sample

    Returns nothing
    """
    focal = spipe.is_focal_crispr(focal=True)

    # guide counts
    fkey = "SFR_FB_CELL_CT" if focal else "SFR_CRSP_CELL_CT"
    parts_df, num = merge_multi_csvs_df(spipe, fkey, samp, combine=combine, axis=1)
    cnt_df = aggregate_col_df(parts_df, "combined", mode="sum").astype(int)
    cnt_df = clean_merged_df_col_names(cnt_df, kill_str="count__", as_float=False)
    spipe.write_df(cnt_df, fkey, samp=samp)

    # Target to cell mapping; All strings, so simply concat and sort
    fkey = "SFR_FB_CELL_MAP_TAB" if focal else "SFR_CRSP_CELL_MAP_TAB"
    parts_df, num = merge_multi_csvs_df(
        spipe, fkey, samp, combine=combine, axis=0, addn_col="bc_wells"
    )
    map_df = parts_df.sort_index()
    spipe.write_df(map_df, fkey, samp=samp)

    # Also save as json guide-to-cell dict
    map_dict = crispr.g2cell_dict_from_df(map_df)
    fkey = (
        "SFR_FB_CELL_MAP_DIC"
        if spipe.is_focal_crispr(focal=True)
        else "SFR_CRSP_CELL_MAP_DIC"
    )
    fname = spipe.filepath(fkey, samp)
    utils.write_json(map_dict, fname)
    spipe.report_run_story(f"Wrote guide-to-cell dictionary {fname}")


def aggregate_col_df(df, agg_col, mode="sum", weights=[], row_lis=None, na="NaN"):
    """Add / replace a column with aggregate of dataframe values (sum, mean)

    agg_col = aggregate column(s); If exists, replace values, else add new
    mode = how to combine; append, sum, mean
    weights = weights for weighted average (mode mean)
    row_lis = subset of rows to modify
    na = value to replace NaN if not None

    Return dataframe
    """
    col_lis = list(df.columns)
    # If aggregate col is new, create, init to na
    if agg_col not in col_lis:
        df[agg_col] = "NA"
        # agg_col will be first col
        col_lis = [agg_col] + col_lis

    # List of cols to aggregate over
    agg_lis = [c for c in col_lis if c != agg_col]

    # If no row subset, all rows
    if row_lis is None:
        row_lis = list(df.index)

    # One row at a time
    for row_lab in row_lis:
        # Try to get values as number list (may be strings?)
        try:
            num_list = df.loc[row_lab, agg_lis].astype(float).values
        # Else value = first item
        except Exception:
            num_list = None
            v = df.loc[row_lab, agg_lis].values[0]

        # If have numbers, combine
        if num_list is not None:
            # Clean up any NaN
            new_list, new_weights = clean_list_nan_vals(
                num_list, weights=weights, na=na
            )

            if len(new_list) == 0:
                v = na
            elif mode == "sum":
                v = np.sum(new_list)
            elif mode == "mean":
                if len(new_weights) > 0:
                    v = np.average(new_list, weights=new_weights)
                else:
                    v = np.average(new_list)
            else:
                story = f"Bad mode in aggregate_col_df: '{mode}'"
                raise ValueError(story)

        df.at[row_lab, agg_col] = v

    # Reorder so agg column is first
    df = df[col_lis]

    return df


def clean_list_nan_vals(num_list, weights=[], drop_na=True, na="NaN"):
    """Remove / replace NaN values in list, adjusting associated weight list too

    Return tuple of lists
    """
    new_list = []
    new_weights = []

    for i, v in enumerate(num_list):
        if np.isnan(v):
            if drop_na:
                continue
            v = na
        if not np.isnan(v):
            new_list.append(v)
            if len(weights):
                new_weights.append(weights[i])

    return new_list, new_weights


def replace_df_nan_vals(df, row_key="number_of_cells", na=0, substring=True):
    """Replace NaN values in dataframe with given value

    Return dataframe
    """
    for row in df.index:
        # Qualify row for replacement
        if row_key:
            if substring:
                keep = row_key in row
            else:
                keep = row_key == row
            if not keep:
                continue
        df.loc[row] = df.loc[row].fillna(na)
    return df


def combine_clean_up(spipe):
    """Clean up any combine-specific things"""
    return True
