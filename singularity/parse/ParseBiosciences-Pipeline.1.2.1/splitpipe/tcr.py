#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import os
import re
import sys
import traceback
import pandas as pd

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import combine as comb
    import pb_stats
    import utils
else:
    from splitpipe import combine as comb
    from splitpipe import pb_stats
    from splitpipe import utils


# ---------------------------------------------------------------------------
def get_tcr_call_igblast_paths():
    """Get paths to TCR call script and IgBlast

    Return tuple of paths
    """
    # Get top level path
    try:
        path = os.path.dirname(__file__)
    except Exception:
        path = "."

    # Specific locations
    call_path = f"{path}/scripts/run_tcr_AK.py"
    igb_path = f"{path}/TCR/IgBlast/"
    db_path = f"{path}/TCR/IgBlast/tr_database/"
    return call_path, igb_path, db_path


def have_tcr_ins(spipe, verb=True):
    """Check if TCR input files exist

    Return tuple status, list of files
    """
    # Combine mode only needs sublibs
    if spipe.is_combine():
        return True, []

    ok = False
    # Check IgBlast path and fasta
    check_files = []
    _, igblast_dir, db_path = get_tcr_call_igblast_paths()
    infastq = spipe.filepath("PF_FASTQ_BC", None)

    check_files.append(infastq)
    check_files.append(igblast_dir)

    # If list is good, we have inputs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    else:
        story = "Bad TCR input file(s) / Installed database(s)"
        spipe.set_problem(story)

    # Check specific genome files
    if ok:
        if not have_igblast_dbs(spipe, verb=True):
            ok = False
            check_files.append("Genome-specific-db-files")

    return ok, check_files


def have_tcr_outs(spipe):
    """Check if TCR output files exist

    Return True/False, list of files
    """
    ok = False
    check_files = []
    for samp in spipe.get_samples():
        # Only check unfiltered stat file
        check_files.append(spipe.filepath("SF_TCR_UF_STAT", samp))

    # If list is good, we have outs
    if utils.check_infile(check_files, verb=False):
        ok = True
    return (ok, check_files)


def run_tcr(spipe):
    """Handle main TCR processing

    Return status
    """
    extra = "combine" if spipe.is_combine() else "analysis"
    step_name = f"tcr ({extra})"
    spipe.report_proc_step(step_name, set_mode=False)
    ok = True
    # Input status and list
    i_ok, i_list = have_tcr_ins(spipe)
    o_ok, o_list = have_tcr_outs(spipe)
    # Run only if don't already have outputs, or if fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for tcr: {bad_list}"
            spipe.set_problem(story)
        else:
            run_tcr_analysis(spipe)

        # If all good, aggregate summary files
        if spipe.no_problems():
            pb_stats.aggregate_asum_csv(spipe)

    else:
        spipe.report_run_story("Skipping tcr and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")

    ok = spipe.no_problems()
    spipe.report_proc_step(step_name, status=ok)
    return ok


def run_tcr_analysis(spipe):
    """Run tcr analysis script, check outputs, process stats

    Return status
    """
    # Catch possible problems
    try:
        ok = run_tcr_analysis_script(spipe)
    except Exception as e:
        ok = False
        # Dump traceback, stdout and log
        traceback.print_exc(file=sys.stdout)
        sys.stdout.flush()
        traceback.print_exc(file=spipe.get_log_fh())
        story = f"run_tcr_analysis_script: {e}"
        spipe.set_problem(story)

    if ok:
        ok = check_tcr_outputs(spipe)
    if ok:
        ok = finish_tcr_stats(spipe)
    return ok


def have_igblast_dbs(spipe, verb=True):
    """Check if specified genome database files are found

    Return bool
    """
    have = False
    db_names = get_igblast_db_list(trunc_at="_")
    tcr_genome = get_tcr_genome_name(spipe)
    if verb:
        story = f"Checking for {tcr_genome} in DB list {db_names}"
        spipe.report_run_story(story)

    have = bool(tcr_genome in db_names)
    return have


def get_tcr_genome_name(spipe):
    """Get TCR genome name based on params

    Return str
    """
    tcr_genome = spipe.get_par_val("tcr_genome", as_str=True)
    if spipe.get_par_val("tcr_use_imgt_db", as_bool=True):
        tcr_genome = "imgt_" + tcr_genome + "_tcr"
    else:
        tcr_genome = tcr_genome + "_tcr"
    return tcr_genome


def get_igblast_db_list(db_ext=".nsq", trunc_at=""):
    """Get list of installed database names

    db_ext = DB file extension to collect
    trunc_at = char to truncate names with (e.g. with '_' 'human_V' >--> 'human')

    Return sorted list of names
    """
    _, _, db_path = get_tcr_call_igblast_paths()
    name_set = set()

    # Tally ".nseq" files in db path
    fnames = utils.fnames_for_path(path=db_path, ext=db_ext)
    for f in fnames:
        # Only name without extension
        name = os.path.basename(f).replace(db_ext, "")
        # Truncating?
        if trunc_at:
            # name = name.split(trunc_at)[0]
            name = trunc_at.join(name.split(trunc_at)[:-1])
        name_set.add(name)

    return sorted(name_set)


def write_sample_def_list(spipe, fname):
    """Write out list of sample defs for TCR script

    Returns number of samples written
    """
    n = 0
    with open(fname, "w") as OUTFILE:
        for samp in spipe.get_samples():
            name = samp.get_name()
            wind_str = ",".join(samp.get_wind_list())
            print(f"{name}\t{wind_str}", file=OUTFILE)
            n += 1
    return n


def write_sublib_list(spipe, fname):
    """Write out list of sublib paths for TCR script

    Returns number of sublibs written
    """
    n = 0
    with open(fname, "w") as OUTFILE:
        for slib in spipe.get_sublibs():
            path = slib.get_path()
            print(path, file=OUTFILE)
            n += 1
    return n


def tcr_parent_tscp_folder_samp(spipe):
    """Find the 'transcriptome folder' sample for the tcr parent case

    Return sample obj
    """
    best = None
    # First look for all-sample
    for samp in spipe.get_samples(non_allwell=False, simple=False):
        best = samp
        break

    # Next look for all-well
    if not best:
        for samp in spipe.get_samples(non_allwell=False, meta=False):
            best = samp
            break

    # Finally, look for sample with the most wells
    if not best:
        max_wells = 0
        for samp in spipe.get_samples(allwell=False):
            if samp.num_wells() > max_wells:
                best = samp
                max_wells = samp.num_wells()
    return best


def run_tcr_analysis_script(spipe):
    """Run tcr analysis script"""
    spipe.report_run_story("Preparing command line for TCR")

    # Paths to script and IgBlast
    script_path, igblast_dir, _ = get_tcr_call_igblast_paths()
    tcr_genome = get_tcr_genome_name(spipe)

    # Get script version info
    ok, path, script_ver, story = utils.get_exec_path_and_ver(script_path, "-V")
    spipe.report_run_story(f"Version info: TCR-script {path} {script_ver}")

    output_dir = spipe.filepath("DIR_TOP", None)
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    trim = spipe.get_par_val("tcr_trim", as_int=True)
    # Cook up input command line
    tcr_command = f"{script_path} \
        --wkdir {output_dir} \
        --igblastdir {igblast_dir} \
        --genome {tcr_genome} \
        --nthreads {nthreads} \
        --trim {trim} \
        "

    # List of samples as <name> <well,index,list>
    slist_file = spipe.filepath("TMP_TEMP1_FILE", None)
    n = write_sample_def_list(spipe, slist_file)
    spipe.report_run_story(f"Specified {n} samples in list file")
    tcr_command += f" --wfile {slist_file}"

    # Combine mode?
    if spipe.is_combine():
        slist_file = spipe.filepath("TMP_TEMP2_FILE", None)
        n = write_sublib_list(spipe, slist_file)
        spipe.report_run_story(f"Specified {n} sublibs in list file")
        tcr_command += f" --sublibraries {slist_file}"
        tcr_command += f" --output_dir {output_dir}"
    # Normal non-combine
    else:
        # Parent?
        parent_dir = spipe.get_par_val("parent_dir", as_path=True)
        if parent_dir:
            # Actually need path to cell file for sample in parent
            psamp = tcr_parent_tscp_folder_samp(spipe)
            fpath = spipe.filepath("SF_NF_CELL", psamp, top_dir=parent_dir)
            tcr_command += f" --cell_file {fpath}"

    # Clean to single spaces
    command = re.sub(" +", " ", tcr_command)
    ok, callout, ex_story = utils.call_exec(command, verb=True)
    # Story to log only (call verb=True already prints outcome)
    spipe.report_run_story(ex_story, to_log=True)

    # Output saved to file
    ver = spipe.get_version()
    fpath = spipe.filepath("PF_TCR_SUBPROC_OUT", None)
    utils.file_header(
        ofile=fpath, story=f"TCR subproc {script_path} {script_ver}", ver=ver
    )
    with open(fpath, "a") as OUTFILE:
        com_str = utils.wrap_exe_com_args(command)
        print("\n# Call:", file=OUTFILE)
        print(com_str, file=OUTFILE)
        print("\n# Output:", file=OUTFILE)
        print(callout, file=OUTFILE)
        spipe.report_run_story(f"TCR subproc output saved to {fpath}")

    if not ok:
        spipe.set_problem("Subprocess failed: " + command)
        spipe.set_problem(ex_story)

    return ok


def check_tcr_outputs(spipe):
    """Process / clean tcr script outputs

    Return number output files found
    """
    #update_tcr_outputs(spipe)

    spipe.report_run_story("Checking TCR output files")
    tot_outs = 0
    for samp in spipe.get_samples():
        n, missing = check_samp_tcr_outputs(spipe, samp, filt=False)
        spipe.report_run_story2(
            f"Sample {samp.get_name()} found {n} TCR output files (unfiltered)"
        )
        tot_outs += n
        if missing:
            story = f"Failed to get {len(missing)} expected TCR outputs; {missing}"
            spipe.set_problem(story)

        if spipe.get_use_case("tcr_parent"):
            n, missing = check_samp_tcr_outputs(spipe, samp, filt=True)
            spipe.report_run_story2(
                f"Sample {samp.get_name()} found {n} TCR output files (filtered)"
            )
            tot_outs += n
            if missing:
                story = f"Failed to get {len(missing)} expected TCR outputs; {missing}"
                spipe.set_problem(story)
                
    return tot_outs


def check_samp_tcr_outputs(spipe, samp, filt=False):
    """Check TCR output files exist for given sample

    Return tuple (number found files, list of missing)
    """
    key_list = ['SF_TCR_F_STAT',
        'SF_TCR_F_AIRR',
        'SF_TCR_F_BC_REPORT',
        'SF_TCR_F_CLON_FREQ',
        'SF_TCR_F_CONTIGS']
    if not filt:
        key_list = [k.replace("_TCR_F_", "_TCR_UF_") for k in key_list]

    n = 0
    missing = []
    for fkey in key_list:
        fname = spipe.filepath(fkey, samp)
        have = utils.check_infile(fname, verb=False, toxic=False)
        if have:
            n += 1
        else:
            missing.append(fname)

    return (n, missing)


def finish_tcr_stats(spipe):
    """Merge tcr stats with other stats and save analysis csv files

    Handles combine mode files if needed

    Return nothing
    """
    # Combine mode
    if spipe.is_combine():
        # Run default combine; May average or sum tcr metrics at this point
        comb.comb_global_files(spipe)
        for samp in spipe.get_samples():
            comb.combine_single_dge_output(spipe, samp)

        # Explicity update tcr in stats files
        # Filtered or not tcr data source depends on parent; github issue #303
        filt = True if spipe.get_use_case("tcr_parent") else False
        finish_tcr_comb_stats(spipe, filt=filt)

    # Non-combine
    else:
        # If have a parent, use filtered
        fkey = "SF_TCR_F_STAT" if spipe.get_use_case("tcr_parent") else "SF_TCR_UF_STAT"

        stat_dict = spipe.get_seq_pipe_stats_dict()
        for samp in spipe.get_samples():
            tcr_df = spipe.read_csv(fkey, samp=samp)
            val_col = tcr_df.columns[0]
            new_dict = tcr_df[val_col].to_dict()
            new_dict.update(stat_dict)
            new_dict["sample_well_count"] = samp.num_wells(all_sublibs=True)
            # Reads per sample
            stats = pb_stats.get_samp_reads_info(spipe, samp, stat_dict)
            new_dict.update(stats)

            # Saves both all and summary stats
            pb_stats.save_analysis_summary(spipe, samp, new_dict)


def finish_tcr_comb_stats(spipe, filt=False):
    """Finish combine mode stats for tcr

    Replace combined tcr values with explict values from combind tcr metrics file
    Writes updated stats files

    Returns nothing
    """
    fkey = "SF_TCR_F_STAT" if filt else "SF_TCR_UF_STAT"

    # Now update with explicitly calculated combined values
    for samp in spipe.get_samples():
        # Replace tcr stats in all stats
        stats_df = spipe.get_stats_data(samp)
        fkey = "SF_TCR_F_STAT" if spipe.get_use_case("tcr_parent") else "SF_TCR_UF_STAT"
        tcr_df = spipe.read_csv(fkey, samp=samp)
        val_col = tcr_df.columns[0]
        # Each new TCR one at a time
        for row in tcr_df.index:
            val = tcr_df.at[row, val_col]
            # print(row, val, sep='\t')
            stats_df.at[row, "combined"] = val

        # Save full and summary again
        pb_stats.save_analysis_summary(spipe, samp, stats_df, col="combined")


############################################################
def update_tcr_outputs(spipe):
    """Revise up first-pass file outputs to updated versions"""
    # Clean up metrics file to simple csv
    for samp in spipe.get_samples():
        update_metrics_file(spipe, samp, filt=False)
        update_metrics_file(spipe, samp, filt=True)


def update_metrics_file(spipe, samp, filt=False):
    """Clean up stat metric file"""
    fkey = "SF_TCR_F_STAT" if filt else "SF_TCR_UF_STAT"
    fname = spipe.filepath(fkey, samp)

    update_fname = fname.replace(".csv", ".txt")

    file_lines = utils.lines_from_fname(
        update_fname, verb=False, comment="#", all_lines=False
    )
    if not file_lines:
        return True

    index = []
    vals = []
    for line in file_lines:
        ind, v = line.split(",")
        if ind.startswith("tcr_"):
            index.append(ind)
            vals.append(v)

    tcr_df = pd.DataFrame(vals, columns=["value"], index=index)
    tcr_df.index.name = "statistic"

    spipe.write_df(tcr_df, fkey, samp=samp)
    return True
