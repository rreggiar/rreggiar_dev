#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


# ---------------------------------------------------------------------------
def have_inqc_ins(spipe, verb=True):
    """Check if inQC input files exist

    Return True/False, list of files
    """
    ok = False
    check_files = []
    # Input fastq files
    fq1 = spipe.get_par_val("fq1", as_path=True)
    fq2 = spipe.get_par_val("fq2", as_path=True)
    if fq1 and fq2:
        check_files.append(fq1)
        check_files.append(fq2)
        # If list is good, we have inputs
        if utils.check_infile(check_files, verb=verb):
            ok = True
    else:
        story = f"Bad fastq file(s): fq1={fq1}, fq2={fq2}"
        spipe.set_problem(story)
    return ok, check_files


def have_fastqc_outs(spipe, prefix=""):
    """Check if fastQC output files exist

    Return True/False, list of files
    """
    ok = False
    file_list = []

    # Look for files in output dir (which may not exist itself)
    fqc_path = spipe.filepath("DIR_FASTQC", None)
    for fname in utils.fnames_for_path(path=fqc_path):
        if fname.endswith(f"{prefix}fastqc.html") or fname.endswith(
            f"{prefix}fastqc.zip"
        ):
            file_list.append(fname)

    if file_list:
        ok = True
    return ok, file_list


def run_inqc(spipe):
    """Handle input QC"""
    spipe.report_proc_step("InQC")
    ok = True
    # Input status and list
    i_ok, i_list = have_inqc_ins(spipe)
    if not i_ok:
        bad_list = utils.bad_infile_list(i_list)
        story = f"Don't have inputs for inQC: {bad_list}"
        spipe.set_problem(story)
        ok = False

    # Simple file stats
    if ok:
        report_in_fastq_stats(spipe)

    # Min length checks
    if ok:
        spipe.report_run_story("Checking fastq read lengths")
        if not fastq_lengths_ok(spipe):
            spipe.set_problem("Input fastq file(s) failed read length requirements")
            ok = False

    # BC match checks
    if ok:
        spipe.report_run_story("Checking perfect barcode fractions")
        if not valid_bc_frac_ok(spipe):
            spipe.set_problem("Input fastq file(s) failed valid barcode fraction")
            ok = False

    # FastQC
    if ok:
        if spipe.get_par_val("inqc_fastqc_run", as_bool=True):
            o_ok, o_list = have_fastqc_outs(spipe)
            if o_ok:
                spipe.report_run_story("Skipping fastQC input checks; Outputs exist")
                spipe.report_run_story2(f"Outputs: {o_list}")
            else:
                ok = run_fastqc(spipe)
        else:
            spipe.report_run_story("Skipping fastQC input checks; Parameter not set")

    ok = spipe.no_problems()
    spipe.report_proc_step("InQC", status=ok)
    return ok


def fastq_lengths_ok(spipe):
    """Check fastq files for minimal length requirements

    Return status
    """
    # R1 and R2 fastq samples
    for read in [1, 2]:
        ok, story = check_seq_good_len(spipe, read)
        if not ok:
            spipe.set_problem(story)
            return False
        spipe.report_run_story2(story)
    return True


def check_seq_good_len(spipe, read):
    """Check sample of reads (in dataframe) for length constraints

    Return tuple (status, story)
    """
    if read == 1:
        min_len = spipe.get_par_val("pre_min_R1_len", as_int=True)
        min_frac = spipe.get_par_val("inqc_min_good_R1len_frac", as_float=True)
    elif read == 2:
        min_len = len(spipe.get_par_val("bc_amp_seq", as_str=True))
        min_frac = spipe.get_par_val("inqc_min_good_R2len_frac", as_float=True)
    else:
        raise ValueError(f"check_seq_good_len given bad read: {read}")

    # Get seq sample
    ok, seq_df = spipe.get_fastq_samp_df(read)
    if not ok:
        story = f"Failed to get R{read} reads sample"
        return False, story

    min_reads = spipe.get_par_val("inqc_min_reads", as_int=True)
    if len(seq_df) < min_reads:
        story = f"Too few R{read} reads: Got {len(seq_df)} < {min_reads} min"
        return False, story

    # Check length fraction
    good = True
    good_frac = frac_df_seq_good_len(seq_df, min_len=min_len)
    # Min fractions of length-passing reads
    if good_frac < min_frac:
        story = f"Too few R{read} reads with OK length: {good_frac} < {min_frac}; min len {min_len}"
        good = False
    else:
        story = f"R{read} fastq read length good fraction: {good_frac}"

    return good, story


def frac_df_seq_good_len(seq_df, col="seq", min_len=0, max_len=0):
    """Return fraction of sequences in dataframe with length in range

    Return fraction
    """
    n_too_long = n_too_short = 0
    n_total = len(seq_df)
    # Series of seq lengths
    s_lens = seq_df[col].str.len()
    if min_len:
        n_too_short = sum(s_lens < min_len)
    if max_len:
        n_too_long = sum(s_lens > max_len)
    n_good = n_total - n_too_short - n_too_long
    return round(n_good / n_total, 4)


def valid_bc_frac_ok(spipe):
    """Check if valid barcode fractions are ok

    Return status
    """
    min_frac = spipe.get_par_val("inqc_min_perf_frac", as_float=True)
    warn_frac = spipe.get_par_val("inqc_warn_perf_frac", as_float=True)
    if (min_frac <= 0) and (warn_frac <= 0):
        return True

    # Read sample and barcodes
    _, fq2_df = spipe.get_fastq_samp_df("r2", split_bc=True)
    n_samp = len(fq2_df)
    bc_info = spipe.get_bc_info()

    # Collect perfect fractions
    fracs = []
    for r in range(3):
        bc_col = f"bc{r+1}"
        bc_set = bc_info.get_bc_seqs(r + 1, as_set=True)
        n_val = len(fq2_df[fq2_df[bc_col].isin(bc_set)])
        fracs.append(round(n_val / n_samp, 4))

    min_frac = min(fracs)
    spipe.report_run_story(f"Perfrect barcode fractions: {fracs}")
    if min_frac < min_frac:
        spipe.set_problem(f"Perfect bc fraction {min_frac} too low")
        return False
    if min_frac < warn_frac:
        spipe.report_warning(f"Perfect bc fraction {min_frac} is low")
    return True


def report_in_fastq_stats(spipe):
    """Report file stats for input fastqs

    Return nothing
    """
    fq1 = spipe.get_par_val("fq1", as_path=True)
    fq2 = spipe.get_par_val("fq2", as_path=True)
    fq1_str = utils.report_file_info_str(fq1)
    fq2_str = utils.report_file_info_str(fq2)
    spipe.report_run_story(f"Input fastq R1: {fq1_str}")
    spipe.report_run_story(f"Input fastq R2: {fq2_str}")


# Template for contaminant file
FASTQC_CON_TEMP = "{pkg_path}/config/fastqc-contaminant_list.txt"


def run_fastqc(spipe):
    """Run fastQC

    Return status
    """
    # FastQC path; Optional so may not have
    xpath = spipe.get_exec("fastqc")
    if not xpath:
        spipe.report_run_story("fastQC not found; Skipping")
        return True

    # Limit number of threads
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    nthreads = min(4, nthreads)

    # Output dir; Needs to exists for fastqc
    fqc_outdir = spipe.filepath("DIR_FASTQC", None)
    if not utils.make_subdir(fqc_outdir, verb=False, toxic=False):
        spipe.set_problem(f"Failed to create fastQC output dir {fqc_outdir}")
        return False
    # Output summary file
    out_sum = spipe.filepath("PF_FASTQC_SUM_OUT", None)

    # Inputs
    fq1 = spipe.get_par_val("fq1", as_path=True)
    fq2 = spipe.get_par_val("fq2", as_path=True)
    # Sample record range; Four lines per rec
    first_rec, last_rec, s_story = spipe.get_par_slice("fastq_samp_slice")
    first_line = first_rec * 4
    last_line = last_rec * 4
    num_lines = last_line - first_line

    story = f"Calling fastQC for inputs; Sample {first_rec} to {last_rec} fastq records"
    spipe.report_run_story(story)

    # contaminants file if found
    path = spipe.get_pkg_path()
    con_path = FASTQC_CON_TEMP.format(pkg_path=path)
    con_arg = ""
    if utils.check_infile(con_path, verb=False):
        # Leading space to concat to full arg list
        con_arg = f" --contaminants {con_path}"

    # Presume fastq may not be gzipped
    cat_call = "zcat" if fq1.endswith(".gz") else "cat"
    # Construct head/tail part of read sampling from first/last records
    samp_call = f"head -n {last_line} | tail -n {num_lines}"
    fqc_call = (
        f"{xpath} stdin --noextract --outdir {fqc_outdir} --threads {nthreads}"
        + con_arg
    )

    # Both R1 and R2 processed
    sum_out_parts = []
    for r, fq in [("R1", fq1), ("R2", fq2)]:
        # Command like:
        # zcat fastq.gz | head -n | tail -n | fastqc stdin
        command = f"{cat_call} {fq} | {samp_call} | {fqc_call}"

        ok, callout, ex_story = utils.call_exec(command, verb=True)
        if not ok:
            spipe.set_problem("Subprocess failed: " + command)
            spipe.set_problem(ex_story)
            return False

        # Get PASS/FAIL story from zip; This command yields summary.txt contents
        command = f"unzip -p {fqc_outdir}stdin_fastqc.zip stdin_fastqc/summary.txt"
        # vok=0 puts no limit on output
        ok, callout, ex_story = utils.call_exec(command, verb=False)
        if not ok:
            spipe.set_problem("Subprocess failed: " + command)
            spipe.set_problem(ex_story)
            return False
        # Collect cleaned up output
        for line in callout.split("\n"):
            if "stdin" in line:
                # Switch around, drop 'stdin' and put R1 or R2 first:
                #   From:   PASS	Basic Statistics	stdin
                #     To:   R1	PASS	Basic Statistics
                parts = [r] + line.split("\t")[:-1]
                sum_out_parts.append(parts)

        # Rename 'stdin' output files R1 or R2
        o_ok, o_list = have_fastqc_outs(spipe, prefix="stdin_")
        for oname in o_list:
            nname = oname.replace("stdin_", f"{r}_")
            o_str = utils.copy_file(oname, nname, move=True, verb=False)
            spipe.report_run_story2(o_str)

    ok = True

    # Save summarized results
    # Possibly report non-pass directly
    show_warn = spipe.get_par_val("inqc_fastqc_show_warn", as_bool=True)
    if show_warn:
        spipe.report_run_story2("Non-PASS fastQC results:")
    with open(out_sum, "w") as OUTFILE:
        ver = spipe.get_version()
        header = utils.file_header(story="Input fastQC test results", ver=ver)
        print(header, file=OUTFILE)
        for parts in sum_out_parts:
            line = "\t".join(parts)
            print(line, file=OUTFILE)
            if show_warn and (parts[1] != "PASS"):
                spipe.report_run_story2(f"fastQC {line}")
        spipe.report_run_story2(f"New file: {out_sum}")

    return ok
