#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Alignment related stuff
#
import re
import pandas as pd
from collections import Counter
import gzip

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
    import crispr
    import maptarg
    import utils
else:
    from splitpipe import bcutils
    from splitpipe import crispr
    from splitpipe import maptarg
    from splitpipe import utils


# ---------------------------------------------------------------------------


def have_preprocess_ins(spipe, verb=True):
    """Check if preprocess input files exist

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


def have_preprocess_outs(spipe, verb=False):
    """Check if preprocess output files exist

    Return True/False, list of files
    """
    ok = False
    fastq_new = spipe.filepath("PF_FASTQ_BC", None)
    check_files = [fastq_new]
    # If list is good, we have inputs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def run_preprocess(spipe):
    """Run preprocessing step to prepare fastq for alignment"""
    spipe.report_proc_step("Preprocess")
    ok = True
    # Input and output status and list
    i_ok, i_list = have_preprocess_ins(spipe)
    o_ok, o_list = have_preprocess_outs(spipe)
    # Run only if don't have output or fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for pre processing: {bad_list}"
            spipe.set_problem(story)
        else:
            preprocess_fastq(spipe)
    else:
        spipe.report_run_story("Skipping preprocess and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")

    ok = spipe.no_problems()
    spipe.report_proc_step("Preprocess", status=ok)
    return ok


# Feedback frequency
PRE_UPDATE = 200000


def preprocess_fastq(spipe):
    """Performs steps before running the alignment.

    Merge processed barcode data from R2 with sequence from R1

    Return status
    """
    fastq1 = spipe.get_par_val("fq1", as_path=True)
    fastq2 = spipe.get_par_val("fq2", as_path=True)
    fastq_new = spipe.filepath("PF_FASTQ_BC_RAW", None)

    tso_seq = spipe.get_par_val("tso_seq", as_str=True)
    tso_len = len(tso_seq)
    tso_read_limit = spipe.get_par_val("tso_read_limit", as_int=True)
    # Min read lengths; R1 is parameter; R2 based on amplicon seq
    r1_min_len = spipe.get_par_val("pre_min_R1_len", as_int=True)
    r2_min_len = len(spipe.get_par_val("bc_amp_seq", as_str=True))

    # Barcode stuff
    bc_info = spipe.get_bc_info()
    bc_max_edist = spipe.get_par_val("bc_edit_dist", as_int=True)
    # List of dicts, zero-pad to index match bc round
    bc_dicts = bc_info.get_bc_dict(as_list=True)
    # Amplicon offsets for poly-N and barcodes
    bc_starts, bc_ends = bc_info.get_amp_starts_ends()

    # Mapping of bc seqs and various things for name
    bc1_seq_to_windzp = bc_info.get_seq_to_windzp(1)
    bc2_seq_to_windzp = bc_info.get_seq_to_windzp(2)
    bc3_seq_to_windzp = bc_info.get_seq_to_windzp(3)
    bc1_seq_to_type = bc_info.get_seq_to_type(1)
    bc1_seq_to_bci = bc_info.get_seq_to_bci(1)
    bc2_seq_to_bci = bc_info.get_seq_to_bci(2)
    bc3_seq_to_bci = bc_info.get_seq_to_bci(3)

    # Sampled read bounds; Only for R1 processing here; R2 loaded as dataframe
    samp_first, samp_last, samp_story = spipe.get_par_slice("fastq_samp_slice")
    n_samp = samp_last - samp_first

    # Get R2 fastq sample, with barcodes split out
    ok, fq2_df = spipe.get_fastq_samp_df("r2", split_bc=True)
    if not ok:
        spipe.set_problem("Failed to get R2 reads sample")
        return False

    # Init counters
    read_stats = {
        "number_of_reads": 0,
        "reads_valid_barcode": 0,
        "reads_too_short": 0,
        "reads_ambig_bc": 0,
        "reads_tso_trim": 0,
        "reads_Q30_sampled": 0,
        "cDNA_Q30": 0,
        "bc_edit_dist_NA": 0,
    }
    for d in range(bc_max_edist + 1):
        key = f"bc_edit_dist_{d}"
        read_stats[key] = 0

    # Q30 for barcodes and polyN
    spipe.report_run_story2("Collecting R2 Q30 stats")
    q30_dict = bcutils.fq2_df_Q30_stats(fq2_df)
    read_stats.update(q30_dict)
    # Tally of read lengths and seqs
    valid_R1_lengths = Counter()
    # Three barcode rounds with 1-based indexing
    valid_bc_counts = [None]
    valid_bc_counts.append(Counter())
    valid_bc_counts.append(Counter())
    valid_bc_counts.append(Counter())
    # Counting barcode combos for two plates?
    if spipe.num_sample_plates() > 1:
        valid_2_counts = Counter()
    else:
        valid_2_counts = None

    # Return dict of perfect bc combos (for rescue) and threshold (to decide if ambig)
    spipe.report_run_story2("Collecting perfect bc counts for correction")
    counts, ct_thresh, bc_stats = bcutils.get_perfect_bc_counts(fq2_df, bc_info)
    read_stats.update(bc_stats)

    # Any problem (e.g. too few seqs) yeilds None
    if not counts:
        story = "Failed to process fq2 sample for for bc correction"
        spipe.set_problem(story)
        return False

    story = f"Got bc count_threshold {ct_thresh}, {len(counts)} counts"
    spipe.report_run_story2(story)

    # Keeping no-barcode reads?
    keep_no_bc = spipe.get_par_val("pre_keep_no_bc_fq", as_bool=True)
    if keep_no_bc:
        nobc1_fname = spipe.filepath("PF_FASTQ_NOBC_R1", None)
        nobc2_fname = spipe.filepath("PF_FASTQ_NOBC_R2", None)
        nobc1_out = open(nobc1_fname, "w")
        nobc2_out = open(nobc2_fname, "w")

    # Read name checking flags
    check_name_match = spipe.get_par_val("pre_check_name_match", as_bool=True)
    check_name_R1R2 = spipe.get_par_val("pre_check_name_R1R2", as_bool=True)

    # Header parsing regular expression
    keep_seq_index = spipe.get_par_val("pre_head_keep_seq_index", as_bool=True)
    if keep_seq_index:
        # regex has space so need raw
        head_regex = spipe.get_par_val("pre_head_regex", as_raw=True)
        HEAD_REGEX = re.compile(head_regex)
        timestamp = utils.ymdm_timestamp()

    # seq_index_set = set()

    # Open input and output fastq files
    spipe.report_run_story2("Correcting barcodes for fastq")
    spipe.report_run_story2(f"Collecting R1 Q30 stats; {n_samp} samples ({samp_story})")
    with gzip.open(fastq1, "rb") as f1_file, gzip.open(fastq2, "rb") as f2_file, open(
        fastq_new, "w"
    ) as fout:
        n = 0
        while True:
            # Each file, header, seq, strand, and quality lines
            head1 = f1_file.readline()
            seql1 = f1_file.readline()
            strl1 = f1_file.readline()
            qual1 = f1_file.readline()
            head2 = f2_file.readline()
            seql2 = f2_file.readline()
            strl2 = f2_file.readline()
            qual2 = f2_file.readline()
            # All done?
            if (not seql1) or (not seql2):
                break
            n += 1

            header1 = head1.decode().strip()
            header2 = head2.decode().strip()
            seq1 = seql1.decode().strip()
            seq2 = seql2.decode().strip()
            strand1 = strl1.decode().strip()
            strand2 = strl2.decode().strip()
            quality1 = qual1.decode().strip()
            quality2 = qual2.decode().strip()

            head_is_good = True
            name1 = name2 = read1 = read2 = index1 = ""
            # Parse and check headers
            if check_name_match:
                if keep_seq_index:
                    mat1 = HEAD_REGEX.match(header1)
                    mat2 = HEAD_REGEX.match(header2)
                    if mat1 and mat2:
                        name1 = mat1[1]
                        name2 = mat2[1]
                        read1 = mat1[2]
                        read2 = mat2[2]
                        index1 = mat1[3]
                        # seq_index_set.add(index1)
                    else:
                        head_is_good = False
                else:
                    name1 = header1.split()[0]
                    name2 = header2.split()[0]
                    read1 = header1.split()[1][0]
                    read2 = header2.split()[1][0]
                # Names and read numbers must match
                if head_is_good:
                    if check_name_R1R2:
                        if (name1 != name2) or (read1 != "1") or (read2 != "2"):
                            head_is_good = False
                    else:
                        if name1 != name2:
                            head_is_good = False

            if not head_is_good:
                story = f"Fastq R1 R2 header mismatch at record {n}\n"
                story += f"Fastq files: {fastq1} {fastq2}\n"
                story += f"Header fq1: {header1}\n"
                story += f"Header fq2: {header2}\n"
                story += f"Names |{name1}| |{name2}|\n"
                story += f"Reads |{read1}| |{read2}|\n"
                # Don't stop
                if True:
                    spipe.report_run_story2(story)
                    story = "Turning off fastq header checking; 'pre_check_name_match'"
                    spipe.report_warning(story)
                    check_name_match = False
                else:
                    spipe.set_problem(story)
                    print(story)
                    return False

            read_stats["number_of_reads"] += 1
            good_reads = no_ambigs = True
            # If either read is too short, pair is bad
            if (len(seq1) < r1_min_len) or (len(seq2) < r2_min_len):
                good_reads = False
                read_stats["reads_too_short"] += 1

            # Get barcode and polyN subseqs from R2 seq
            if good_reads:
                polyN = seq2[bc_starts[0] : bc_ends[0]]
                rbc1 = seq2[bc_starts[1] : bc_ends[1]]
                rbc2 = seq2[bc_starts[2] : bc_ends[2]]
                rbc3 = seq2[bc_starts[3] : bc_ends[3]]
                # No ambig bases allowed (Only check for N?)
                if ("N" in polyN) or ("N" in rbc1) or ("N" in rbc2) or ("N" in rbc3):
                    no_ambigs = False
                    read_stats["reads_ambig_bc"] += 1

            # Get corrected barcodes (if possible)
            if good_reads and no_ambigs:
                # Returns three bc seq strings, or empty strings if no match
                # Also returns maximum edit distance and sum of these
                bc1, bc2, bc3, ed_max, ed_sum = bcutils.correct_barcodes(
                    rbc1, rbc2, rbc3, counts, ct_thresh, bc_dicts, bc_max_edist
                )
                # Edit dist stats
                if ed_max < 0:
                    ed_max = "NA"
                ed_key = f"bc_edit_dist_{ed_max}"
                read_stats[ed_key] += 1

            # Process if good, no ambigs and have all barcodes
            if good_reads and no_ambigs and bc1 and bc2 and bc3:
                read_stats["reads_valid_barcode"] += 1

                # Jump past TSO if at start of read, else may not align
                #   (STAR may not allow 5' non-alignments ?)
                TSO_location = seq1.find(tso_seq)
                if 0 <= TSO_location < tso_read_limit:
                    seq1 = seq1[TSO_location + tso_len :]
                    quality1 = quality1[TSO_location + tso_len :]
                    read_stats["reads_tso_trim"] += 1

                # Save length after any trimming
                valid_R1_lengths[len(seq1)] += 1
                valid_bc_counts[1][bc1] += 1
                valid_bc_counts[2][bc2] += 1
                valid_bc_counts[3][bc3] += 1
                if valid_2_counts is not None:
                    valid_2_counts[(bc1, bc2)] += 1

                # Build string encoding
                # Well (cell) ids
                w1 = bc1_seq_to_windzp[bc1]
                w2 = bc2_seq_to_windzp[bc2]
                w3 = bc3_seq_to_windzp[bc3]
                well_id = f"{w1}_{w2}_{w3}"
                # Primer type; Round 1 only
                ptype = bc1_seq_to_type[bc1]
                # Barcode int indexes
                i1 = bc1_seq_to_bci[bc1]
                i2 = bc2_seq_to_bci[bc2]
                i3 = bc3_seq_to_bci[bc3]
                bc_bci = f"{i1}_{i2}_{i3}"
                # Corrected barcode seqs (not raw, from R2)
                # bc_str = f"{rbc1}_{rbc2}_{rbc3}"
                bc_str = f"{bc1}_{bc2}_{bc3}"
                # v0.9.4 header line as: w1_w2_w3__pt__i1_i2_i3__bc1_bc2_bc3__polyN
                # Combined well indexes, primer type, bc indexes, raw bc seqs, poly N
                name_parts = [well_id, ptype, bc_bci, bc_str, polyN]
                if keep_seq_index:
                    index_str = f"{timestamp}_{index1}"
                    name_parts.append(index_str)
                head = "@" + "__".join(name_parts)

                print(head, file=fout)
                print(seq1, file=fout)
                print(strand1, file=fout)
                print(quality1, file=fout)

            # Not good / don't have three good barcodes
            else:
                if keep_no_bc:
                    print(header1, file=nobc1_out)
                    print(seq1, file=nobc1_out)
                    print(strand1, file=nobc1_out)
                    print(quality1, file=nobc1_out)
                    print(header2, file=nobc2_out)
                    print(seq2, file=nobc2_out)
                    print(strand2, file=nobc2_out)
                    print(quality2, file=nobc2_out)

            # Collect stats regardless of read regardless bc-match, ambigs
            # Only "good_read" (i.e. long enough) and within sampling window
            if good_reads and (samp_first <= n <= samp_last):
                read_stats["reads_Q30_sampled"] += 1
                # Only read1 qual; Ignore last char (line return)
                # read_stats['cDNA_Q30'] += np.mean([ord(c)>62 for c in quality1[:-1]])
                read_stats["cDNA_Q30"] += bcutils.seq_qual_score(quality1[:-1])

            if (n % PRE_UPDATE) == 0:
                per = utils.report_percent_str(read_stats["reads_valid_barcode"], den=n)
                spipe.report_run_story3(f"Processed {n} fastq records; Valid bc {per}")

    # More than one seq index?
    # if len(seq_index_set) > 1:
    #    spipe.report_warning(f"Found multiple sequencing indexes: {len(seq_index_set)}")

    # Close no barcode fastq if saving these
    if keep_no_bc:
        nobc1_out.close()
        nobc2_out.close()
        spipe.report_run_story2(f"No-barcode R1 fastq: {nobc1_fname}")
        spipe.report_run_story2(f"No-barcode R2 fastq: {nobc2_fname}")

    # Feedback
    n_reads = read_stats["number_of_reads"]
    n_valbc = read_stats["reads_valid_barcode"]
    n_short = read_stats["reads_too_short"]
    n_amb = read_stats["reads_ambig_bc"]
    per_bc = utils.report_percent_str(n_valbc, n_reads)
    per_sh = utils.report_percent_str(n_short, n_reads)
    per_am = utils.report_percent_str(n_amb, n_reads)
    story = f"Preprocess total {n_reads} reads, {n_valbc} valid BC {per_bc}"
    spipe.report_run_story2(story)
    story = f"Preprocess too short {n_short} {per_sh}, with ambig {n_amb} {per_am}"
    spipe.report_run_story2(story)

    # Save stats to files
    pre_write_seq_pipe_stats(
        spipe, read_stats, valid_R1_lengths, valid_bc_counts, valid_2_counts
    )

    preprocess_clean_up(spipe)
    return n


def pre_write_seq_pipe_stats(
    spipe, read_stats, valid_R1_lengths, valid_counts, valid_2_counts=None, round_to=4
):
    """Write preprocess step sequencing and pipeline stats to file

    Returns nothing
    """
    seq_stats = {}
    pipe_stats = {}
    pipe_stats["number_of_reads"] = read_stats["number_of_reads"]
    seq_stats["reads_Q30_sampled"] = read_stats["reads_Q30_sampled"]

    # Denominator = Q30 sample count
    q30_den = read_stats["reads_Q30_sampled"]
    if q30_den < 1:
        q30_den = 1
    # Copy various things into different outputs
    for k, val in read_stats.items():
        if (k == "reads_Q30_sampled") or (k == "number_of_reads"):
            continue
        # Q30 values rounded; Normalize for cDNA; Barcodes done already
        if "_Q30" in k:
            if k == "cDNA_Q30":
                val = round(val / q30_den, round_to)
            else:
                val = round(val, round_to)
            seq_stats[k] = val
        elif k.startswith("bc_edit_") or k.startswith("reads_") or ("perfrac" in k):
            pipe_stats[k] = val

    # Cook up dataframes and save
    seq_df = pd.DataFrame.from_dict(seq_stats, orient="index", columns=["value"])
    seq_df.index.name = "statistic"
    spipe.write_df(seq_df, "PF_STAT_SEQ")

    pipe_df = pd.DataFrame.from_dict(pipe_stats, orient="index", columns=["value"])
    pipe_df.index.name = "statistic"
    # Dump all lines to log; n_rows=0
    spipe.write_df(pipe_df, "PF_STAT_PIPE", n_rows=0)

    # Valid R1 read lengths
    new_df = pd.DataFrame.from_dict(valid_R1_lengths, orient="index", columns=["count"])
    new_df.index.name = "length"
    spipe.write_df(new_df, "PF_BCC_READ_LENS")

    # Valid bc counts
    bcutils.write_val_bc_counts(spipe, valid_counts, valid_2_counts=valid_2_counts)


def preprocess_clean_up(spipe):
    """Handle clean up"""
    # Compress barcode corrected fastq
    fastq_new = spipe.filepath("PF_FASTQ_BC_RAW", None)
    spipe.report_run_story2(f"Compressing {fastq_new}")
    comp_lev = spipe.get_par_val("gzip_compresslevel", as_int=True)
    utils.check_and_gzip_file(fastq_new, comp_lev=comp_lev)
    # Don't need to hold sampled read info
    spipe.free_fastq_samp_df()


# --------------- Align mode ------------------------------------------------


def have_align_ins(spipe, verb=True):
    """Check if alignment inputs exist

    Return True/False, list of files
    """
    ok = False
    # Fastq and stat files, genome dir or guides
    infastq = spipe.filepath("PF_FASTQ_BC", None)
    stats = spipe.filepath("PF_STAT_PIPE", None)
    check_files = [infastq, stats]

    if not spipe.is_focal_crispr(map_direct=True):
        genome_dir = spipe.get_par_val("genome_dir", as_path=True)
        check_files.append(genome_dir)

    # If list is good, we have ins
    if utils.check_infile(check_files, verb=verb):
        ok = True

    # Need to check min number of reads
    if ok:
        al_min = spipe.get_par_val("align_min_start_reads", as_int=True)
        stat_dict = spipe.get_seq_pipe_stats_dict(verb=False)
        n_reads = int(stat_dict["number_of_reads"])
        if n_reads < al_min:
            ok = False
            story = f"Number of reads {n_reads} too few for alignment (min {al_min})"
            check_files = [story]

    return ok, check_files


def have_align_outs(spipe, verb=False):
    """Check if alignment outputs exist

    Return True/False, list of files
    """
    ok = False
    # bam file
    bam = spipe.filepath("PF_BAM_SORTED", None)
    check_files = [bam]
    # If list is good, we have outs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def run_align(spipe):
    """Run alignment step"""
    spipe.report_proc_step("Alignment")
    ok = True

    # Input and output status and list
    i_ok, i_list = have_align_ins(spipe)
    o_ok, o_list = have_align_outs(spipe)

    # Run only if don't already have outputs, or if fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for align: {bad_list}"
            spipe.set_problem(story)
        else:
            # Get read count and make sure enough
            stat_dict = spipe.get_seq_pipe_stats_dict()
            n_reads = stat_dict.get("reads_valid_barcode", 0)
            min_reads = spipe.get_par_val("align_min_start_reads", as_int=True)
            if n_reads < min_reads:
                story = f"Too few valid-barcode reads {n_reads} to align; Need at least {min_reads}"
                spipe.set_problem(story)
            else:
                if spipe.is_focal_crispr(map_direct=True):
                    ok = crispr.map_fastq_guides(spipe)
                else:
                    ok = maptarg.proc_map_targ_counts(spipe)
                    if ok:
                        ok = run_star(spipe)
                # Only clean if good; If issue, leave outputs
                if ok:
                    align_clean_up(spipe)
    else:
        spipe.report_run_story("Skipping alignment and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")

    ok = spipe.no_problems()
    spipe.report_proc_step("Alignment", status=ok)
    return ok


def get_STAR_align_call(spipe, clean_ss=True):
    """Get command to call STAR for alignment

    Return string
    """
    infastq = spipe.filepath("PF_FASTQ_BC", None)
    outpref = spipe.filepath("PF_FASTQ_BC_PREF", None)
    genome_dir = spipe.get_par_val("genome_dir", as_path=True)
    nthreads = spipe.get_par_val("nthreads", as_int=True)

    # Multi-mapper filter limit
    filt_mm_max = spipe.get_par_val("align_star_outFilterMultimapNmax", as_int=True)
    xpath = spipe.get_exec("STAR")

    # Cook up input command line then exec
    # --outSAMtype BAM Unsorted for bam output
    # --readFilesCommand zcat   For gzipped barcode fasta input
    star_command = f"{xpath} \
        --genomeDir {genome_dir} \
        --runThreadN {nthreads} \
        --readFilesIn {infastq} \
        --outFileNamePrefix {outpref} \
        --outSAMtype BAM Unsorted \
        --readFilesCommand zcat \
        --outFilterMultimapNmax {filt_mm_max} \
        "

    # Keeping unmapped reads too?
    if spipe.get_par_val("align_star_keep_no_map", as_bool=True):
        # This keeps unmapped reads 'within' the normal bam output
        star_command += " --outSAMunmapped Within"

    # Any extra args?
    exargs = spipe.get_par_val("align_star_extra_args", as_str=True)
    if exargs:
        star_command += " "
        star_command += exargs.replace(",", " ")

    # Clean to single space?
    if clean_ss:
        star_command = re.sub(" +", " ", star_command)

    return star_command


def run_star(spipe):
    """Align reads using STAR"""
    command = get_STAR_align_call(spipe)

    ok, callout, ex_story = utils.call_exec(command, verb=True)
    # Story to log only (call verb=True already prints outcome)
    spipe.report_run_story(ex_story, to_log=True)

    # Output saved to file
    ver = spipe.get_version()
    fpath = spipe.filepath("PF_STAR_ALIGN_OUT", None, genome_dir="output_dir")
    utils.file_header(ofile=fpath, story="STAR align call", ver=ver)
    with open(fpath, "a") as OUTFILE:
        com_str = utils.wrap_exe_com_args(command)
        print("\n# Call:", file=OUTFILE)
        print(com_str, file=OUTFILE)
        print("\n# Output:", file=OUTFILE)
        print(callout, file=OUTFILE)
        spipe.report_run_story(f"STAR output saved to {fpath}")

    if ok:
        parse_star_output(spipe)
    else:
        spipe.set_problem("Subprocess failed: " + command)
        spipe.set_problem(ex_story)

    return ok


def parse_star_output(spipe):
    """Parse output from STAR to collect stats and save to file"""
    n_input = n_unique = n_multi = n_toomany = 0

    fname = spipe.filepath("PF_FASTQ_FIN_OUT", None)
    # Explictly open to crash on fail
    INFILE = open(fname)
    for line in INFILE:
        if line.find("Number of input reads") >= 0:
            n_input = int(line.split("|")[1])
        elif line.find("Uniquely mapped reads number") >= 0:
            n_unique = int(line.split("|")[1])
        elif line.find("Number of reads mapped to multiple loci") >= 0:
            n_multi = int(line.split("|")[1])
        elif line.find("Number of reads mapped to too many loci") >= 0:
            n_toomany = int(line.split("|")[1])
    INFILE.close()

    perc = utils.report_percent_str(n_unique, den=n_input)
    story = f"Alignment result: {n_unique} reads of {n_input} map uniquely; {perc}"
    spipe.report_run_story(story)

    # Update stats file; Load data frame, update rows, write back out
    pipe_df = spipe.read_csv("PF_STAT_PIPE", names="statistic,value", verb=False)
    pipe_df.loc["reads_align_input", "value"] = n_input
    pipe_df.loc["reads_align_unique", "value"] = n_unique
    pipe_df.loc["reads_align_multimap", "value"] = n_multi
    pipe_df.loc["reads_too_many_loci", "value"] = n_toomany
    pipe_df.index.name = "statistic"
    # Dump all lines to log; n_rows=0
    spipe.write_df(pipe_df, "PF_STAT_PIPE", n_rows=0)


def align_clean_up(spipe):
    """Handle clean up"""
    # Compress barcode corrected fastq
    fname = spipe.filepath("PF_FASTQ_SJ_TAB", None)
    spipe.report_run_story2(f"Compressing {fname}")
    comp_lev = spipe.get_par_val("gzip_compresslevel", as_int=True)
    utils.check_and_gzip_file(fname, comp_lev=comp_lev)


def have_postprocess_ins(spipe, verb=True):
    """Check if postprocess input files exist

    Return True/False, list of files
    """
    ok = False
    # (unsorted) bam files
    bam = spipe.filepath("PF_BAM_RAW", None)
    check_files = [bam]
    # If list is good, we have ins
    if utils.check_infile(check_files, verb=verb):
        ok = True

    # Need to check min mapping fraction
    if ok:
        warn_mfrac = spipe.get_par_val("post_warn_map_frac", as_float=True)
        min_mfrac = spipe.get_par_val("post_min_map_frac", as_float=True)
        stat_dict = spipe.get_seq_pipe_stats_dict(verb=False)
        mfrac = stat_dict["reads_align_unique"] / stat_dict["reads_align_input"]
        # Under min = fail
        if mfrac < min_mfrac:
            ok = False
            story = f"Mapping fraction {mfrac} too small (min {min_mfrac})"
            check_files = [story]
        # Warn?
        elif mfrac < warn_mfrac:
            spipe.report_warning(f"Low mapping fraction: {mfrac}")

    return ok, check_files


def have_postprocess_outs(spipe, verb=False):
    """Check if postprocess output files exist

    Return True/False, list of files
    """
    ok = False
    # bam file
    bam = spipe.filepath("PF_BAM_SORTED", None)
    check_files = [bam]
    # If list is good, we have outs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def run_postprocess(spipe):
    """Run post-alignment processing step

    Return status
    """
    spipe.report_proc_step("Postprocess")
    ok = True
    # Input and output status and list
    i_ok, i_list = have_postprocess_ins(spipe)
    o_ok, o_list = have_postprocess_outs(spipe)
    # Run only if don't have output or fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for post processing: {bad_list}"
            spipe.set_problem(story)
            ok = False
        else:
            ok = sort_sam(spipe)
    else:
        spipe.report_run_story("Skipping postprocess and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")

    ok = spipe.no_problems()
    spipe.report_proc_step("Postprocess", status=ok)
    return ok


def sort_sam(spipe):
    """Sort samfile by header (i.e. cell_barcodes)

    Return status
    """
    xpath = spipe.get_exec("samtools")
    inbam = spipe.filepath("PF_BAM_RAW", None)
    outbam = spipe.filepath("PF_BAM_SORTED", None)
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    # Temp prefix from output bam file base name
    # From samtools:
    #   -T PREFIX  Write temporary files to PREFIX.nnnn.bam
    tprefix = spipe.filepath("PF_BAM_SORTED_PREF", None)

    # Sort on name (-n)
    command = f"{xpath} sort -n --threads {nthreads} -T {tprefix} -o {outbam} {inbam}"

    # Any extra args?
    exargs = spipe.get_par_val("post_samtools_extra_args", as_str=True)
    if exargs:
        command += " "
        command += exargs.replace(",", " ")

    ok, callout, ex_story = utils.call_exec(command, verb=True)
    # Story to log only
    spipe.report_run_story(ex_story, to_log=True)
    if not ok:
        spipe.set_problem("Subprocess failed: " + command)
        spipe.set_problem(ex_story)
        return False

    return True
