#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import gzip
import numpy as np
import pandas as pd
import pysam

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


# ---------------------------------------------------------------------------
def get_genome_name(spipe):
    """Get focal/crispr specific genome_name"""
    return "gRNA"


def load_guide_info(spipe, add_gene_info=True):
    """Get guide info structure (like gene info) from file

    Return dict
    """
    new_dict = {}

    guide_fname = spipe.get_par_val("crsp_guides", as_path=True)
    spipe.report_run_story(f"Loading and preparing guides from {guide_fname}")
    guide_df = get_guide_df(guide_fname)
    if not isinstance(guide_df, pd.DataFrame):
        spipe.set_problem(f"Failed to load guides: {guide_fname}")
        return None

    if not check_guide_data_df(spipe, guide_df):
        spipe.set_problem(f"Problem with guides: {guide_fname}")
        return None

    new_dict["crsp_use_star"] = spipe.get_par_val("crsp_use_star", as_bool=True)
    new_dict["guide_source"] = guide_fname
    new_dict["guide_number"] = len(guide_df)

    genome = get_genome_name(spipe)
    new_dict["genome"] = genome
    # Add "gene" id and genome-prefix id to guide dataframe and save
    guide_df["gene_id"] = [f"gene_id_{x+1}" for x in list(guide_df.index)]
    guide_df["gpx_guide_name"] = f"{genome}_" + guide_df["guide_name"]
    new_dict["guide_df"] = guide_df

    # Mapping dicts
    new_dict["guide_seq_to_index"] = dict(
        zip(guide_df["guide_sequence"], guide_df.index)
    )
    new_dict["guide_name_to_gene"] = dict(
        zip(guide_df["guide_name"], guide_df["target_gene"])
    )
    new_dict["gpx_guide_name_to_gene_id"] = dict(
        zip(guide_df["gpx_guide_name"], guide_df["gene_id"])
    )
    new_dict["gpx_guide_name_to_guide_name"] = dict(
        zip(guide_df["gpx_guide_name"], guide_df["guide_name"])
    )

    # Mapping offset coords
    mc_dict = get_guide_map_coords(spipe, guide_df)
    if not mc_dict:
        spipe.set_problem(f"Problem getting guide coords: {guide_fname}")
        return None
    new_dict.update(mc_dict)

    # List of mismatch dict, one per guide length
    max_mm = spipe.get_par_val("crsp_match_hamming_dist", as_int=True)
    new_dict["guide_max_mm"] = max_mm
    mm_dict_list = []
    for glen in new_dict["guide_len_list"]:
        g_len_df = guide_df[guide_df["guide_len"] == glen]
        mm_dict_list.append(get_guide_match_dict(g_len_df, max_mm=max_mm))
    new_dict["guide_mm_list"] = mm_dict_list

    # Add 'genome' info fields
    if add_gene_info:
        new_dict = guide_add_gene_info(spipe, new_dict)

    return new_dict


def get_guide_map_coords(spipe, guide_df):
    """Get length and coordinate info for guide mapping to reads

    Return dict
    """
    new_dict = {}

    # Guide and full template length list, large to small
    new_dict["guide_len_list"] = sorted(guide_df["guide_len"].unique(), reverse=True)
    new_dict["full_len_list"] = sorted(
        guide_df["full_seq"].str.len().unique(), reverse=True
    )

    # If offset given, use that for start, else find
    g_off = spipe.get_par_val("crsp_guide_offset", as_int=True)
    if g_off < 1:
        if spipe.is_focal_crispr(map_direct=True):
            g_off = find_crsp_guide_R1_offset(spipe, guide_df)
        else:
            g_off = 1
    # 1-based coord; Less = problem
    if g_off < 1:
        return None

    # Read coords; Start = single value, Ends = large to small guide list
    new_dict["guide_start"] = g_off - 1
    new_dict["guide_end_list"] = [
        x + new_dict["guide_start"] for x in new_dict["guide_len_list"]
    ]
    return new_dict


def find_crsp_guide_R1_offset(spipe, guide_df):
    """Attempt to find crisrp guide offset into R1

    Return 1-based int
    """
    # Unique prefix seqs
    p_list = list(guide_df["prefix"].unique())
    if len(p_list) > 1:
        spipe.set_problem(
            f"Guide prefix cannot differ if no offset; Found {len(p_list)}"
        )
        return 0

    # Get 3' part of prefix
    pre_seq = p_list[0]
    ssize = spipe.get_par_val("crsp_pref_samp_len", as_int=True)
    if ssize > len(pre_seq):
        spipe.set_problem(
            f"Guide prefix len {len(pre_seq)} < crsp_pref_samp_len {ssize}"
        )
        return 0

    pre_seq = pre_seq[-ssize:]
    spipe.report_run_story(f"Looking for R1 guide offset with prefix 3' end {pre_seq}")
    ok, seq_df = spipe.get_fastq_samp_df("R1")
    if not ok:
        spipe.set_problem("Failed to get R1 read sample")
        return 0

    # List of find hit frequency then get top count
    find_vals = seq_df["seq"].apply(lambda s: s.find(pre_seq)).value_counts()
    n_top = find_vals.max()
    # Pass pref seq and shift to 1-based coord
    offset = find_vals.idxmax() + len(pre_seq) + 1

    f_score = round(n_top / len(seq_df), 3)
    spipe.report_run_story(f"Found guide offset: {offset}, score {f_score}")

    return offset


def guide_add_gene_info(spipe, guide_info):
    """Create gene_info fields from guide_info

    Return guide_info
    """
    guide_df = guide_info["guide_df"]

    guide_df["gene_id"] = [f"gene_id_{x+1}" for x in list(guide_df.index)]

    # Make all genes
    all_genes = pd.DataFrame()
    all_genes["gene_id"] = guide_df["gene_id"]
    all_genes["gene_name"] = guide_df["guide_name"]
    all_genes["genome"] = get_genome_name(spipe)

    guide_info["all_genes"] = all_genes

    guide_info["genome_list"] = [guide_info["genome"]]

    guide_info["gene_id_to_idx"] = {v: k for k, v in dict(all_genes["gene_id"]).items()}
    guide_info["gene_id_to_name"] = dict(
        zip(guide_df["gene_id"], guide_df["guide_name"])
    )

    genome = get_genome_name(spipe)
    guide_info["gene_id_to_genome"] = {g: genome for g in guide_df["gene_id"]}

    return guide_info


def get_guide_info_specs(guide_info):
    """Get settings for focal / crsipr (e.g. for json)

    Return dict
    """
    ndict = {}
    if guide_info:
        ndict["guide_source"] = guide_info["guide_source"]
        ndict["guide_number"] = guide_info["guide_number"]
        ndict["crsp_use_star"] = guide_info["crsp_use_star"]

        # STAR doesn't use these
        if not guide_info["crsp_use_star"]:
            ndict["guide_max_mm"] = int(guide_info["guide_max_mm"])
            # Cast for json (supports simple data types)
            ndict["guide_len_list"] = [int(x) for x in guide_info["guide_len_list"]]
            ndict["guide_start"] = int(guide_info["guide_start"])
            ndict["guide_end_list"] = [int(x) for x in guide_info["guide_end_list"]]
    return ndict


def get_guide_df(fname):
    """Read in and set guides from file

    Return dataframe
    """
    # Format may/may not have header... So load without and check
    guide_df = utils.read_csv(fname, header=None)
    if not isinstance(guide_df, pd.DataFrame):
        return None

    # Check if first line is header
    row_words = list(guide_df.iloc[0])
    found, _ = utils.list_find_match(row_words, "guide_name", nocase=True)
    if found:
        guide_df.columns = row_words
        guide_df = guide_df[1:]

    # Standardize at least first columns
    all_cols = list(guide_df.columns)
    first_cols = [
        "guide_name",
        "prefix",
        "guide_sequence",
        "suffix",
        "target_gene",
    ]
    guide_df.columns = first_cols + all_cols[len(first_cols) :]

    # Add guide len and full length seq from parts (e.g. for STAR)
    guide_df["guide_len"] = guide_df["guide_sequence"].str.len()
    guide_df["full_seq"] = (
        guide_df["prefix"] + guide_df["guide_sequence"] + guide_df["suffix"]
    )
    return guide_df


def check_guide_data_df(spipe, guide_df, check_pairwise=True):
    """Check guide data(frame) is all good

    Return status
    """
    # Make sure names are unique
    n_dup, story = utils.check_df_col_dup_counts(guide_df, "guide_name")
    if n_dup:
        spipe.set_problem(f"Guide names not unique; {n_dup} duplicate name(s)")
        # Preface newline to get table of dups on new lines
        spipe.set_problem("\n" + story)
        return False

    # Make sure only sequences and guides are unique
    only_char = "ACGT"
    n_bad, story = utils.check_df_col_str_content(
        guide_df, "guide_sequence", only_char=only_char
    )
    if n_bad:
        spipe.set_problem(f"Guide sequence(s) include bad chars (not {only_char})")
        # Pref to get table of dups on new lines
        spipe.set_problem("\n" + story)
        return False
    n_dup, story = utils.check_df_col_dup_counts(guide_df, "guide_sequence")
    if n_dup:
        spipe.set_problem(f"Guides not unique; {n_dup} duplicate seq(s)")
        # Pref to get table of dups on new lines
        spipe.set_problem("\n" + story)
        return False

    # Prefixes all have to be same len and can't differ if no offset
    # Only matters for direct mapping (not STAR)
    if spipe.is_focal_crispr(map_direct=True):
        p_list = list(guide_df["prefix"].unique())
        if spipe.get_par_val("crsp_guide_offset", as_int=True) < 1:
            if len(p_list) > 1:
                story = f"Guide prefix cannot differ if no offset; Found {len(p_list)}"
                spipe.set_problem(story)
                return False
        else:
            min_len = min([len(s) for s in p_list])
            max_len = max([len(s) for s in p_list])
            if max_len != min_len:
                story = (
                    f"Guide prefix length cannot differ; Found {min_len} to {max_len}"
                )
                spipe.set_problem(story)
                return False

    # Check of pairwise similarity conflicts; Any issues reported by check
    if check_pairwise:
        if not check_guide_pairwise(spipe, guide_df):
            return False

    spipe.report_run_story(f"Loaded {len(guide_df)} guide sequences; Good to use")
    story = guide_len_coord_story(guide_df)
    spipe.report_run_story(story)
    return True


def guide_len_coord_story(guide_df):
    """Get story describing guide lengths and coords

    Return str
    """
    p_str = df_col_str_len_story(guide_df, "prefix")
    g_str = df_col_str_len_story(guide_df, "guide_sequence")
    f_str = df_col_str_len_story(guide_df, "full_seq")
    story = f"Prefix length {p_str}, guide lengths {g_str}, full template seqs {f_str}"
    return story


def df_col_str_len_story(df, col, sep="-"):
    """Get str with length(s) for strings in given dataframe column

    return str
    """
    len1 = df[col].str.len().min()
    len2 = df[col].str.len().max()
    if len1 == len2:
        story = f"{len1}"
    else:
        story = f"{len1}{sep}{len2}"
    return story


def check_guide_pairwise(spipe, guide_df):
    """Check that guides don't have pairwise conflicts

    Return status
    """
    ok = True
    mm_dist = spipe.get_par_val("crsp_match_hamming_dist", as_int=True)

    # Self pairwise difference; Flag if <= target mismatch level
    spipe.report_run_story(
        f"Checking {len(guide_df)} guide pair-wise mismatch distsances; Min {mm_dist}"
    )
    rep_list = check_guide_mm_dists(guide_df, mm_dist)
    if rep_list:
        story = f"Guides differ <= target mismatch distance {mm_dist}\n"
        story += "\n".join(rep_list)
        spipe.set_problem(story)
        ok = False

    return ok


def hamming_dist(f_seq, s_seq, max_mm=0):
    """Get Hamming distance between two (same length) sequences

    If different lengths, only compares up to end of shorter seq

    Return int
    """
    diffs = 0
    for f_base, s_base in zip(f_seq, s_seq):
        if f_base != s_base:
            diffs += 1
            if max_mm and (diffs >= max_mm):
                break
    return diffs


def check_guide_mm_dists(guide_df, flag_mm):
    """Check guide-guide mismatch distances for too-similar conflicts

    guide_df = guides
    flag_mm = min mismatch Hamming dist to flag

    Collects list of strings reporting violations

    Return list
    """
    rep_list = []

    s_list = list(guide_df["guide_sequence"])
    n_list = list(guide_df["guide_name"])

    for i in range(len(s_list) - 1):
        f_seq = s_list[i]
        for j in range(i + 1, len(s_list)):
            s_seq = s_list[j]
            # Limit Hamming dist check max to min allowed plus 1
            hamdis = hamming_dist(f_seq, s_seq, max_mm=(flag_mm + 1))
            # Report if within min mismatch distance
            if hamdis <= flag_mm:
                f_name = n_list[i]
                s_name = n_list[j]
                rep_list.append(
                    f"Guide Hamming distance [{i}] {f_name} and [{j}] {s_name} = {hamdis}"
                )
    return rep_list


def gen_mm_seqs(seq, max_mm=1, pref=""):
    """Generate mismatch sequences based on given seq

    max_mm = max mismatches to yield
    pref = prefix sequence; Used for recursion

    Return list of sequences
    """
    if len(seq) < max_mm:
        return []

    # Self sequence (0 mismatch); Initially prefix is empty
    seq_list = [pref + seq]
    # 3 x N mismatches
    if max_mm > 0:
        for i in range(len(seq)):
            up_seq = seq[:i]
            dn_seq = seq[i + 1 :]
            for b in "ACGT":
                if b == seq[i]:
                    continue
                mm_seq = pref + up_seq + b + dn_seq
                seq_list.append(mm_seq)
                # More mismatches? Recur
                if (max_mm > 1) and dn_seq:
                    seq_list += gen_mm_seqs(
                        dn_seq, max_mm=max_mm - 1, pref=pref + up_seq + b
                    )
    return seq_list


def get_mm_cigar(f_seq, s_seq):
    """Get mismatch annotated CIGAR string for two seqs

    Assume seqs same lengths and only check ident; No InDels

    Return str
    """
    cigar = ""
    mat_len = mm_len = 0

    for f_base, s_base in zip(f_seq, s_seq):
        if f_base == s_base:
            mat_len += 1
            if mm_len:
                cigar += f"{mm_len}X"
                mm_len = 0
        else:
            mm_len += 1
            if mat_len:
                cigar += f"{mat_len}M"
                mat_len = 0
    # Final part
    if mat_len:
        cigar += f"{mat_len}M"
    if mm_len:
        cigar += f"{mm_len}X"

    return cigar


def get_guide_match_dict(guide_df, max_mm=1):
    """Get crispr guide read-direct-matching dict

    guide_df = dataframe with guide data
    max_mm = maximum mismatches to consider

    Dict keys = perfect and mismatch guide sequences;
         vals = [perfect seq, Hamming distance, cigar string]

    Return dict
    """
    mat_dict = {}
    for seq in guide_df["guide_sequence"]:
        seq_list = gen_mm_seqs(seq, max_mm=max_mm)
        for mm_seq in seq_list:
            cigar = get_mm_cigar(seq, mm_seq)
            hamdis = hamming_dist(seq, mm_seq)
            mat_dict[mm_seq] = [seq, hamdis, cigar]

    return mat_dict


# ------------------------- guides as genome / STAR stuff ------------------


def prep_guide_genome(spipe):
    """Prepare 'genome' guide genome for use (i.e. mapping or STAR)

    Return status
    """
    if not spipe.is_focal_crispr(crispr=True):
        return True

    ok = True
    # Direct mapping
    if spipe.is_focal_crispr(map_direct=True):
        guide_info = spipe.get_guide_info()
        spipe.write_df(guide_info["all_genes"], "PF_ALL_GENES")
    # STAR
    elif spipe.get_par_val("crsp_use_star", as_bool=True):
        ok = run_mkref_star_script(spipe)
    else:
        spipe.set_problem("Prep guide genome not direct, not STAR")
        ok = False
    return ok


def get_mkref_script_par(spipe):
    """Get script and parameter file paths to run mkref for focal / crispr

    Generates parameter file to be used

    Return status, script and par paths
    """
    ok = True
    # Top level package path
    path = spipe.get_pkg_path()
    call_path = f"{path}/scripts/mkref_fb_genome.sh"

    # Write parameter file
    par_path = spipe.filepath("PF_CRSP_MKREF_PAR", None)
    try:
        with open(par_path, "w") as OUTFILE:
            # Hardcoded for STAR
            mkref_star_genomeSAindexNbases = 3
            # crsp_use_star = spipe.get_par_val("crsp_use_star", as_bool=True)
            crsp_match_hamming_dist = spipe.get_par_val(
                "crsp_match_hamming_dist", as_int=True
            )

            print(
                f"mkref_star_genomeSAindexNbases\t{mkref_star_genomeSAindexNbases}",
                file=OUTFILE,
            )
            print(f"crsp_match_hamming_dist\t{crsp_match_hamming_dist}", file=OUTFILE)
            # print(f"crsp_use_star\t{crsp_use_star}", file=OUTFILE)
            # Don't need memory tracing
            print("log_memory_use False", file=OUTFILE)
        spipe.report_run_story2(f"Wrote focal/crispr mkref par file: {par_path}")

    except Exception as e:
        story = f"Failed to write focal/crispr parameter file; {e}"
        spipe.set_problem(story)
        ok = False

    return ok, call_path, par_path


def run_mkref_star_script(spipe):
    """Call 'split-pipe mkref' via script to create guide genome database

    Call mkref for STAR or direct matching, dictated by generated parameter file
    If new db is created, set associated parameters in spipe

    Return status
    """
    spipe.report_run_story("Will create focal/crispr genome with STAR")
    # Script info and cook up command line
    ok, script_path, par_path = get_mkref_script_par(spipe)

    if ok:
        # Script usage: <genome_name> <guides> <output_dir> <pars>
        gname = get_genome_name(spipe)
        guides = spipe.get_par_val("crsp_guides", as_path=True)
        output_dir = spipe.filepath("DIR_PROC_FB_GENO", None)

        command = " ".join([script_path, gname, guides, output_dir, par_path])

        ok, callout, ex_story = utils.call_exec(command, verb=True)
        # Story to log only (call verb=True already prints outcome)
        spipe.report_run_story(ex_story, to_log=True)

    # If script worked, update with new genome
    if ok:
        # Report numbers for targeted or not...
        spipe._set_par_val("genome_dir", spipe.filepath("DIR_PROC_FB_GENO", None))
        spipe._set_up_genome()
        gene_info = spipe.get_gene_info()
        all_genes = gene_info["all_genes"]
        targ_genes = gene_info.get("target_genes", [])
        n_genes = len(all_genes)
        n_targs = len(targ_genes)
        story = f"Focal barcoding / crispr with {n_genes} genes, {n_targs} target genes"
        spipe.report_run_story(story)
    else:
        spipe.set_problem("Failed to create focal/crispr STAR db")

    return ok


def mkref_prep_db(spipe):
    """Prepare focal/crispr mkref database

    This is the top level call from mkref for focal/crispr use case
    Handles both STAR and direct mapping cases

    Return status
    """
    ok = True

    # If using STAR, make gfasta from guides
    if spipe.get_par_val("crsp_use_star", as_bool=True):
        spipe.report_run_story("Focal/crispr STAR; mkref prep")

        # load and set guide info struct
        guide_info = load_guide_info(spipe)
        if not guide_info:
            spipe.set_problem("No fb/crispr guides for fasta/genome")
        spipe._set_guide_info(guide_info)

        guide_df = guide_info["guide_df"]

        gname = get_genome_name(spipe)
        fasta_fname = spipe.filepath("PF_CRSP_FASTA", None)
        n_rec = make_guide_fasta(guide_df, fasta_fname)
        spipe.report_run_story(f"Wrote fb/crispr fasta with {n_rec} records")
        if n_rec:
            # Set gfasta params. This is a list of lists; set just this one
            spipe._set_par_val("gfasta", [[gname, fasta_fname]])
        else:
            spipe.set_problem("No fb/crispr fasta records saved")
            ok = False

    # Direct mapping
    else:
        spipe.report_run_story("Focal/crispr direct mapping; No mkref steps")
        ok = False

    return ok


def make_guide_fasta(guide_df, fasta_fname):
    """Make on-the-fly guide fasta file for crispr STAR db

    Return number of fasta records written
    """
    head_list = list(">" + guide_df["guide_name"])
    seq_list = list(guide_df["full_seq"])

    n_rec = 0
    with open(fasta_fname, "w") as OUTFILE:
        for head, seq in zip(head_list, seq_list):
            print(head, file=OUTFILE)
            print(seq, file=OUTFILE)
        n_rec = len(head_list)
    return n_rec


# ------------------- guide 'genome' direct mapping -------------------------


def get_guide_bam_header(spipe):
    """Get bam header info for guides

    Return json like str
    """
    guide_info = spipe.get_guide_info()
    guide_df = guide_info["guide_df"]
    name_list = list(guide_df["gpx_guide_name"])
    seq_list = list(guide_df["guide_sequence"])

    # Program name, version, command line
    pname = "split-pipe"
    version = str(spipe.get_version(numonly=True))
    comline = spipe.get_command_line(wrap=False)

    # Header version and sort order
    header = {
        "HD": {"VN": version, "SO": "unsorted"},
    }

    # Each guide = reference chrom; Name and length
    ch_list = []
    for name, seq in zip(name_list, seq_list):
        ch_rec = {"LN": len(seq), "SN": name}
        ch_list.append(ch_rec)
    header["SQ"] = ch_list

    # Program info
    header["PG"] = [{"ID": pname, "PN": pname, "VN": version, "CL": comline}]

    return header


def g2cell_dict_from_df(df, g_col="guide", c_col="bc_wells"):
    """Get guide / gene to cell mapping dict from dataframe

    Return dict
    """
    # If guides are index, create col with these
    add_col = False
    if df.index.name == g_col:
        df[g_col] = df.index
        add_col = True

    new_dict = {}
    for gene in sorted(df[g_col].unique()):
        new_dict[gene] = sorted(list(df[df[g_col] == gene][c_col]))

    # If added col, drop it
    if add_col:
        df.drop(columns=[g_col], inplace=True)

    return new_dict


# Feedback frequency
DMAP_UPDATE = 200000


def map_fastq_guides(spipe, max_reads=0):
    """Map guides in fastq via direct match

    Save mappings to bam, and stats to file

    Return number of matches
    """
    spipe.report_run_story("Starting guide direct mapping")

    guide_info = spipe.get_guide_info()
    guide_seq_to_index = guide_info["guide_seq_to_index"]
    # Lists sorted large to small, unique guide lengths
    guide_end_list = guide_info["guide_end_list"]
    guide_len_list = guide_info["guide_len_list"]
    guide_mm_list = guide_info["guide_mm_list"]
    # Single start for everything; Initial end = max end = first in list
    g_start = guide_info["guide_start"]
    cig_pref = f"{g_start}S"
    max_g_end = guide_end_list[0]

    story = f"Read guide coords: start={g_start} end={max_g_end} ({guide_end_list})"
    spipe.report_run_story2(story)

    # Input and output files
    bch_fq = spipe.filepath("PF_FASTQ_BC", None)
    bam_fname = spipe.filepath("PF_BAM_RAW", None)

    # Define a minimal header for the BAM file (you may need to customize this)
    header = get_guide_bam_header(spipe)

    a = pysam.AlignedSegment()
    n_reads = n_match = 0

    # Output is a BAM file object
    with gzip.open(bch_fq, "rb") as INFILE, pysam.AlignmentFile(
        bam_fname, "wb", header=header
    ) as OFILE:
        while True:
            n_reads += 1
            if max_reads and (n_reads >= max_reads):
                break

            if (n_reads % DMAP_UPDATE) == 0:
                per = utils.report_percent_str(n_match, den=n_reads)
                spipe.report_run_story3(
                    f"Processed {n_reads} fastq records; Match guide {per}"
                )

            # Each file, header, seq, strand, and quality lines
            head = INFILE.readline()
            seql = INFILE.readline()
            _ = INFILE.readline()
            qual = INFILE.readline()

            if not seql:
                break

            header = head.decode().strip()
            seq = seql.decode().strip()
            # strand = strl.decode().strip()
            quality = qual.decode().strip()

            # Get starting guide window; Max of all guide lengths
            rguide = seq[g_start:max_g_end]

            # Each length has map dict; Large to small
            for i, g_len in enumerate(guide_len_list):
                # Trim to length
                rguide = rguide[0:g_len]
                gmatch = guide_mm_list[i].get(rguide)
                if gmatch:
                    g_end = guide_end_list[i]
                    break
            if not gmatch:
                continue

            n_match += 1

            # Unpack
            perfect_guide, num_mm, cig_str = gmatch

            a.query_name = header[1:]
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array(quality)
            a.reference_start = 0
            a.mapping_quality = 255
            a.cigarstring = cig_pref + cig_str + f"{len(seq)-g_end}S"
            a.reference_id = guide_seq_to_index[perfect_guide]

            # Tags like STAR
            # NH = number of hits (always 1)
            # NM = edit dist
            # HI = query hit index (always 1, as 1 hit)
            # AS = alignment score
            """
            tag_list = [
                ("NH", 1),
                ("HI", 1),
                ("AS", len(perfect_guide) - int(num_mm)),
                ("NM", int(num_mm)),
            ]
            a.set_tags(tag_list)
            """
            a.set_tag("NH", 1)
            a.set_tag("HI", 1)
            a.set_tag("AS", len(perfect_guide) - int(num_mm))
            a.set_tag("NM", int(num_mm))

            # gene_id = guide_seq_to_gene_id[gmatch[0]]
            # gene_name = guide_seq_to_gene[gmatch[0]]
            # geno_region = 1

            """
            set_read_tags(
                a,
                gene_id,
                gene_name,
                geno_region,
                bc_wind,
                bc_bcis,
                bc_seqs,
                polyN_seq,
                samp_name,
            )
            """

            OFILE.write(a)

    per = utils.report_percent_str(n_match, den=n_reads)
    spipe.report_run_story3(
        f"Alignment result: {n_match} reads of {n_reads} map: {per}"
    )

    spipe.report_run_story(f"Mapping saved to {bam_fname}")
    spipe.report_run_story(f"Total reads {n_reads}, guide matches {n_match}")
    update_pipe_stats(spipe, n_reads, n_match)

    return n_match


def update_pipe_stats(spipe, n_reads, n_match):
    """Update mapping stats to file"""
    pipe_df = spipe.read_csv("PF_STAT_PIPE", names="statistic,value", verb=False)
    # Assume perfect case with direct mapping
    pipe_df.loc["reads_align_input", "value"] = n_reads
    pipe_df.loc["reads_align_unique", "value"] = n_match
    pipe_df.loc["reads_align_multimap", "value"] = 0
    pipe_df.loc["reads_too_many_loci", "value"] = 0
    pipe_df.index.name = "statistic"
    # Dump all lines to log; n_rows=0
    spipe.write_df(pipe_df, "PF_STAT_PIPE", n_rows=0)


# ------------------- DGE outputs -------------------------------------------


def write_guide_cell_map_files(spipe, samp, dge_var_dict, thresh):
    """Write guide to cell mapping files for focal / CRISPR

    dge_var_dict = dict with dge stuff; barcode indexes, cells, matrix

    Returns nothing
    """
    spipe.report_run_story(f"Saving guide cell mapping and counts; Threshold {thresh}")

    # Matrix, cells and genes
    mtx = dge_var_dict["dge_matrix"]
    bc_wells = dge_var_dict["bc_indexes"]
    gene_names = spipe.get_gene_list()

    # Collect guide-cell mapping and counts
    df_list = []
    cnt_list = []
    for i, gene in enumerate(gene_names):
        # Single guide
        df = pd.DataFrame.sparse.from_spmatrix(
            mtx[:, i], index=bc_wells, columns=[gene]
        )
        # Reset index so cells ('bc_wells') is a column
        df = df[df[gene] >= thresh].reset_index()
        df["guide"] = gene
        df_list.append(df[["guide", "index"]])
        cnt_list.append([gene, len(df)])

    # Save as dataframes
    map_df = pd.concat(df_list, axis=0)
    map_df.columns = ["guide", "bc_wells"]
    map_df.set_index("guide", inplace=True)
    fkey = (
        "SFR_FB_CELL_MAP_TAB"
        if spipe.is_focal_crispr(focal=True)
        else "SFR_CRSP_CELL_MAP_TAB"
    )
    spipe.write_df(map_df, fkey, samp=samp)

    # Also save as json guide-to-cell dict
    map_dict = g2cell_dict_from_df(map_df)
    fkey = (
        "SFR_FB_CELL_MAP_DIC"
        if spipe.is_focal_crispr(focal=True)
        else "SFR_CRSP_CELL_MAP_DIC"
    )
    fname = spipe.filepath(fkey, samp)
    utils.write_json(map_dict, fname)
    spipe.report_run_story(f"Wrote guide-to-cell dictionary {fname}")

    # Simple count table
    cnt_df = pd.DataFrame(cnt_list, columns=["guide", "count"])
    cnt_df["count"] = cnt_df["count"].astype(int)
    cnt_df.set_index("guide", inplace=True)
    fkey = "SFR_FB_CELL_CT" if spipe.is_focal_crispr(focal=True) else "SFR_CRSP_CELL_CT"
    spipe.write_df(cnt_df, fkey, samp=samp)


#### ---------------- transcript and guide filtering  -------------------------------------


def get_up_score(lst, minima):
    scr = 0
    previous = 0
    for i in range(len(lst)):
        if i == 0:
            previous = lst[i]
        else:
            if lst[i] > previous and lst[i] > minima:
                scr += 1
            previous = lst[i]
    return scr


def get_threshold_from_multimodal(
    indf, feature="read", rng=5, rt_default=3, tscp_default=1.5
):
    # test = tscp_grna[tscp_grna["count"]>1]
    if feature == "read":
        hist_values = indf.groupby(["count"]).size().index.tolist()
        hist_count = indf.groupby(["count"]).size().tolist()
        min_dict = {}
        for i in range(len(hist_count)):
            if i < rng:
                compare_lst = hist_count[0 : i + rng]
            else:
                compare_lst = hist_count[i - rng : i + rng]
            if min(compare_lst) == hist_count[i]:
                minima = hist_count[i]
                if (sum(hist_count[0:i]) / sum(hist_count)) < 0.9:
                    score = get_up_score(hist_count[i : len(hist_count)], minima)
                    min_dict[hist_values[i]] = score

        if len(min_dict) == 0:
            local_minima = rt_default
        else:
            local_minima = max(min_dict, key=min_dict.get)

    elif feature == "tscp":
        # hist_values = indf.groupby(["log_tscp_cnt"]).size().index.tolist()
        # hist_count = indf.groupby(["log_tscp_cnt"]).size().tolist()
        hist_values = np.histogram(indf["log_tscp_cnt"], bins="fd")[1][1:]
        hist_count = np.histogram(indf["log_tscp_cnt"], bins="fd")[0]
        hist_values = hist_values[hist_count != 0]
        hist_count = hist_count[hist_count != 0]
        min_dict = {}
        for i in range(len(hist_count)):
            if i < rng:
                compare_lst = hist_count[0 : i + rng]
            else:
                compare_lst = hist_count[i - rng : i + rng]
            if min(compare_lst) == hist_count[i]:
                minima = hist_count[i]
                if (sum(hist_count[0:i]) / sum(hist_count)) < 0.9:
                    score = get_up_score(hist_count[i : len(hist_count)], minima)
                    min_dict[hist_values[i - 1]] = score

        if len(min_dict) == 0:
            local_minima = tscp_default
        else:
            local_minima = max(min_dict, key=min_dict.get)
    return local_minima, min_dict


def crispr_get_cutoff(
    tscp_cnt,
    grnas=[],
    rt_default=3,
    tscp_default=1,
    manual=False,
    manual_read_cutoff=1,
    manual_tscp_cutoff=1,
):
    if len(grnas) > 0:
        tscp_grna = tscp_cnt[tscp_cnt.gene_name.isin(grnas)]
    else:
        tscp_grna = tscp_cnt
    # implement minima algo to get read cutoff
    if manual is True:
        read_cutoff = manual_read_cutoff
    else:
        read_cutoff, read_cutoff_dict = get_threshold_from_multimodal(
            tscp_grna, feature="read", rt_default=rt_default
        )

    tscp_read_filtered = tscp_grna[tscp_grna["count"] >= read_cutoff]
    # implement minima algo to get transcript cutoff
    tscp_read_filtered["bc_wells_gene"] = (
        tscp_read_filtered["bc_wells"] + "___" + tscp_read_filtered["gene"]
    )
    tscp_indf = tscp_read_filtered.groupby("bc_wells_gene").agg({"polyN": ["size"]})
    tscp_indf.columns = ["tscp_cnt"]
    tscp_indf["log_tscp_cnt"] = np.log2(tscp_indf.tscp_cnt)
    if manual is True:
        tscp_cutoff = np.log2(manual_tscp_cutoff)
    else:
        tscp_cutoff, tscp_cutoff_dict = get_threshold_from_multimodal(
            tscp_indf, feature="tscp", tscp_default=tscp_default
        )
    # tscp_indf_filtered = tscp_indf[tscp_indf["log_tscp_cnt"] >= tscp_cutoff]
    # tscp_file_grna = tscp_filtered.set_index(['bc_wells', 'gene_name'])
    # filtered_tscp_file_grna = tscp_file_grna.loc[tscp_indf_filtered.index.tolist()]
    mask_set = set(tscp_indf[tscp_indf["log_tscp_cnt"] >= tscp_cutoff].index)
    filtered_tscp_file = tscp_read_filtered[
        tscp_read_filtered["bc_wells_gene"].isin(mask_set)
    ]

    if manual is True:
        tscp_return = manual_tscp_cutoff
    else:
        tscp_return = 2**tscp_cutoff

    return filtered_tscp_file, read_cutoff, tscp_return
