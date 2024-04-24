#!/usr/bin/env python
# coding: utf-8

import os
import sys
import multiprocessing as mp
from multiprocessing import Process
import pandas as pd
from collections import defaultdict
import numpy as np
import subprocess
import copy
from collections import Counter
import time
import glob
import re
import argparse
from itertools import compress


# 2023-04-17 RTK; Bumped version from Ashwath's original
# 2023-11-11 RTK; Bumped verion from Ashwath's (v1.3 >--> v1.4)
# 2023-12-07 RTK; v1.5 Rename sample .txt outputs to .csv
VER_NUM = "v1.5"

VERSION = os.path.basename(__file__) + " " + VER_NUM


def main():
    parser = argparse.ArgumentParser(
        description="Processes TCR barcode_head.fastq.gz file"
    )
    parser.add_argument(
        "--wkdir", help="Path to sample directory within 'alignment' folder"
    )
    parser.add_argument(
        "--igblastdir",
        help="Path to igblast directory. This is the directory where bin,tr_fasta,tr_database folders are located",
    )
    parser.add_argument(
        "--genome",
        help="Genome to use as a reference. It has to be either human or mouse.",
    )
    parser.add_argument(
        "--kmer",
        type=int,
        default=21,
        help="Length of the kmer to use in TCR assembly (default: 21)",
    )
    parser.add_argument(
        "--read_threshold",
        type=int,
        default=2,
        help="Minimum number of reads required for TCR assembly (default: 2)",
    )
    parser.add_argument(
        "--trim",
        type=int,
        default=30,
        help="Number of bases to trim from reads before TCR assembly (default: 30)",
    )
    parser.add_argument(
        "--gfile",
        help="Get sample group specifications from well name file; See --explain_gfile for format",
    )
    parser.add_argument(
        "--wfile",
        help="Get sample group specifications from barcode list file. It is a tab separated file with sample name in first column and comma separated barcode indices in second column.",
    )
    parser.add_argument(
        "--build_contigs_at_cell_level",
        action="store_true",
        help="Trigger TCR assembly at the cell level. This will not consider transcripts when building contigs. By default, this is 'off' which means transcript information is used in TCR assembly.(default: off)",
    )
    parser.add_argument(
        "--cell_file",
        help="Path to the transcriptome cell metadata file. This option triggers the filtering of cells that are not present in the transcriptome data. If this path is not given, filtered files will not be generated",
    )
    parser.add_argument(
        "--sublibraries",
        help="File with paths to output directories of each sublibrary (Combine mode only). This also triggers combine mode. ",
    )
    parser.add_argument("--output_dir", help="Path to out directory for combine mode.")
    parser.add_argument(
        "--nthreads", type=int, default=1, help="Number of threads to use (default = 1)"
    )
    parser.add_argument(
        "--keep_temp_files",
        action="store_true",
        help="Does not delete temporary files generated from the pipeline. Warning: Can take up a lot of space. (default : off)",
    )
    parser.add_argument(
        "--explain_gfile",
        default=0,
        action="store_true",
        help="Explain the details about the gfile",
    )
    parser.add_argument("-V", "--version", action="version", version=VERSION)
    args = parser.parse_args()

    # argparse
    if len(sys.argv) < 2:
        print("No arguments given")
        return True
    if args.explain_gfile:
        explain_gfile_details()
        return True
    nthreads = args.nthreads
    sampfile = ""
    sample_flag = 0
    sample_file_type = ""
    if not args.gfile:
        if not args.wfile:
            sample_flag = 0
        else:
            sampfile = check_infile(args.wfile, verb=True)
            sample_flag = 1
            sample_file_type = "wfile"
    else:
        sampfile = check_infile(args.gfile, verb=True)
        sample_flag = 1
        sample_file_type = "gfile"

    if args.sublibraries:
        if args.output_dir:
            print("Running combine mode")
            tmp = run_combine(
                args.sublibraries,
                args.output_dir,
                sample_flag,
                sampfile,
                sample_file_type,
                nthreads,
            )
            return True
        else:
            print("No output directory provided for combine mode")
            return True

    if not args.wkdir:
        print("Working directory not provided")
        return True
    if not args.igblastdir:
        print("Igblast directory not provided")
        return True
    if not args.genome:
        print("No genome provided")
        return True

    # get arguments
    wkdir = args.wkdir
    genome = args.genome
    v_seqment_fasta = args.igblastdir + "/tr_fasta/" + genome + "_V"
    j_seqment_fasta = args.igblastdir + "/tr_fasta/" + genome + "_J"
    polyN_contig_flag = 1
    if args.build_contigs_at_cell_level:
        polyN_contig_flag = 0
    read_threshold = args.read_threshold
    igblast_dir = args.igblastdir
    trim = args.trim
    kmer_len = args.kmer
    cleanup = 1
    if args.keep_temp_files:
        cleanup = 0

    # Check if barcode file has entries
    chk_barcode_file = os.path.isfile(wkdir + "/process/barcode_head.fastq.gz")
    if chk_barcode_file:
        chk_barcode_flag = (
            os.stat(wkdir + "/process/barcode_head.fastq.gz").st_size == 0
        )
        if chk_barcode_flag:
            print("barcode_head.fastq.gz is empty. Stopping run.")
            return True
    else:
        print("barcode_head.fastq.gz doesnt exist. Stopping run.")
        return True
    # set IGDATA if its not set for igblast
    if os.getenv("IGDATA") is None:
        print("IGDATA not set, setting it to " + igblast_dir)
        os.environ["IGDATA"] = igblast_dir
    elif os.path.samefile(os.getenv("IGDATA"), igblast_dir) == False:
        print("IGDATA not set to igblast directory. Setting it to " + igblast_dir)
        os.environ["IGDATA"] = igblast_dir

    # set TMPDIR
    tmpdir = wkdir + "/tmp/"
    mk_cmd = "mkdir " + tmpdir
    c = subprocess.getstatusoutput(mk_cmd)
    print("Setting TMPDIR to " + tmpdir)
    os.environ["TMPDIR"] = tmpdir

    # Check parent directory
    transcriptome_file = ""
    tscp_flag = 0
    if args.cell_file:
        tscp_flag = 1
        print(
            "Transcriptome file provided. Filtered output folder with files with cells from transcriptome will be generated."
        )
        transcriptome_file = args.cell_file
    else:
        print(
            "Transcriptome file not provided. No filtered output folder will be generated."
        )

    if sample_flag == 0:
        print("No sample group file provided")
        print("No sample level summary will be generated")

    # Run TCR analysis
    mk_cmd = "mkdir " + wkdir + "/process/tcr/"
    c = subprocess.getstatusoutput(mk_cmd)
    print("Sorting input fastq file")
    c = sort_fastq(
        infastq=wkdir + "/process/barcode_head.fastq.gz",
        ofastq=wkdir + "/process/tcr/barcode_head_sorted.fastq",
        odir=wkdir + "/process/tcr/",
    )

    # Split fastq
    print("Splitting input fastq file")
    filenames, total_reads, all_reads_per_cell = split_fq_new(
        nthreads=nthreads,
        infastq=wkdir + "/process/tcr/barcode_head_sorted.fastq",
        odir=wkdir + "/process/tcr",
    )

    # PolyN correction, assembly and annotation
    Pros = []
    print("Starting TCR assembly")

    for i in range(nthreads):
        p = Process(
            target=run_tcr,
            args=(
                i,
                wkdir + "/process/tcr",
                filenames[i],
                genome,
                v_seqment_fasta,
                j_seqment_fasta,
                polyN_contig_flag,
                read_threshold,
                igblast_dir,
                kmer_len,
                trim,
            ),
        )
        Pros.append(p)
        p.start()

    for t in Pros:
        t.join()

    # combine igblast outputs
    igblast_all, igblast_condensed = igblast_combine(wkdir=wkdir + "/process/tcr")

    # Generate summary files
    mk_cmd = "mkdir -p " + wkdir + "/process/scratch/TCR_unfiltered/"
    c = subprocess.getstatusoutput(mk_cmd)

    print("Generating sublibrary unfiltered summary files")
    prod_dub_info = generate_summary_files_threaded(
        igblast_condensed,
        opfix=wkdir + "/process/scratch/TCR_unfiltered/",
        polyN_contig_flag=polyN_contig_flag,
        nthreads=nthreads,
        proc_dir=wkdir + "/process/tcr",
    )

    print("Generating sublibrary unfiltered metrics file")
    summary_df = make_summary_file(
        indir=wkdir + "/process/tcr",
        total_reads=total_reads,
        total_reads_per_cell=all_reads_per_cell,
        prod_dub_info=prod_dub_info,
    )

    tmp = generate_metrics_file(
        total_reads,
        summary_df,
        proc_dir=wkdir + "/process",
        opfix=wkdir + "/process/scratch/TCR_unfiltered/tcr_unfilt",
        isSample=False,
        filtered=False,
    )

    # Generate sample level files
    if sample_flag == 1:
        print("Generating unfiltered sample summary files")
        sample_dict = get_sample_wells(sampfile, intype=sample_file_type)
        tmp = get_sample_level_files(
            igblast_condensed,
            sample_dict,
            main_odir=wkdir,
            polyN_contig_flag=polyN_contig_flag,
            summary_df=summary_df,
            total_reads=total_reads,
            filtered=False,
            nthreads=nthreads,
        )

    # Generate filtered files if needed
    if tscp_flag == 1:
        print("Generating sublibrary filtered summary files")
        mk_cmd = "mkdir -p " + wkdir + "/process/scratch/TCR_filtered/"
        c = subprocess.getstatusoutput(mk_cmd)

        igblast_condensed_filtered, pass_cells = get_tscp_cells_only(
            igblast_condensed, transcriptome_file
        )
        prod_dub_info = generate_summary_files_threaded(
            igblast_condensed_filtered,
            opfix=wkdir + "/process/scratch/TCR_filtered/",
            polyN_contig_flag=polyN_contig_flag,
            nthreads=nthreads,
            proc_dir=wkdir + "/process/tcr",
        )

        # Generate metrics for filtered.
        # pass_cells = list(set(igblast_condensed_filtered['bc_wind'].tolist()))
        print("Generating sublibrary filtered cells metrics file")
        tmp = generate_metrics_file(
            total_reads,
            summary_df,
            proc_dir=wkdir + "/process",
            opfix=wkdir + "/process/scratch/TCR_filtered/tcr_filt",
            pass_cells=pass_cells,
            isSample=False,
            filtered=True,
        )

        if sample_flag == 1:
            print("Generating filtered sample summary files")
            sample_dict = get_sample_wells(sampfile, intype=sample_file_type)
            tmp = get_sample_level_files(
                igblast_condensed_filtered,
                sample_dict,
                main_odir=wkdir,
                polyN_contig_flag=polyN_contig_flag,
                summary_df=summary_df,
                total_reads=total_reads,
                pass_cells=pass_cells,
                filtered=True,
                nthreads=nthreads,
            )

    if cleanup == 1:
        print("Cleaning up tmp directories")
        rm_cmd = "rm -r " + wkdir + "/process/tcr/"
        c = subprocess.getstatusoutput(rm_cmd)
        rm_cmd = "rm -r " + tmpdir
        c = subprocess.getstatusoutput(rm_cmd)


###############################################################
# Functions from split-pipe
# Well index parts
def fqrec_to_bc_winds(r):
    return r.split("__")[0]


def fqrec_to_bc1_wind(r):
    return fqrec_to_bc_winds(r).split("_")[0]


def fqrec_to_bc2_wind(r):
    return fqrec_to_bc_winds(r).split("_")[1]


def fqrec_to_bc3_wind(r):
    return fqrec_to_bc_winds(r).split("_")[2]


def fqrec_to_ptype(r):
    return r.split("__")[1]


# Barcode index parts
def fqrec_to_bc_bcis(r):
    return r.split("__")[2]


def fqrec_to_bc1_bci(r):
    return fqrec_to_bc_bcis(r).split("_")[0]


def fqrec_to_bc2_bci(r):
    return fqrec_to_bc_bcis(r).split("_")[1]


def fqrec_to_bc3_bci(r):
    return fqrec_to_bc_bcis(r).split("_")[2]


# Barcode sequence parts
def fqrec_to_bc_seqs(r):
    return r.split("__")[3]


def fqrec_to_bc1_seq(r):
    return r.split("__")[3].split("_")[0]


def fqrec_to_bc2_seq(r):
    return r.split("__")[3].split("_")[1]


def fqrec_to_bc3_seq(r):
    return r.split("__")[3].split("_")[2]


# Poly-N part
def fqrec_to_polyN(r):
    return r.split("__")[4]


###############################################################
##Splitting fastq functions
# Get chunk filepath for fastq splitting.
def chunk_fastq_filepath(fastq=False, chunk=None):
    if fastq:
        fq_fname = "barcode_head_sorted.fastq"  # hard coded standard output from splitpipe pre processing.
        if chunk is None:
            fname = fq_fname
        else:
            fname = f"{fq_fname.split('.fastq')[0]}.chunk{chunk}.fastq"
            fpath = fname
    else:
        story = "chunk_fastq_filepath called with fastq set to false"
        raise ValueError(story)
    return fpath


# Functions for sorting and splitting barcode_head.fastq.gz file.
# Sorting function
def sort_fastq(
    infastq="barcode_head.fastq.gz", ofastq="barcode_head_sorted.fastq", odir="./"
):
    # sort_cmd = """zcat {0} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {1}""".format(infastq,ofastq)
    # step by step is somehow faster than piping.
    zcat_tmp = odir + "/zcat_tmp.txt"
    sort_tmp = odir + "/sort_tmp.txt"

    zcat_cmd = """zcat {0} | paste - - - - > {1}""".format(infastq, zcat_tmp)
    c = subprocess.getstatusoutput(zcat_cmd)
    sort_cmd = """sort -k1,1 -t " " {0} > {1}""".format(zcat_tmp, sort_tmp)
    c = subprocess.getstatusoutput(sort_cmd)
    tr_cmd = """cat {0} | tr "\t" "\n" > {1}""".format(sort_tmp, ofastq)
    c = subprocess.getstatusoutput(tr_cmd)
    cleanup_cmd = """rm {0} {1}""".format(zcat_tmp, sort_tmp)
    c = subprocess.getstatusoutput(cleanup_cmd)


# Helper
def file_len(filename):
    i = -1
    with open(filename) as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


# Split fastq into "n" fastq files. This is allow for parallel processing.
def split_fq_new(nthreads, infastq="barcode_head_sorted.fastq", odir="./"):
    # Init outputs for each thread
    fq_chunks = {}
    filenames = []
    for i in range(nthreads):
        chunk_fname = chunk_fastq_filepath(fastq=True, chunk=(i + 1))
        filenames.append(odir + "/" + chunk_fname)
        fq_chunks[i] = open(odir + "/" + chunk_fname, "w")
    # Get the total reads
    total_reads = file_len(infastq) / 4
    reads_per_chunk = int(np.ceil(total_reads / nthreads))
    print(f"Reads per chunk {reads_per_chunk}")
    c = 0
    prev_cell_barcode = ""
    reads = []
    reads_seq = []
    reads_qual = []
    total_reads_per_cell = {}  # track number of reads per cell for metrics
    # Open file
    fileHandler = open(infastq, "r")
    while True:
        # Get next line from file
        aline = fileHandler.readline()
        # If line is empty then end of file reached
        if not aline:
            break
        aline = aline.rstrip()
        reads.append(aline[1:])
        cell_barcode = fqrec_to_bc_winds(aline[1:])
        read_sq = fileHandler.readline()
        reads_seq.append(read_sq.rstrip())
        tmp = fileHandler.readline()
        read_qual = fileHandler.readline()
        reads_qual.append(read_qual.rstrip())
        if cell_barcode != prev_cell_barcode:
            chunk = int(np.floor(c / reads_per_chunk))
            # Check to keep last chunk right
            if prev_cell_barcode in total_reads_per_cell:
                print(
                    "Duplicate cell entry when counting reads per cell in split_fq_new "
                    + prev_cell_barcode
                )
            elif prev_cell_barcode != "":
                total_reads_per_cell[prev_cell_barcode] = len(reads) - 1

            if chunk >= nthreads:
                chunk = nthreads - 1

            for r in range(len(reads[:-1])):
                fq_chunks[chunk].write("@" + reads[r] + "\n")
                fq_chunks[chunk].write(reads_seq[r] + "\n")
                fq_chunks[chunk].write("+\n")
                fq_chunks[chunk].write(reads_qual[r] + "\n")
            reads = reads[-1:]
            reads_seq = reads_seq[-1:]
            reads_qual = reads_qual[-1:]
        prev_cell_barcode = cell_barcode
        c += 1
    # print leftover
    if c > 0:
        if prev_cell_barcode in total_reads_per_cell:
            print(
                "Duplicate cell entry when counting reads per cell in split_fq_new "
                + prev_cell_barcode
            )
        elif prev_cell_barcode != "":
            total_reads_per_cell[prev_cell_barcode] = len(reads)
        chunk = int(np.floor(c / reads_per_chunk))
        if chunk >= nthreads:
            chunk = nthreads - 1
        for r in range(len(reads)):
            fq_chunks[chunk].write("@" + reads[r] + "\n")
            fq_chunks[chunk].write(reads_seq[r] + "\n")
            fq_chunks[chunk].write("+\n")
            fq_chunks[chunk].write(reads_qual[r] + "\n")
    for i in range(nthreads):
        fq_chunks[i].close()
    return filenames, total_reads, total_reads_per_cell


###############################################################
# # polyN correction functions

# Adapted Existing functions from splitpipe for poly N correction.

# Helper functions.
bases = list("ACGT")


def single_mut_seqs(seq, n):
    """Return a list of all sequences with a 1 nt mutation (no ins/del)."""
    mut_seqs = []
    for i in range(n):
        for b in bases:
            if b != seq[i]:
                mut_seqs.append(seq[:i] + b + seq[i + 1 :])
    return mut_seqs


def hamdist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2
    Max dif of 2
    """
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
        if diffs > 1:
            break
    return diffs


# Collapse polyN. Adapted the concept from splitpipe.
def collapse_polyN_fastq(kmer_list, counts):
    """Collapses polyN that are within 1 hamming dist of each other.
    Collapses and returns a dict with collapsed mapping

    Returns dict original polyN to new mapped polyN. Selection is done based on count.
    """
    # TODO; This could be sped up, maybe quite a bit
    kmer_len = len(kmer_list[0])
    kmer_seqs_dict = dict(zip(kmer_list, np.zeros(len(kmer_list))))
    ham_dict = defaultdict(list)
    seq_len = len(kmer_list)
    kmer_map_dict = defaultdict(list)
    # Collect lists of within-Hamming-dist seqs for each seq
    # If there are many seqs, create a dict with all one-offs
    if seq_len > 50:
        for seq in kmer_list:
            mut_seqs = single_mut_seqs(seq, kmer_len)
            for ms in mut_seqs:
                try:
                    kmer_seqs_dict[ms]
                except:
                    pass
                else:
                    ham_dict[ms].append(seq)
    # Else direct seq-to-seq comparison
    # (if only single seq, seq_len==1, this does nothing)
    else:
        for i in range(seq_len):
            for j in range(i + 1, seq_len):
                s1, s2 = kmer_list[i], kmer_list[j]
                ham_dist = hamdist(s1, s2)
                if ham_dist <= 1:
                    ham_dict[s1].append(s2)
                    ham_dict[s2].append(s1)

    # Create dict of [seq] = count
    kmer_counts = dict(zip(kmer_list, counts))
    # Now try to merge in one-off seqs (and counts) together
    for kmer in kmer_list:
        cur_count = kmer_counts[kmer]
        for hm in ham_dict[kmer]:
            # If find one-off with more counts, merge kmer into one-off and delete
            try:
                if kmer_counts[hm] > cur_count:
                    kmer_counts[hm] += cur_count
                    del kmer_counts[kmer]
                    kmer_map_dict[kmer] = hm
                    break
            except KeyError:
                pass
    return kmer_counts, kmer_map_dict


# Wrapper function for the polyN helper functions. Returns a dataframe with corrected polyN.
def get_corrected_polyN(
    bc_wind_list,
    bc_ptype_list,
    bc_bci_list,
    bc_seq_list,
    polyN_list,
    seq_list,
    qual_list,
):
    df = pd.DataFrame(
        {
            "bc_wells": bc_wind_list[:-1],
            "ptype": bc_ptype_list[:-1],
            "bc_bcis": bc_bci_list[:-1],
            "cell_barcode": bc_seq_list[:-1],
            "polyN": polyN_list[:-1],
            "read_sq": seq_list[:-1],
            "read_qual": qual_list[:-1],
        }
    )
    polyN_count_df = df.groupby(["bc_bcis", "polyN"]).size().reset_index()
    polyN_count_df.columns = ["bc_bcis", "polyN", "count"]
    polyN_dict = polyN_count_df.groupby(["bc_bcis"]).apply(
        lambda x: pd.Series(collapse_polyN_fastq(list(x["polyN"]), list(x["count"])))
    )

    df["polyN_corrected"] = None
    for bcis in set(polyN_count_df["bc_bcis"]):
        new_polyN = df.loc[df["bc_bcis"] == bcis]["polyN"].apply(
            lambda x: polyN_dict[1].loc[bcis][x] if x in polyN_dict[1].loc[bcis] else x
        )
        df.loc[df["bc_bcis"] == bcis, "polyN_corrected"] = new_polyN
    df["header"] = df[
        ["bc_wells", "ptype", "bc_bcis", "cell_barcode", "polyN_corrected"]
    ].apply("__".join, axis=1)
    return df


# Printing fastq files with corrected polyN
def print_corrected_polyN_fastq(df, ofile):
    df = df.sort_values("polyN_corrected")

    def print_row(h, s, q, ofile):
        ofile.write("@" + h + "\n" + s + "\n+\n" + q + "\n")

    print_out = [
        print_row(x[0], x[1], x[2], ofile)
        for x in zip(df["header"], df["read_sq"], df["read_qual"])
    ]


# Umbrella function to run polyN correction on each chunk fastq file. Returns corrected polyN fastq file.


def run_polyN(chunk_file):
    fqfile = open(chunk_file, "r")
    o_filename = chunk_file.split(".fastq")[0] + ".polyN.corrected.fastq"
    o_fqfile = open(o_filename, "w")  # output file
    # init
    bc_wind_list = []
    bc_ptype_list = []
    bc_bci_list = []
    bc_seq_list = []
    polyN_list = []
    seq_list = []
    qual_list = []
    cell = 0
    next_cell = False
    while True:
        aline = fqfile.readline()
        # If line is empty then end of file reached
        if not aline:
            break
        aline = aline.rstrip()
        read_name = aline[1:]
        read_sq = fqfile.readline().rstrip()
        tmp = fqfile.readline().rstrip()
        read_qual = fqfile.readline().rstrip()

        # Get barcode, polyN information
        bc_wind = fqrec_to_bc_winds(read_name)
        bc_ptype = fqrec_to_ptype(read_name)
        bc_bcis = fqrec_to_bc_bcis(read_name)
        bc_seqs = fqrec_to_bc_seqs(read_name)
        polyN_seq = fqrec_to_polyN(read_name)

        # Append to existing list
        bc_wind_list.append(bc_wind)
        bc_ptype_list.append(bc_ptype)
        bc_bci_list.append(bc_bcis)
        bc_seq_list.append(bc_seqs)
        polyN_list.append(polyN_seq)
        seq_list.append(read_sq)
        qual_list.append(read_qual)
        cell += 1

        if cell > 1:
            # Compare well indexes to check if next cell
            #  (Not bc seqs, which are uncorrected so may differ by chance)
            #  (Not bci barcode indexes, which may map to same well==cell )
            if bc_wind_list[-1] != bc_wind_list[-2]:
                next_cell = True  # Process reads from one cell

            if next_cell:
                # polyN correction and printing new fastq files.
                poly_corr_df = get_corrected_polyN(
                    bc_wind_list,
                    bc_ptype_list,
                    bc_bci_list,
                    bc_seq_list,
                    polyN_list,
                    seq_list,
                    qual_list,
                )
                print_corrected_polyN_fastq(poly_corr_df, o_fqfile)

                # Track number of reads in the cell
                num_reads = len(bc_wind_list[:-1])

                # Reset lists, keeping the first read from the new cell
                bc_wind_list = [bc_wind]
                bc_ptype_list = [bc_ptype]
                bc_bci_list = [bc_bcis]
                bc_seq_list = [bc_seqs]
                polyN_list = [polyN_seq]
                seq_list = [read_sq]
                qual_list = [read_qual]
                cell = 1
                next_cell = False

    # Capture last cell
    if cell > 1:
        # Adding tmp so the function doesnt miss the last read.
        bc_wind_list.append("tmp")
        bc_ptype_list.append("tmp")
        bc_bci_list.append("tmp")
        bc_seq_list.append("tmp")
        polyN_list.append("tmp")
        seq_list.append("tmp")
        qual_list.append("tmp")
        # polyN correction and printing new fastq files.
        poly_corr_df = get_corrected_polyN(
            bc_wind_list,
            bc_ptype_list,
            bc_bci_list,
            bc_seq_list,
            polyN_list,
            seq_list,
            qual_list,
        )
        print_corrected_polyN_fastq(poly_corr_df, o_fqfile)
        num_reads = len(bc_wind_list[:-1])

    fqfile.close()
    o_fqfile.close()
    return o_filename


###############################################################
# # Contig formation functions


# Function to get edges and read information from a sequence
def get_edgeinfo_from_sequence(edges, edges_r, rcnt, sequence, k=3):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence. Also returns a dictionary with
    all edges.
    """
    read_only_edges = []
    # tmp_done = {}
    seqlen = len(sequence)
    for i in range(0, len(sequence)):
        kmer = sequence[i : i + k]
        next_kmer = sequence[i + 1 : i + k + 1]
        if len(kmer) != k:
            continue
        if len(next_kmer) != k:
            continue
        edge_key = kmer + ":" + next_kmer
        if edge_key not in read_only_edges:
            read_only_edges.append(edge_key)
            if edge_key in edges:
                edges[edge_key] += 1
                edges_r[edge_key].append((rcnt, seqlen))
            else:
                edges[edge_key] = 1
                edges_r[edge_key] = [(rcnt, seqlen)]

    return edges, edges_r, read_only_edges


# Get only kmers and edges
def get_kmer_and_edges_from_sequence(kmers, edges, sequence, k=3):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence. Also returns a dictionary with
    all edges.
    """
    tmp_done = (
        {}
    )  # Using temp dict to track edges. This prevents repeat regions getting very high edge support
    # from one read. This makes it so that each edge can only get one support from one read.
    for i in range(0, len(sequence)):
        kmer = sequence[i : i + k]
        next_kmer = sequence[i + 1 : i + k + 1]
        if len(kmer) != k:
            continue

        # count occurrence of this kmer in sequence
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
        if len(next_kmer) != k:
            continue
        edge_key = kmer + ":" + next_kmer
        if edge_key not in tmp_done:
            tmp_done[edge_key] = 1
            if edge_key in edges:
                edges[edge_key] += 1
            else:
                edges[edge_key] = 1

    return kmers, edges


# Helper functions for contig extension from debruijin graph.
def get_fw_info_greedy(kmer, edges):
    bases = ["A", "T", "G", "C"]
    next_kmers = [kmer[1:] + base for base in bases]
    fw_list = []
    fw_vals = []
    for km in next_kmers:
        key = kmer + ":" + km
        if key in edges:
            if len(fw_list) == 0:
                fw_list.append(key)
                fw_vals.append(edges[key])
            else:
                if edges[key] > fw_vals[0]:
                    fw_list[0] = key
                    fw_vals[0] = edges[key]
    return fw_list, fw_vals


def get_bw_info_greedy(kmer, edges):
    bases = ["A", "T", "G", "C"]
    prev_kmers = [base + kmer[:-1] for base in bases]
    bw_list = []
    bw_vals = []
    for km in prev_kmers:
        key = km + ":" + kmer
        if key in edges:
            if len(bw_list) == 0:
                bw_list.append(key)
                bw_vals.append(edges[key])
            else:
                if edges[key] > bw_vals[0]:
                    bw_list[0] = key
                    bw_vals[0] = edges[key]
    return bw_list, bw_vals


# Contig forward function
def get_fw_contig_greedy(contig, edge, edges, idx, done, min_edge_support=0):
    kmer1 = edge.split(":")[0]
    kmer2 = edge.split(":")[1]
    fw_list, fw_vals = get_fw_info_greedy(kmer2, edges)
    contig[idx].append(kmer2[-1])
    tmp = contig[idx].copy()
    if edge in done:
        return contig
    else:
        done[edge] = 1
        if len(fw_list) == 0:
            return contig
        if len(fw_list) == 1:
            # Needs testing
            if fw_vals[0] > min_edge_support:
                get_fw_contig_greedy(contig, fw_list[0], edges, idx, done)
            else:
                return contig
        if len(fw_list) > 1:
            print("Error more than one edge returned")


# Contig reverse function
def get_bw_contig_greedy(contig, edge, edges, idx, done, min_edge_support=0):
    kmer1 = edge.split(":")[0]
    kmer2 = edge.split(":")[1]
    bw_list, bw_vals = get_bw_info_greedy(kmer1, edges)
    contig[idx].append(kmer1[0])
    tmp = contig[idx].copy()
    if edge in done:
        return contig
    else:
        done[edge] = 1
        if len(bw_list) == 0:
            return contig
        if len(bw_list) == 1:
            # NOTE: Needs testing
            if bw_vals[0] > min_edge_support:
                get_bw_contig_greedy(contig, bw_list[0], edges, idx, done)
            else:
                return contig
        if len(bw_list) > 1:
            print("Error more than one edge returned")


# Function to get contig support. Easier than adapting into the recursive function. It
# calculates read support by using the overlap threshold. Only reads that have > threshold are
# counted as a support for the contig. Hence, it is possible to get zero support for some contigs.


def get_contig_support(contig, edges_r, edges, k, read_overlap_threshold=0.9):
    edge_scores = []
    all_edges = []
    for i in range(0, len(contig)):
        kmer = contig[i : i + k]
        next_kmer = contig[i + 1 : i + k + 1]
        if len(kmer) != k:
            continue
        if len(next_kmer) != k:
            continue
        edge_key = kmer + ":" + next_kmer
        if edge_key in edges_r:
            edge_scores.append(edges[edge_key])
            all_edges.extend(edges_r[edge_key])
            # Checks. Needs cleanup.
            tt = len(list(set(edges_r[edge_key])))
            if tt < len(edges_r[edge_key]):
                print(
                    "Read represented multiple times in an edge in get_contig_support"
                )
            if len(edges_r[edge_key]) != edges[edge_key]:
                print("Disagreement between edges and edges_r in get_contig_support")
        else:
            print(
                "Error " + edge_key + " not found in dictionary in get_contig_support"
            )

    set_all_edges = Counter(all_edges)
    read_cnt = 0
    used_reads = []
    for rd in set_all_edges:
        if (set_all_edges[rd] + k) / rd[1] >= read_overlap_threshold:
            read_cnt += 1
            used_reads.append(rd[0])
    max_score = max(edge_scores)
    min_score = min(edge_scores)
    avg_score = sum(edge_scores) / len(edge_scores)
    return max_score, min_score, avg_score, read_cnt, used_reads


# Wrapper function to get support information. Also prints out the fasta file with the contigs.
def get_support_wrapper(cell_contigs, cell_info, edges_r, edges, o_handle, kmer_len):
    header = cell_info
    num = 1
    all_reads_list = []
    for cntg in cell_contigs.keys():
        max_score, min_score, avg_score, read_cnt, reads_list = get_contig_support(
            cntg, edges_r, edges, k=kmer_len
        )
        o_handle.write(
            ">"
            + header
            + "__CTG"
            + str(num)
            + "__"
            + str(max_score)
            + "__"
            + str(min_score)
            + "__"
            + str(avg_score)
            + "__"
            + str(read_cnt)
            + "\n"
        )
        o_handle.write(cntg + "\n")
        num += 1
        all_reads_list.extend(reads_list)

    return len(list(set(all_reads_list)))


# Function for building contigs for each block. Returns all the contigs and edge information for the block.
def build_contig(
    v_edges,
    j_edges,
    bc_wind_list,
    bc_bci_list,
    bc_seq_list,
    seq_list,
    polyN_list,
    kmer_len=21,
):
    # Get edges information for the reads.
    edges = {}
    edges_r = {}
    rcnt = 0

    reads_overlapping_vj = 0
    for seq in seq_list:
        rcnt += 1
        edges, edges_r, read_only_edges = get_edgeinfo_from_sequence(
            edges, edges_r, rcnt, seq, k=kmer_len
        )
        for read_edge in read_only_edges:
            if read_edge in v_edges:
                reads_overlapping_vj += 1
                break
            elif read_edge in j_edges:
                reads_overlapping_vj += 1
                break

    # init
    all_contigs = {}
    track_all_edges = {}  # Tracker so we dont repeat

    # Get the most common edges and go down the list.
    edge_counter = Counter(edges).most_common()
    for edge_items in edge_counter:
        edge = edge_items[0]
        edge_count = edge_items[1]
        if edge in v_edges:  # Check if edge is in known v segment
            if edge in track_all_edges:  # Tracker so we dont repeat
                continue

            # Forward extension
            contig = [[]]
            done = {}
            get_fw_contig_greedy(contig, edge, edges, 0, done)
            track_all_edges = dict(track_all_edges, **done)
            # Re initialize for backward contig extension
            done = {}
            contig_bw = [[]]
            get_bw_contig_greedy(contig_bw, edge, edges, 0, done)
            track_all_edges = dict(track_all_edges, **done)

            for bw in contig_bw:
                bw.reverse()  # flip the reverse extended contigs

            for fw in contig:
                for bw in contig_bw:
                    cntg = "".join(bw[:-1]) + edge.split(":")[0] + "".join(fw)
                    if cntg not in all_contigs:
                        all_contigs[cntg] = 1
                    else:
                        all_contigs[cntg] += 1
        else:
            continue
    return all_contigs, edges, edges_r, reads_overlapping_vj


# Wrapper function for building contigs.
def get_contigs(
    polyN_file, v_edges, j_edges, polyN_contig_flag, read_threshold, kmer_len, trim=0
):
    fqfile = open(polyN_file, "r")
    bc_wind_list = []
    bc_ptype_list = []
    bc_bci_list = []
    bc_seq_list = []
    polyN_list = []
    seq_list = []
    qual_list = []
    cell = 0
    next_cell = False
    total_reads_used = {}
    o_fasta = polyN_file.split(".fastq")[0] + ".contigs.fa"
    o_handle = open(o_fasta, "w")  # output contig fasta file.
    total_reads_overlapping_vj = {}
    while True:
        aline = fqfile.readline()
        # If line is empty then end of file reached
        if not aline:
            break
        aline = aline.rstrip()
        read_name = aline[1:]
        read_sq = fqfile.readline().rstrip()
        if trim > 0:
            # trim first trim bases of each sequence
            read_sq = read_sq[trim:]
        tmp = fqfile.readline().rstrip()
        read_qual = fqfile.readline().rstrip()

        # Get info
        bc_wind = fqrec_to_bc_winds(read_name)
        bc_ptype = fqrec_to_ptype(read_name)
        bc_bcis = fqrec_to_bc_bcis(read_name)
        bc_seqs = fqrec_to_bc_seqs(read_name)
        polyN_seq = fqrec_to_polyN(read_name)

        # Append
        bc_wind_list.append(bc_wind)
        bc_ptype_list.append(bc_ptype)
        bc_bci_list.append(bc_bcis)
        bc_seq_list.append(bc_seqs)
        polyN_list.append(polyN_seq)
        seq_list.append(read_sq)
        qual_list.append(read_qual)

        cell += 1

        if polyN_contig_flag == 0:  # only build contigs for unique well barcode
            check_key = bc_wind_list
        else:  # build contig for unique polyN within each well barcode.
            check_key = [i + "\t" + j for i, j in zip(bc_wind_list, polyN_list)]
        if cell > 1:
            # Compare well indexes
            #  (Not bc seqs, which are uncorrected so may differ by chance)
            #  (Not bci barcode indexes, which may map to same well==cell )
            if check_key[-1] != check_key[-2]:
                next_cell = True
                # Process reads from one cell
            if next_cell:
                # NOTE: filtering based on number of reads.
                if len(check_key) - 1 >= read_threshold:
                    cell_contigs, edges, edges_r, reads_overlapping_vj = build_contig(
                        v_edges,
                        j_edges,
                        bc_wind_list[:-1],
                        bc_bci_list[:-1],
                        bc_seq_list[:-1],
                        seq_list[:-1],
                        polyN_list[:-1],
                        kmer_len=kmer_len,
                    )
                    # Reset lists, keeping the first read from the new cell
                    if polyN_contig_flag == 0:  # NOTE : dont include polyN in cell info
                        cell_info = "__".join(
                            [
                                bc_wind_list[-2],
                                bc_ptype_list[-2],
                                bc_bci_list[-2],
                                bc_seq_list[-2],
                                "NA",
                            ]
                        )
                    else:
                        cell_info = "__".join(
                            [
                                bc_wind_list[-2],
                                bc_ptype_list[-2],
                                bc_bci_list[-2],
                                bc_seq_list[-2],
                                polyN_list[-2],
                            ]
                        )
                    reads_used_count = get_support_wrapper(
                        cell_contigs,
                        cell_info,
                        edges_r,
                        edges,
                        o_handle,
                        kmer_len=kmer_len,
                    )

                    if bc_wind_list[-2] in total_reads_used:
                        total_reads_used[bc_wind_list[-2]] += reads_used_count
                        total_reads_overlapping_vj[
                            bc_wind_list[-2]
                        ] += reads_overlapping_vj
                    else:
                        total_reads_used[bc_wind_list[-2]] = reads_used_count
                        total_reads_overlapping_vj[
                            bc_wind_list[-2]
                        ] = reads_overlapping_vj

                    tr = len(bc_wind_list[:-1])
                    bc_wind_list = [bc_wind]
                    bc_ptype_list = [bc_ptype]
                    bc_bci_list = [bc_bcis]
                    bc_seq_list = [bc_seqs]
                    polyN_list = [polyN_seq]
                    seq_list = [read_sq]
                    qual_list = [read_qual]
                    cell = 1
                    next_cell = False
                else:
                    bc_wind_list = [bc_wind]
                    bc_ptype_list = [bc_ptype]
                    bc_bci_list = [bc_bcis]
                    bc_seq_list = [bc_seqs]
                    polyN_list = [polyN_seq]
                    seq_list = [read_sq]
                    qual_list = [read_qual]
                    cell = 1
                    next_cell = False

    if cell >= 1:
        if polyN_contig_flag == 0:  # NOTE : only build contigs for unique well barcode
            check_key = bc_wind_list
        else:  # NOTE: build contig for unique polyN within each well barcode.
            check_key = [i + "\t" + j for i, j in zip(bc_wind_list, polyN_list)]
        if len(check_key) - 1 >= read_threshold:
            cell_contigs, edges, edges_r, reads_overlapping_vj = build_contig(
                v_edges,
                j_edges,
                bc_wind_list,
                bc_bci_list,
                bc_seq_list,
                seq_list,
                polyN_list,
                kmer_len=kmer_len,
            )
            if polyN_contig_flag == 0:  # NOTE : dont include polyN in cell info
                cell_info = "__".join(
                    [
                        bc_wind_list[-1],
                        bc_ptype_list[-1],
                        bc_bci_list[-1],
                        bc_seq_list[-1],
                        "NA",
                    ]
                )
            else:
                cell_info = "__".join(
                    [
                        bc_wind_list[-1],
                        bc_ptype_list[-1],
                        bc_bci_list[-1],
                        bc_seq_list[-1],
                        polyN_list[-1],
                    ]
                )
            reads_used_count = get_support_wrapper(
                cell_contigs, cell_info, edges_r, edges, o_handle, kmer_len=kmer_len
            )

            if bc_wind_list[-2] in total_reads_used:
                total_reads_used[bc_wind_list[-1]] += reads_used_count
                total_reads_overlapping_vj[bc_wind_list[-1]] += reads_overlapping_vj
            else:
                total_reads_used[bc_wind_list[-1]] = reads_used_count
                total_reads_overlapping_vj[bc_wind_list[-1]] = reads_overlapping_vj

            tr = len(bc_wind_list)
    fqfile.close()
    o_handle.close()
    return o_fasta, total_reads_used, total_reads_overlapping_vj


###############################################################
# # Igblast functions
# Function to run igblast.
# Needs cleanup for the temp fix.
def run_igblast(igblast_dir, genome, query_file, out_file):
    org_name = genome.split("_")[-2]
    igblast_cmd = """{0}/bin/igblastn \
    -germline_db_V {0}/tr_database/{1}_V \
    -num_alignments_V 1 \
    -germline_db_D {0}/tr_database/{1}_D \
    -num_alignments_D 1 \
    -germline_db_J {0}/tr_database/{1}_J \
    -num_alignments_J 1 \
    -c_region_db {0}/tr_database/{1}_C \
    -num_alignments_C 1 \
    -auxiliary_data {0}/optional_file/{1}_gl.aux \
    -num_threads 16 \
    -outfmt 19 \
    -domain_system imgt -ig_seqtype TCR -organism {2} \
    -query {3} \
    -out {4}""".format(
        igblast_dir, genome, org_name, query_file, out_file
    )

    s = subprocess.getstatusoutput(igblast_cmd)


# Parse igblast output and add relevant columns
def parse_igblast_out(igblastfile, polyN_contig_flag=1):
    cols_list = [
        "sequence_id",
        "sequence",
        "locus",
        "rev_comp",
        "productive",
        "v_call",
        "d_call",
        "j_call",
        "c_call",
        "sequence_alignment",
        "germline_alignment",
        "cdr1",
        "cdr2",
        "junction",
        "junction_aa",
        "v_cigar",
        "d_cigar",
        "j_cigar",
        "v_identity",
        "j_identity",
        "complete_vdj",
    ]
    annotated_tcr = pd.read_csv(igblastfile, sep="\t", usecols=cols_list, dtype=object)
    # Get cell barcodes and n10 from header
    annotated_tcr = annotated_tcr.assign(
        bc_wind=annotated_tcr["sequence_id"].str.split("__").str[0]
    )
    annotated_tcr = annotated_tcr.assign(
        n10=annotated_tcr["sequence_id"].str.split("__").str[4]
    )
    annotated_tcr = annotated_tcr.assign(
        cell_barcode_n10=annotated_tcr["bc_wind"] + "__" + annotated_tcr["n10"]
    )
    # Only keep productive TCRs:
    annotated_tcr = annotated_tcr.query("productive=='T'")

    # Additional filter to make sure the CDR3 aa starts with a C and is atleast 5 AA long. we use junction aa here.
    annotated_tcr = annotated_tcr[
        annotated_tcr["junction_aa"].str.contains("^C[A-Z]{4,}", na=False)
    ]

    # Label the max_edge_count from the header
    annotated_tcr = annotated_tcr.assign(
        contignum=annotated_tcr["sequence_id"].str.split("__").str[5]
    )
    annotated_tcr = annotated_tcr.assign(
        max_edge_count=annotated_tcr["sequence_id"].str.split("__").str[6]
    )
    annotated_tcr = annotated_tcr.assign(
        min_edge_count=annotated_tcr["sequence_id"].str.split("__").str[7]
    )
    annotated_tcr = annotated_tcr.assign(
        avg_edge_count=annotated_tcr["sequence_id"].str.split("__").str[8]
    )
    annotated_tcr = annotated_tcr.assign(
        read_support=annotated_tcr["sequence_id"].str.split("__").str[9]
    )

    annotated_tcr.max_edge_count = annotated_tcr.max_edge_count.astype(float)
    annotated_tcr.min_edge_count = annotated_tcr.min_edge_count.astype(float)
    annotated_tcr.avg_edge_count = annotated_tcr.avg_edge_count.astype(float)
    annotated_tcr.read_support = annotated_tcr.read_support.astype(float)
    # Set a minimum number of reads per contig proportional to the max number of edge counts (arbitrarily set to 1/10)
    annotated_tcr = annotated_tcr.rename(columns={"sequence_id": "sequence_id_long"})
    if annotated_tcr.shape[0] == 0:
        annotated_tcr = annotated_tcr.assign(sequence_id=[])
    elif polyN_contig_flag == 1:
        annotated_tcr["sequence_id"] = annotated_tcr[
            ["bc_wind", "n10", "contignum"]
        ].agg("__".join, axis=1)
    else:
        annotated_tcr["sequence_id"] = annotated_tcr[["bc_wind", "contignum"]].agg(
            "__".join, axis=1
        )
    return annotated_tcr


# Collapsing igblast output helper functions.
def get_top_cdr(igblast_df, polyN_contig_flag=1, mode="n10"):
    if igblast_df.shape[0] == 0:
        return igblast_df
    else:
        if mode == "n10":
            # for n10/transcripts, choose contig based on number of times a cdr3aa is there and max_edge_count
            top_cdr3 = (
                igblast_df.groupby(["junction_aa"])["max_edge_count"]
                .aggregate(["size", "max"])
                .sort_values(by=["size", "max"], ascending=False)
                .index.values.tolist()[0]
            )
            igblast_df = igblast_df[igblast_df["junction_aa"] == top_cdr3]
            igblast_df = igblast_df.reset_index()
            chosen_one = igblast_df.iloc[[0]]
            return chosen_one
        if mode == "cell":
            # for cells, choose contigs based on transcript support and read support. Ignores transcript support if contigs built at cell level and not transcript level.
            n10_summ = (
                igblast_df.groupby(["junction_aa"])["read_support"]
                .aggregate(["size", "sum"])
                .sort_values(by=["size", "sum"], ascending=False)
            )
            return_df = pd.DataFrame()
            for i in range(n10_summ.shape[0]):
                igblast_df_tmp = igblast_df[
                    igblast_df["junction_aa"] == n10_summ.index.values.tolist()[i]
                ]
                igblast_df_tmp = igblast_df_tmp.sort_values(
                    by=["read_support", "max_edge_count", "avg_edge_count"],
                    ascending=False,
                )
                igblast_df_tmp = igblast_df_tmp.reset_index()
                chosen_one = igblast_df_tmp.iloc[[0]]
                if polyN_contig_flag == 1:
                    chosen_one = chosen_one.assign(
                        read_support=list(
                            n10_summ.loc[chosen_one["junction_aa"]]["sum"]
                        )
                    )
                    chosen_one = chosen_one.assign(
                        n10_supp=list(n10_summ.loc[chosen_one["junction_aa"]]["size"])
                    )
                return_df = pd.concat([return_df, chosen_one])
            return return_df


# wrapper
def collapse_wrapper(igblast_df, polyN_contig_flag=1, mode="n10"):
    return_df = pd.DataFrame()
    trb_chosen = get_top_cdr(
        igblast_df[igblast_df["locus"] == "TRB"], polyN_contig_flag, mode=mode
    )
    return_df = pd.concat([return_df, trb_chosen])
    tra_chosen = get_top_cdr(
        igblast_df[igblast_df["locus"] == "TRA"], polyN_contig_flag, mode=mode
    )
    return_df = pd.concat([return_df, tra_chosen])
    return return_df


# Function to collapse and choose representative contigs at a polyN/n10/transcript level
def collapse_cdr3_for_polyN(igblast_comb, read_support_threshold=1):
    # Remove rows with NA for junction_aa
    igblast_comb = igblast_comb[igblast_comb["junction_aa"].notna()]
    all_ids = list(set(igblast_comb["bc_wind"]))  # get all barcodes
    cell_igblast_condensed = pd.DataFrame(columns=igblast_comb.columns)
    for ids in all_ids:
        cell_igblast = igblast_comb[igblast_comb["bc_wind"] == ids]
        all_n10 = list(set(cell_igblast["n10"]))
        # condense contigs from each n10
        for polyn in all_n10:
            cell_n10_igblast = cell_igblast[cell_igblast["n10"] == polyn]
            cell_n10_igblast = cell_n10_igblast[
                cell_n10_igblast["read_support"] >= read_support_threshold
            ]  # filter out any contigs with read support less than threshold
            cell_n10_igblast = cell_n10_igblast.sort_values(
                by=["max_edge_count", "read_support", "avg_edge_count"], ascending=False
            )
            chosen_ones = collapse_wrapper(cell_n10_igblast, mode="n10")
            cell_igblast_condensed = pd.concat([cell_igblast_condensed, chosen_ones])
    return cell_igblast_condensed


# Function to collapse and choose representative contigs at a cell level
def collapse_cdr3_for_cell(n10_df, polyN_contig_flag=1, read_support_threshold=2):
    # Filter contigs based on read support
    median_read_support = n10_df["read_support"].median()
    read_cutoff = max(
        median_read_support / 10, read_support_threshold
    )  # NOTE: Fine tune this in future.
    n10_df_filtered = n10_df[n10_df["read_support"] >= read_cutoff]
    all_ids = list(set(n10_df_filtered["bc_wind"]))
    all_igblast = pd.DataFrame(columns=n10_df.columns)
    for ids in all_ids:
        cell_igblast = n10_df_filtered[n10_df_filtered["bc_wind"] == ids]
        chosen_ones = collapse_wrapper(
            cell_igblast, polyN_contig_flag=polyN_contig_flag, mode="cell"
        )
        all_igblast = pd.concat([all_igblast, chosen_ones])
    return all_igblast


# Combine igblast files
# Just concatenates input igblast files. Handles both condensed and not condensed files.
def igblast_combine(wkdir):
    all_igblast_files = glob.glob(wkdir + "/*chunk*.igblast.productive.tsv")
    cols_list = [
        "sequence_id",
        "sequence",
        "locus",
        "rev_comp",
        "productive",
        "v_call",
        "d_call",
        "j_call",
        "c_call",
        "sequence_alignment",
        "germline_alignment",
        "cdr1",
        "cdr2",
        "junction",
        "junction_aa",
        "v_cigar",
        "d_cigar",
        "j_cigar",
        "v_identity",
        "j_identity",
        "complete_vdj",
    ]
    dtype_list = {}
    for i in cols_list:
        dtype_list[i] = object
    igblast_comb = pd.DataFrame()
    if len(all_igblast_files) == 0:
        print("No igblast files")
    else:
        for ifile in all_igblast_files:
            annotated_tcr = pd.read_csv(ifile, sep="\t", dtype=dtype_list)
            igblast_comb = pd.concat([igblast_comb, annotated_tcr])

    all_igblast_files = glob.glob(wkdir + "/*chunk*.igblast.condensed.tsv")
    igblast_condensed = pd.DataFrame()
    if len(all_igblast_files) == 0:
        print("No igblast files")
    else:
        for ifile in all_igblast_files:
            annotated_tcr = pd.read_csv(ifile, sep="\t", dtype=dtype_list)
            igblast_condensed = pd.concat([igblast_condensed, annotated_tcr])
    return igblast_comb, igblast_condensed


###############################################################
# # Extra functions
# Get edges information for V sequences.
def get_vj_edges(vj_seqment_fasta, kmer_len):
    vj_seqment_seqs = []
    with open(vj_seqment_fasta) as inFile:
        seqblock = []
        for aline in inFile.readlines():
            aline = aline.rstrip()
            if aline.startswith(">"):
                if seqblock:
                    concat_seq = "".join(seqblock)
                    concat_seq2 = (
                        concat_seq.upper()
                    )  # NOTE: if working with igblast processed fasta, the "." are already removed.
                    seqblock = []
                    vj_seqment_seqs.append(concat_seq2)
            else:
                seqblock.append(aline)

        if seqblock:
            concat_seq = "".join(seqblock)
            concat_seq2 = concat_seq.replace(".", "").upper()
            vj_seqment_seqs.append(concat_seq2)

    # get V kmers and edges
    kmers = {}
    edges = {}
    for seq in vj_seqment_seqs:
        kmers, edges = get_kmer_and_edges_from_sequence(kmers, edges, seq, k=kmer_len)

    return kmers, edges


###############################################################
# Overall wrapper.
def run_tcr(
    cnter,
    odir,
    chunk_file,
    genome,
    v_seqment_fasta,
    j_segment_fasta,
    polyN_contig_flag,
    read_threshold,
    igblast_dir,
    kmer_len,
    trim=0,
):
    # Run polyN correction
    polyN_file = run_polyN(chunk_file)

    # Get V gene edges
    v_kmers, v_edges = get_vj_edges(v_seqment_fasta, kmer_len)
    j_kmers, j_edges = get_vj_edges(j_segment_fasta, kmer_len)

    # Build contigs
    contig_fasta, total_reads_used, total_reads_overlapping_vj = get_contigs(
        polyN_file=polyN_file,
        v_edges=v_edges,
        j_edges=j_edges,
        polyN_contig_flag=polyN_contig_flag,
        read_threshold=read_threshold,
        kmer_len=kmer_len,
        trim=trim,
    )
    # store tmp files with some metrics
    met_file = odir + "/reads_used_" + str(cnter) + "_thread.txt"
    with open(met_file, "w") as ofile:
        ofile.write("barcode\treads_used\treads_vj\n")
        for key, value in total_reads_used.items():
            ofile.write(
                key
                + "\t"
                + str(value)
                + "\t"
                + str(total_reads_overlapping_vj[key])
                + "\n"
            )
    # Igblast
    o_tsv = contig_fasta.split(".fa")[0] + ".igblast.tsv"
    c = run_igblast(
        igblast_dir=igblast_dir, genome=genome, query_file=contig_fasta, out_file=o_tsv
    )
    # parse igblast outputs
    annotated_tcr = parse_igblast_out(o_tsv, polyN_contig_flag)
    annofile = contig_fasta.split(".fa")[0] + ".igblast.productive.tsv"
    # write parsed output to file
    annotated_tcr.to_csv(annofile, sep="\t")
    # condense igblast outputs to choose best contigs.
    if polyN_contig_flag == 1:
        n10_collapse = collapse_cdr3_for_polyN(
            annotated_tcr, read_support_threshold=read_threshold
        )
        igblast_condensed = collapse_cdr3_for_cell(
            n10_collapse, read_support_threshold=read_threshold
        )
    else:
        igblast_condensed = collapse_cdr3_for_cell(
            annotated_tcr,
            polyN_contig_flag=polyN_contig_flag,
            read_support_threshold=read_threshold,
        )
    # write condensed file out
    igfile = contig_fasta.split(".fa")[0] + ".igblast.condensed.tsv"
    igblast_condensed.to_csv(igfile, sep="\t")


###############################################################
# # Summary file generation functions
# Generates the AIRR file output.
def make_airr_file(igblast_in, ofile, polyN_contig_flag=1):
    if polyN_contig_flag == 1:
        cols_list = [
            "sequence_id",
            "sequence",
            "locus",
            "rev_comp",
            "productive",
            "v_call",
            "d_call",
            "j_call",
            "c_call",
            "sequence_alignment",
            "germline_alignment",
            "cdr1",
            "cdr2",
            "junction",
            "junction_aa",
            "v_cigar",
            "d_cigar",
            "j_cigar",
            "v_identity",
            "j_identity",
            "bc_wind",
            "complete_vdj",
            "read_support",
            "n10_supp",
        ]
    else:
        cols_list = [
            "sequence_id",
            "sequence",
            "locus",
            "rev_comp",
            "productive",
            "v_call",
            "d_call",
            "j_call",
            "c_call",
            "sequence_alignment",
            "germline_alignment",
            "cdr1",
            "cdr2",
            "junction",
            "junction_aa",
            "v_cigar",
            "d_cigar",
            "j_cigar",
            "v_identity",
            "j_identity",
            "bc_wind",
            "complete_vdj",
            "read_support",
        ]
    igblast_sub = igblast_in[cols_list]

    igblast_sub = igblast_sub.rename(
        columns={
            "junction": "cdr3",
            "junction_aa": "cdr3_aa",
            "read_support": "read_count",
            "bc_wind": "cell_barcode",
        }
    )
    if polyN_contig_flag == 1:
        igblast_sub = igblast_sub.rename(columns={"n10_supp": "transcript_count"})
    igblast_sub.to_csv(ofile, index=False, sep="\t")


# Function to check if cells are possible doublets. At this level, the check is done to see if the cell has
# more than 2 TRB or TRA and will be labelled as doublet if so.
def doublet_check(indf, nth, odir, polyN_contig_flag=1):
    all_ids = list(set(indf["bc_wind"]))
    clon_df = pd.DataFrame(
        columns=[
            "barcode",
            "TRA_primary",
            "TRB_primary",
            "TRA_secondary",
            "TRB_secondary",
            "isDoublet",
        ]
    )
    for ids in all_ids:
        idf_sub = indf.query("bc_wind == @ids")
        if polyN_contig_flag == 1:
            idf_sub = idf_sub.sort_values(by=["n10_supp"], ascending=False)
            if "level_0" in idf_sub.columns:
                idf_sub = idf_sub.drop(["level_0"], axis=1)
        else:
            idf_sub = idf_sub.sort_values(by=["read_support"], ascending=False)
        trb_all = idf_sub[idf_sub["locus"] == "TRB"].reset_index()
        tra_all = idf_sub[idf_sub["locus"] == "TRA"].reset_index()
        isdub = 0
        if trb_all.shape[0] > 2 or tra_all.shape[0] > 2:
            isdub = 1
        trb1 = "NA"
        trb2 = "NA"
        tra1 = "NA"
        tra2 = "NA"
        if trb_all.shape[0] > 0:
            if trb_all.shape[0] == 1:
                trb1 = trb_all.iloc[0][["junction_aa"]].values.flatten().tolist()[0]
            if trb_all.shape[0] > 1:
                trb1 = trb_all.iloc[0][["junction_aa"]].values.flatten().tolist()[0]
                trb2 = trb_all.iloc[1][["junction_aa"]].values.flatten().tolist()[0]
        if tra_all.shape[0] > 0:
            if tra_all.shape[0] == 1:
                tra1 = tra_all.iloc[0][["junction_aa"]].values.flatten().tolist()[0]
            if tra_all.shape[0] > 1:
                tra1 = tra_all.iloc[0][["junction_aa"]].values.flatten().tolist()[0]
                tra2 = tra_all.iloc[1][["junction_aa"]].values.flatten().tolist()[0]
        id_dict = {
            "barcode": [ids],
            "TRA_primary": [tra1],
            "TRB_primary": [trb1],
            "TRA_secondary": [tra2],
            "TRB_secondary": [trb2],
            "isDoublet": [isdub],
        }
        id_df = pd.DataFrame(data=id_dict)
        clon_df = pd.concat([clon_df, id_df])
    tmpofile = odir + "/clon_df_" + str(nth) + ".tsv"

    clon_df.to_csv(tmpofile, sep="\t", index=False)


# Helper function for getting clonotype information. It chooses the top clonotype for each TRA and TRB.
def condense_clonotype(indf):
    tra_list = indf.TRA_primary.tolist()
    trb_list = indf.TRB_primary.tolist()
    out_tra = []
    out_trb = []
    opair = {}
    out_size = []
    for i in range(len(tra_list)):
        tra = tra_list[i]
        trb = trb_list[i]
        if tra != "NA":
            if trb != "NA":
                p = tra + "-" + trb
                if p not in opair:
                    out_tra.append(tra)
                    out_trb.append(trb)
                    out_size.append(0)
                    opair[p] = 1

    odict = {"TRA_primary": out_tra, "TRB_primary": out_trb, "size": out_size}
    odf = pd.DataFrame(data=odict)
    return odf


# CURRENTLY NOT USED : Generates clonotype frequency file and dataframe with same information.
# Logic :
# 1. All primary clonotypes are sorted by the number of occurences. Those with "NA" for TRA or TRB are pushed to the bottom.
# 2. It is collapsed such that TRA and TRBs can only be there once since the same TRA or TRB cannot have multiple pairs.
#      This list is used as reference list of clonotypes.
# 3. The primary TRA and TRB information for each cell is used to assign it to a clonotype as follows.
#    a. If the cell has the exact same pair as a reference, it is considered a "match".
#         NOTE: matches can occur when the cell has "NA" for one of TRA/TRB.
#         This will happen when the reference also has a singleton clonotype.
#    b. If the cell has only a TRA or TRB and it matches a reference clonotype but is missing the partner TR, it is marked as "inferred"
#    c. Similarly, if a cell has a TRA or TRB which is NOT in the reference (happens when its not top match with any partner).
#         The missing one is treated as "NA" and check for "inferred".
#    d. If a cell has TRA that matches a reference and a TRB that matches a different reference, it is marked as a doublet.
def old_count_productive_metrics(tra, trb, clon_conf):
    prod_a = 0
    prod_b = 0
    prod_pair = 0
    inf_prod_pair = 0
    if clon_conf == "match":
        if (tra != "NA") & (trb != "NA"):
            prod_pair = 1
            prod_a = 1
            prod_b = 1
        elif (tra != "NA") & (trb == "NA"):
            prod_a = 1
        elif (tra == "NA") & (trb != "NA"):
            prod_b = 1
    # elif clon_conf == "inferred":
    #    if (tra != "NA") & (trb != "NA"):
    #        #print("Error : Possbile doublet being passed as inferred")
    #        inf_prod_pair =1
    #        prod_a = 1
    #        prod_b = 1
    #    elif (tra != "NA") & (trb == "NA"):
    #        prod_a = 1
    #        inf_prod_pair = 1
    #    elif (tra == "NA") & (trb != "NA"):
    #        prod_b = 1
    #        inf_prod_pair = 1
    elif clon_conf == "NA":
        if (tra != "NA") & (trb != "NA"):
            prod_a = 1
            prod_b = 1
            prod_pair = 1
        elif (tra != "NA") & (trb == "NA"):
            prod_a = 1
        elif (tra == "NA") & (trb != "NA"):
            prod_b = 1
    return prod_pair, inf_prod_pair, prod_a, prod_b


def count_productive_metrics(tra, trb):
    prod_a = 0
    prod_b = 0
    prod_pair = 0
    inf_prod_pair = 0
    if (tra != "NA") & (trb != "NA"):
        prod_pair = 1
        prod_a = 1
        prod_b = 1
    elif (tra != "NA") & (trb == "NA"):
        prod_a = 1
    elif (tra == "NA") & (trb != "NA"):
        prod_b = 1
    return prod_pair, inf_prod_pair, prod_a, prod_b


def clonotype_summary(clon_in, freq_track, nth, odir):
    prod_dub_info = {}
    clon_info = pd.DataFrame(
        columns=["barcode", "isDoublet", "clonotype_id", "clonotype_confidence"]
    )
    get_primary_rep = freq_track.copy(deep=True)
    cnter = 0
    time_s = time.time()
    # Track metrics
    for ids in clon_in.barcode.values.tolist():
        cnter += 1
        prod_pair = 0
        inf_prod_pair = 0
        prod_a = 0
        prod_b = 0
        tmp = clon_in.query("barcode == @ids")
        if tmp.shape[0] > 1:
            tmp = tmp.iloc[0].to_frame().transpose()
        clon_id = "NA"
        clon_conf = "NA"
        tra = tmp["TRA_primary"].tolist()[0]
        trb = tmp["TRB_primary"].tolist()[0]
        isdub = tmp["isDoublet"].tolist()[0]
        prod_pair, inf_prod_pair, prod_a, prod_b = count_productive_metrics(tra, trb)
        checkdf = (get_primary_rep["TRA_primary"].isin([tra])) & (
            get_primary_rep["TRB_primary"].isin([trb])
        )
        if isdub == 1:
            clon_conf = "NA"
        elif sum(checkdf) == 1:
            clon_conf = "match"
            get_primary_rep.loc[checkdf, "size"] = get_primary_rep[checkdf]["size"] + 1
            clon_id = get_primary_rep.loc[checkdf, "clonotype_id"].tolist()[0]

        elif sum(checkdf) > 1:
            print("Too many top clonotype matches")
        id_dict = {
            "barcode": [ids],
            "isDoublet": [isdub],
            "clonotype_id": [clon_id],
            "clonotype_confidence": [clon_conf],
        }
        id_df = pd.DataFrame(data=id_dict)
        clon_info = pd.concat([clon_info, id_df])

        if ids in prod_dub_info:
            print("Duplicate " + ids + " in clonotype_summary")
        else:
            prod_dub_info[ids] = [prod_a, prod_b, prod_pair, inf_prod_pair, isdub]
    get_primary_rep.sort_values(by=["size"], ascending=[False])
    get_primary_rep = get_primary_rep.rename(columns={"size": "count"})
    get_primary_rep = get_primary_rep[
        ["TRA_primary", "TRB_primary", "clonotype_id", "count"]
    ]

    ofreq = odir + "/clon_freq_" + str(nth) + ".tsv"
    ocloninfo = odir + "/clon_info_" + str(nth) + ".tsv"
    oprod_dub_info = odir + "/prod_dub_" + str(nth) + ".tsv"
    get_primary_rep.to_csv(ofreq, index=False, sep="\t")
    clon_info.to_csv(ocloninfo, index=False, sep="\t")
    with open(oprod_dub_info, "w") as outFile:
        for key, value in prod_dub_info.items():
            outFile.write(key + "\t" + "\t".join(map(str, value)) + "\n")


# Function to output the report.tsv file.
def get_letter_well(index):
    wells = np.arange(96, dtype=int).reshape(8, 12)
    wells = wells + 1
    row_number_to_letter = dict(zip([i for i in range(8)], list("ABCDEFGH")))
    loc1 = np.where(wells == index)[0][0]
    loc2 = np.where(wells == index)[1][0]
    return row_number_to_letter[loc1] + str(loc2 + 1)


def make_report_file(igblast_in, clon_df, odir, i, polyN_contig_flag=1):
    ofile = odir + "/barcode_report_" + str(i) + ".tsv"
    ffile = odir + "/contigs_fa_" + str(i) + ".fa"
    # Get all bc_wind ids.
    all_ids = list(set(igblast_in["bc_wind"]))
    # outputfile
    out_file = open(ofile, "w")
    fa_file = open(ffile, "w")
    header = []
    if polyN_contig_flag == 1:
        header = [
            "Barcode",
            "TRA_V",
            "TRA_D",
            "TRA_J",
            "TRA_C",
            "TRA_cdr3_aa",
            "TRA_read_count",
            "TRA_transcript_count",
            "TRB_V",
            "TRB_D",
            "TRB_J",
            "TRB_C",
            "TRB_cdr3_aa",
            "TRB_read_count",
            "TRB_transcript_count",
            "secondary_TRA_V",
            "secondary_TRA_D",
            "secondary_TRA_J",
            "secondary_TRA_C",
            "secondary_TRA_cdr3_aa",
            "secondary_TRA_read_count",
            "secondary_TRA_transcript_count",
            "secondary_TRB_V",
            "secondary_TRB_D",
            "secondary_TRB_J",
            "secondary_TRB_C",
            "secondary_TRB_cdr3_aa",
            "secondary_TRB_read_count",
            "secondary_TRB_transcript_count",
            "isMultiplet",
            "bc1_well",
            "bc2_well",
            "bc3_well",
            "bc1_wind",
            "bc2_wind",
            "bc3_wind",
        ]
        #'assigned_clonotype_id',
        #'clonotype_assignment_type']
    else:
        header = [
            "Barcode",
            "TRA_V",
            "TRA_D",
            "TRA_J",
            "TRA_C",
            "TRA_cdr3_aa",
            "TRA_read_support",
            "TRB_V",
            "TRB_D",
            "TRB_J",
            "TRB_C",
            "TRB_cdr3_aa",
            "TRB_read_support",
            "secondary_TRA_V",
            "secondary_TRA_D",
            "secondary_TRA_J",
            "secondary_TRA_C",
            "secondary_TRA_cdr3_aa",
            "secondary_TRA_read_support",
            "secondary_TRB_V",
            "secondary_TRB_D",
            "secondary_TRB_J",
            "secondary_TRB_C",
            "secondary_TRB_cdr3_aa",
            "secondary_TRB_read_support",
            "isMultiplet",
            "bc1_well",
            "bc2_well",
            "bc3_well",
            "bc1_wind",
            "bc2_wind",
            "bc3_wind",
        ]
        #'assigned_clonotype_id',
        #'clonotype_assignment_type']
    out_file.write("\t".join(header) + "\n")
    cnter = 0
    time_s = time.time()
    for ids in all_ids:
        cnter += 1
        if (cnter % 10000) == 0:
            print(str(cnter) + " done")
            tt = time.time() - time_s
            print(tt)
            time_s = time.time()
        tmp = igblast_in.query("bc_wind == @ids")
        cols_list = []
        if polyN_contig_flag == 1:
            cols_list = [
                "v_call",
                "d_call",
                "j_call",
                "c_call",
                "junction_aa",
                "read_support",
                "n10_supp",
            ]
            if "level_0" in tmp.columns:
                tmp = tmp.drop(["level_0"], axis=1)
            tmp_top = tmp.sort_values(by=["n10_supp"], ascending=False)
        else:
            cols_list = [
                "v_call",
                "d_call",
                "j_call",
                "c_call",
                "junction_aa",
                "read_support",
            ]
            tmp_top = tmp.sort_values(by=["read_support"], ascending=False)
        trb1 = "\t".join(["NA"] * len(cols_list))
        tra1 = "\t".join(["NA"] * len(cols_list))
        trb2 = "\t".join(["NA"] * len(cols_list))
        tra2 = "\t".join(["NA"] * len(cols_list))
        seq_ids = []
        trb_df = tmp_top[tmp_top["locus"] == "TRB"].reset_index()
        tra_df = tmp_top[tmp_top["locus"] == "TRA"].reset_index()
        if trb_df.shape[0] > 0:
            if trb_df.shape[0] == 1:
                trb1 = "\t".join(
                    map(str, trb_df.iloc[0][cols_list].values.flatten().tolist())
                )
                seq_ids.append(trb_df.iloc[0]["sequence_id"])
            if trb_df.shape[0] > 1:
                trb1 = "\t".join(
                    map(str, trb_df.iloc[0][cols_list].values.flatten().tolist())
                )
                trb2 = "\t".join(
                    map(str, trb_df.iloc[1][cols_list].values.flatten().tolist())
                )
                seq_ids.append(trb_df.iloc[0]["sequence_id"])
                seq_ids.append(trb_df.iloc[1]["sequence_id"])
        if tra_df.shape[0] > 0:
            if tra_df.shape[0] == 1:
                tra1 = "\t".join(
                    map(str, tra_df.iloc[0][cols_list].values.flatten().tolist())
                )
                seq_ids.append(tra_df.iloc[0]["sequence_id"])
            if tra_df.shape[0] > 1:
                tra1 = "\t".join(
                    map(str, tra_df.iloc[0][cols_list].values.flatten().tolist())
                )
                tra2 = "\t".join(
                    map(str, tra_df.iloc[1][cols_list].values.flatten().tolist())
                )
                seq_ids.append(tra_df.iloc[0]["sequence_id"])
                seq_ids.append(tra_df.iloc[1]["sequence_id"])

        clon_df_bc = clon_df.query("barcode == @ids")
        if clon_df_bc.shape[0] > 1:
            print("too many clon entries")
        if len(clon_df_bc["isDoublet"].tolist()) == 0:
            print(clon_df_bc)
        isdub = clon_df_bc["isDoublet"].tolist()[0]
        clon_id = clon_df_bc["clonotype_id"].tolist()[0]
        clon_conf = clon_df_bc["clonotype_confidence"].tolist()[0]

        # barcode info
        bc_split = ids.split("_")
        bc_wells_info = []
        for i in [0, 1, 2]:
            well_info = get_letter_well(int(bc_split[i]))
            bc_wells_info.append(well_info)
        well_data = "\t".join(bc_wells_info)
        bc_data = "\t".join(bc_split[0:3])

        outline = "\t".join(
            [ids, tra1, trb1, tra2, trb2, str(isdub), well_data, bc_data]
        )
        out_file.write(outline + "\n")

        if isdub == 0:
            for ii in seq_ids:
                fa_file.write(ii + "\n")
    out_file.close()
    fa_file.close()


# Function to output fasta file with the contigs. Does not output doublets and only outputs primary chains.
def make_fasta_file(igblast_in, fasta_seq_ids, ofile, polyN_contig_flag=1):
    out_file = open(ofile, "w")
    cols_list = []
    if polyN_contig_flag == 1:
        cols_list = [
            "sequence_id",
            "locus",
            "v_call",
            "d_call",
            "j_call",
            "c_call",
            "junction_aa",
            "read_support",
            "n10_supp",
        ]
    else:
        cols_list = [
            "sequence_id",
            "locus",
            "v_call",
            "d_call",
            "j_call",
            "c_call",
            "junction_aa",
            "read_support",
        ]

    notdub = igblast_in[igblast_in["sequence_id"].isin(fasta_seq_ids)]
    info = notdub.assign(
        header=notdub[cols_list].apply(
            lambda x: ">" + "|".join(x.dropna().astype(str)), axis=1
        )
    )
    info[["header", "sequence"]].to_csv(ofile, sep="\n", header=False, index=False)


def make_summary_file(indir, total_reads, total_reads_per_cell, prod_dub_info):
    summary_df = pd.DataFrame.from_dict(
        total_reads_per_cell, orient="index", columns=["reads_per_cell"]
    )
    if total_reads != summary_df.reads_per_cell.sum():
        print("Total reads not matching sum of reads per cell in make_summary_file.")

    cols_list = [
        "sequence_id",
        "sequence",
        "locus",
        "productive",
        "junction",
        "junction_aa",
    ]
    filenames = glob.glob(indir + "/*.igblast.tsv")

    summary_df["isTCR"] = 0
    for infile in filenames:
        annotated_tcr = pd.read_csv(infile, sep="\t", usecols=cols_list, dtype=object)
        # Get cell barcodes and n10 from header
        annotated_tcr = annotated_tcr.assign(
            bc_wind=annotated_tcr["sequence_id"].str.split("__").str[0]
        )
        annotated_tcr = annotated_tcr.assign(
            n10=annotated_tcr["sequence_id"].str.split("__").str[4]
        )
        annotated_tcr = annotated_tcr.assign(
            cell_barcode_n10=annotated_tcr["bc_wind"] + "__" + annotated_tcr["n10"]
        )
        # Change made so that only cells with productive contigs will get labeled as "T cells"
        # Only keep productive TCRs:
        annotated_tcr = annotated_tcr.query("productive=='T'")

        # Additional filter to make sure the CDR3 aa starts with a C and is atleast 5 AA long. we use junction aa here.
        annotated_tcr = annotated_tcr[
            annotated_tcr["junction_aa"].str.contains("^C[A-Z]{4,}", na=False)
        ]

        # track cells with any TCR and VDJ transcripts per cell for metrics
        cells_with_TCR = pd.unique(
            annotated_tcr["bc_wind"].sort_values(ascending=False)
        ).tolist()
        summary_df.loc[cells_with_TCR, "isTCR"] = 1
        number_of_vdj_transcripts_per_cell = (
            annotated_tcr.groupby(["bc_wind", "n10"], as_index=False)
            .size()
            .groupby(["bc_wind"])
            .size()
        )
        summary_df.loc[
            number_of_vdj_transcripts_per_cell.index, "VDJ_transcripts_per_cell"
        ] = number_of_vdj_transcripts_per_cell
    summary_df["VDJ_transcripts_per_cell"] = summary_df[
        "VDJ_transcripts_per_cell"
    ].fillna(0)
    # Get reads that are used to build contigs and reads that overlap VJ for cells
    sub_reads_used = 0
    vj_reads = 0
    reads_used_files = glob.glob(indir + "/reads_used_*_thread.txt")
    reads_used_df = pd.DataFrame()
    for infile in reads_used_files:
        tmp_df = pd.read_csv(infile, sep="\t", index_col="barcode")
        reads_used_df = pd.concat([reads_used_df, tmp_df])
    summary_df = summary_df.join(reads_used_df)

    # get productive info and doublet info
    prod_dub_df = pd.DataFrame.from_dict(
        prod_dub_info,
        orient="index",
        columns=["prod_a", "prod_b", "prod_pair", "inf_prod_pair", "isdub"],
    )
    summary_df = summary_df.join(prod_dub_df)

    # fill na
    summary_df = summary_df.fillna(0)

    #
    out_sum_file = indir + "/../tcr_barcode_summary.txt"
    summary_df.to_csv(
        out_sum_file, sep="\t", header=True, index=True, index_label="barcode"
    )
    return summary_df


def generate_summary_files(indf, opfix, polyN_contig_flag):
    make_airr_file(
        indf,
        ofile=opfix + "tcr_annotation_airr.tsv",
        polyN_contig_flag=polyN_contig_flag,
    )
    clon_df = doublet_check(indf, polyN_contig_flag=polyN_contig_flag)
    get_primary_rep, clon_info, prod_dub_info = clonotype_summary(
        clon_df, ofile=opfix + "clonotype_frequency.tsv"
    )

    fasta_seq_ids = make_report_file(
        indf,
        clon_info,
        ofile=opfix + "barcode_report.tsv",
        polyN_contig_flag=polyN_contig_flag,
    )
    make_fasta_file(
        indf,
        fasta_seq_ids,
        ofile=opfix + "tcr_contigs.fa",
        polyN_contig_flag=polyN_contig_flag,
    )
    return prod_dub_info


def combine_clonotype_freq_from_threads(output_dir, nthreads):
    gpr_dict = {}
    for i in range(nthreads):
        tmpifile = output_dir + "/clon_freq_" + str(i) + ".tsv"
        with open(tmpifile, "r") as inFile:
            for aline in inFile.readlines():
                aline = aline.rstrip()
                if re.match("^TRA_primary.*", aline) is None:
                    vals = aline.split("\t")
                    key = vals[0] + "\t" + vals[1] + "\t" + vals[2]
                    if key in gpr_dict:
                        gpr_dict[key] += int(vals[3])
                    else:
                        gpr_dict[key] = int(vals[3])
    gpr_dict2 = {}
    gpr_dict2["TRA"] = []
    gpr_dict2["TRB"] = []
    gpr_dict2["clonotype_id"] = []
    gpr_dict2["count"] = []
    for key, value in gpr_dict.items():
        spl = key.split("\t")
        gpr_dict2["TRA"].append(spl[0])
        gpr_dict2["TRB"].append(spl[1])
        gpr_dict2["clonotype_id"].append(spl[2])
        gpr_dict2["count"].append(value)

    gpr_df = pd.DataFrame(data=gpr_dict2)
    gpr_df = gpr_df.assign(frequency=gpr_df["count"] / sum(gpr_df["count"]))
    gpr_df2 = gpr_df[gpr_df["count"] != 0]
    gpr_df2 = gpr_df2.rename(columns={"TRA_primary": "TRA", "TRB_primary": "TRB"})
    gpr_df2 = gpr_df2.sort_values(by=["count"], ascending=[False])
    return gpr_df2


def generate_summary_files_threaded(indf, opfix, polyN_contig_flag, nthreads, proc_dir):
    make_airr_file(
        indf,
        ofile=opfix + "tcr_annotation_airr.tsv",
        polyN_contig_flag=polyN_contig_flag,
    )
    # Create a tmp directory for threaded outputs and delete at end
    output_dir = proc_dir + "/tmp/"
    mk_cmd = "mkdir -p " + output_dir
    c = subprocess.getstatusoutput(mk_cmd)
    # Run doublet check in threaded manner
    ##run
    all_ids = list(set(indf["bc_wind"]))
    # Check if number is less than threads
    if len(all_ids) < nthreads:
        nthreads = 1
    per_thread = int(round(len(all_ids) / nthreads, 0))
    Pros = []
    for i in range(nthreads):
        if i == nthreads - 1:
            tmp_ids = all_ids[i * per_thread : len(all_ids)]
        else:
            tmp_ids = all_ids[i * per_thread : (i + 1) * per_thread]
        tmp_airr = indf[indf["bc_wind"].isin(tmp_ids)]
        p = Process(
            target=doublet_check, args=(tmp_airr, i, output_dir, polyN_contig_flag)
        )
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()
    ##combine outs
    clon_df = pd.DataFrame(
        columns=[
            "barcode",
            "TRA_primary",
            "TRB_primary",
            "TRA_secondary",
            "TRB_secondary",
            "isDoublet",
        ]
    )
    for i in range(nthreads):
        tmpifile = output_dir + "/clon_df_" + str(i) + ".tsv"
        tmpdf = pd.read_csv(tmpifile, sep="\t")
        clon_df = pd.concat([clon_df, tmpdf])
    clon_df = clon_df.fillna("NA")

    # Run clonotyping in threaded manner
    ##get files set up
    get_primary_rep = clon_df.groupby(
        ["TRA_primary", "TRB_primary"], as_index=False
    ).size()
    get_primary_rep["NA"] = get_primary_rep.applymap(lambda x: "NA" == str(x)).any(
        axis=1
    )  # So only full pairs are on top
    get_primary_rep = get_primary_rep.sort_values(
        by=["NA", "size"], ascending=[True, False]
    )
    get_primary_rep = condense_clonotype(get_primary_rep)
    get_primary_rep = get_primary_rep.assign(
        clonotype_id=["clonotype_" + str(i) for i in range(1, len(get_primary_rep) + 1)]
    )
    ##run
    per_thread = int(round(clon_df.shape[0] / nthreads, 0))
    Pros = []
    for i in range(nthreads):
        if i == nthreads - 1:
            tmp_clon_df = clon_df.iloc[i * per_thread : clon_df.shape[0]]
        else:
            tmp_clon_df = clon_df.iloc[i * per_thread : (i + 1) * per_thread]
        p = Process(
            target=clonotype_summary, args=(tmp_clon_df, get_primary_rep, i, output_dir)
        )
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()
    ##combine outs
    combined_clon_freq = combine_clonotype_freq_from_threads(output_dir, nthreads)
    ofile = opfix + "clonotype_frequency.tsv"
    combined_clon_freq.to_csv(ofile, index=False, sep="\t")
    ##combine clon_info
    clon_info = pd.DataFrame(
        pd.DataFrame(
            columns=["barcode", "isDoublet", "clonotype_id", "clonotype_confidence"]
        )
    )
    for i in range(nthreads):
        tmpifile = output_dir + "/clon_info_" + str(i) + ".tsv"
        tmpdf = pd.read_csv(tmpifile, sep="\t")
        clon_info = pd.concat([clon_info, tmpdf])
    clon_info = clon_info.fillna("NA")
    ##combine prod_dub
    prod_dub = {}
    for i in range(nthreads):
        tmpifile = output_dir + "/prod_dub_" + str(i) + ".tsv"
        with open(tmpifile, "r") as inFile:
            for aline in inFile.readlines():
                aline = aline.rstrip()
                vals = aline.split("\t")
                if vals[0] in prod_dub:
                    print("Duplicate " + vals[0])
                else:
                    prod_dub[vals[0]] = [
                        int(vals[1]),
                        int(vals[2]),
                        int(vals[3]),
                        int(vals[4]),
                        int(vals[5]),
                    ]

    # Run report file in threaded manner
    ##run
    all_ids = list(set(indf["bc_wind"]))
    per_thread = int(round(len(all_ids) / nthreads, 0))
    Pros = []
    for i in range(nthreads):
        if i == nthreads - 1:
            tmp_ids = all_ids[i * per_thread : len(all_ids)]
        else:
            tmp_ids = all_ids[i * per_thread : (i + 1) * per_thread]
        tmp_airr = indf[indf["bc_wind"].isin(tmp_ids)]
        p = Process(
            target=make_report_file,
            args=(tmp_airr, clon_info, output_dir, i, polyN_contig_flag),
        )
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()
    ##combine outs
    header = 0
    ofile = open(opfix + "barcode_report.tsv", "w")
    for i in range(nthreads):
        tmpifile = output_dir + "/barcode_report_" + str(i) + ".tsv"
        with open(tmpifile, "r") as ifile:
            for aline in ifile.readlines():
                aline = aline.rstrip()
                if (header == 0) and (re.match("^Barcode.*", aline) is not None):
                    header = 1
                    ofile.write(aline + "\n")
                elif re.match("^Barcode.*", aline) is None:
                    ofile.write(aline + "\n")
    fasta_seq_ids = []
    for i in range(nthreads):
        tmpifile = output_dir + "/contigs_fa_" + str(i) + ".fa"
        with open(tmpifile, "r") as ifile:
            for aline in ifile.readlines():
                aline = aline.rstrip()
                fasta_seq_ids.append(aline)
    # Make fasta
    make_fasta_file(
        indf,
        fasta_seq_ids,
        ofile=opfix + "tcr_contigs.fa",
        polyN_contig_flag=polyN_contig_flag,
    )
    # Clean up
    rm_cmd = "rm -r " + output_dir
    c = subprocess.getstatusoutput(rm_cmd)
    return prod_dub


def generate_metrics_file(
    total_reads,
    summary_df,
    opfix,
    sample_wells=[],
    pass_cells=[],
    isSample=False,
    filtered=False,
    proc_dir="",
    total_unfilt_rds=0,
    isCombine=False,
    total_t_cells=0,
):
    # v1.5 RTK
    # out_file = open(opfix + "_metrics.txt", "w")
    out_file = open(opfix + "_metrics.csv", "w")
    number_of_unfiltered_reads_all = 0
    if isCombine:
        number_of_unfiltered_reads_all = total_unfilt_rds
    else:
        df = pd.read_csv(proc_dir + "/pre_align_stats.csv")
        number_of_unfiltered_reads_all = float(
            df[df["statistic"] == "number_of_reads"]["value"]
        )

    total_calc_reads_all = summary_df["reads_per_cell"].sum()
    # read summary
    if isSample:
        summary_df = summary_df.assign(bc_wind_1=summary_df.index.str.split("_").str[0])
        summary_df = summary_df.astype({"bc_wind_1": "int32"})
        summary_df = summary_df[summary_df["bc_wind_1"].isin(sample_wells)]
        total_calc_reads = summary_df["reads_per_cell"].sum()
        samp_ratio = total_calc_reads / total_calc_reads_all
        number_of_unfiltered_reads = round(
            samp_ratio * number_of_unfiltered_reads_all, 0
        )
        sub_reads_used = summary_df["reads_used"].sum()
        vj_reads = summary_df["reads_vj"].sum()
    else:
        number_of_unfiltered_reads = number_of_unfiltered_reads_all
        total_calc_reads = summary_df["reads_per_cell"].sum()
        sub_reads_used = summary_df["reads_used"].sum()
        vj_reads = summary_df["reads_vj"].sum()

    percent_used = round((float(sub_reads_used) / float(total_calc_reads)) * 100, 2)
    percent_vj = round((float(vj_reads) / float(total_calc_reads)) * 100, 2)
    out_file.write("#Read summary\n")
    out_file.write("tcr_n_unfiltered_reads," + str(number_of_unfiltered_reads) + "\n")
    out_file.write("tcr_n_valid_bc_reads," + str(total_calc_reads) + "\n")
    out_file.write("tcr_valid_bc_reads_VDJ_assembly," + str(sub_reads_used) + "\n")
    out_file.write("tcr_valid_bc_reads_VDJ_overlapping," + str(vj_reads) + "\n")

    # Cell summary
    if filtered == True:
        summary_df = summary_df[summary_df.index.isin(pass_cells)]

    total_cells_with_TCR = summary_df["isTCR"].sum()
    if isCombine:
        if total_t_cells > 0:
            total_cells_with_TCR = total_t_cells
    mean_reads_per_cell = round(
        summary_df[summary_df["isTCR"] == 1]["reads_per_cell"].mean(), 0
    )
    median_vdj = summary_df[summary_df["isTCR"] == 1][
        "VDJ_transcripts_per_cell"
    ].median()

    out_file.write("\n#Cell summary\n")
    out_file.write("tcr_number_of_cells," + str(total_cells_with_TCR) + "\n")
    out_file.write("tcr_mean_reads_per_cell," + str(mean_reads_per_cell) + "\n")
    out_file.write("tcr_median_VDJ_tscp_per_cell," + str(median_vdj) + "\n")

    # prod summary
    tot_proda = summary_df["prod_a"].sum()
    tot_prodb = summary_df["prod_b"].sum()
    tot_prod_pair = summary_df["prod_pair"].sum()

    tra = round(tot_proda / total_cells_with_TCR, 2)
    trb = round(tot_prodb / total_cells_with_TCR, 2)
    pair = round(tot_prod_pair / total_cells_with_TCR, 2)

    out_file.write("tcr_tcells_with_productive_TRA," + str(tot_proda) + "\n")
    out_file.write("tcr_tcells_with_productive_TRB," + str(tot_prodb) + "\n")
    out_file.write("tcr_tcells_with_productive_pair," + str(tot_prod_pair) + "\n")

    # Doublet
    tot_dub = summary_df["isdub"].sum()
    p_dub = round(tot_dub / total_cells_with_TCR, 2)
    out_file.write("tcr_tcells_n_potential_multiplet," + str(tot_dub) + "\n")

    out_file.close()


# # Sample level summaries
# Using functions written by RK for ATAC data in the fastq_ATAC_ano_sep.py script.
# Well definition regex match patterns
# Single well, Colon:range, Dash-range; Get letter and number separate
WELL_1_REGEX = re.compile(r"([ABCDEFGH])([1-9]|1[012])$")
WELL_C_REGEX = re.compile(r"([ABCDEFGH])([1-9]|1[012]):([ABCDEFGH])([1-9]|1[012])$")
WELL_D_REGEX = re.compile(r"([ABCDEFGH])([1-9]|1[012])-([ABCDEFGH])([1-9]|1[012])$")


def parse_96wells(s):
    """Parse well specification string into well indexes

    Return list of well indexes; As 1-based int
    """
    wells = np.arange(96, dtype=int).reshape(8, 12)
    row_letter_to_number = dict(zip(list("ABCDEFGH"), [i for i in range(8)]))
    sub_wells = []
    # Try will fail on non-match that is processed regardless
    try:
        blocks = s.upper().split(",")
        for b in blocks:
            if ":" in b:
                vals = unpack_well_match(WELL_C_REGEX.match(b), 4)
                s_row = row_letter_to_number[vals[0]]
                s_col = vals[1] - 1
                e_row = row_letter_to_number[vals[2]]
                e_col = vals[3]
                sub_wells += list(wells[s_row : e_row + 1, s_col:e_col].flatten())
            elif "-" in b:
                vals = unpack_well_match(WELL_D_REGEX.match(b), 4)
                s_row = row_letter_to_number[vals[0]]
                s_col = vals[1] - 1
                e_row = row_letter_to_number[vals[2]]
                e_col = vals[3] - 1
                sub_wells += list(
                    np.arange(wells[s_row, s_col], wells[e_row, e_col] + 1)
                )
            else:
                vals = unpack_well_match(WELL_1_REGEX.match(b), 2)
                s_row = row_letter_to_number[vals[0]]
                s_col = vals[1] - 1
                sub_wells += [wells[s_row, s_col]]
        sub_wells = list(np.unique(sub_wells))
    except Exception as e:
        # print(e)
        pass

    # Convert to 1-base
    sub_wells = [i + 1 for i in sub_wells]
    return sub_wells


def unpack_well_match(mat, num):
    """Unpack well regex match with two or four items [letter, well, [letter, well]]"""
    vals = []
    if num == 2:
        vals = [mat[1], int(mat[2])]
    elif num == 4:
        vals = [mat[1], int(mat[2]), mat[3], int(mat[4])]
    return vals


# NOTE: Does not check if guessed kit makes sense
def get_sample_wells(infile, intype="gfile"):  # bc_info removed
    """Get output well/cell subset collection

    Return dictionary with well ids for each sample
    """
    if intype == "gfile":
        # Collect subset defs; Sample wells or cells
        def_list = []
        # Samples will only look at first (one) barcode
        with open(infile) as INFILE:
            for line in INFILE:
                # Ignore '#' commented out / blank lines
                parts = line.strip().split("#")[0].split()
                # Def should have <name> <wells/cells>
                if len(parts) >= 2:
                    def_list.append(parts[:2])

        # Should have something
        if not def_list:
            return None

        # Collect objects with parsed well/cell identifiers
        # NOTE no max wells check for now
        # max_wells = bc_info['num_wells'][1]
        oset = {}
        for name, s_def in def_list:
            bc_lis = parse_96wells(s_def)
            if not bc_lis:
                print(f"Problem parsing {name} '{s_def}' as groups (wells)")
                return None
            # Too many?
            # if max(bc_lis) > max_wells:
            # print(f"Problem parsing {name} '{s_def}'; Kit set for {max_wells} wells ({max(bc_lis)})")
            # return None

            oset[name] = bc_lis
    elif intype == "wfile":
        oset = {}
        with open(infile) as INFILE:
            for line in INFILE:
                # Ignore '#' commented out / blank lines
                parts = line.strip().split("#")[0].split()
                # Def should have <name> <wells/cells>
                if len(parts) >= 2:
                    name = parts[0]
                    bc_lis2 = parts[1].split(",")
                    bc_lis = [int(i) for i in bc_lis2]
                    oset[name] = bc_lis
    return oset


def get_sample_level_files(
    igblast_df,
    sample_dict,
    main_odir,
    polyN_contig_flag,
    summary_df,
    total_reads,
    nthreads,
    pass_cells=[],
    filtered=False,
    total_unfilt_rds=0,
    isCombine=False,
):
    track_samples = {}
    for sample_name, sample_wells in sample_dict.items():
        if sample_name in track_samples:
            print("Duplicate sample name. Skipping.")
            continue
        track_samples[sample_name] = 1
        igblast_df_samp = igblast_df.assign(
            bc_wind_1=igblast_df["bc_wind"].str.split("_").str[0]
        )
        igblast_df_samp = igblast_df_samp.astype({"bc_wind_1": "int32"})
        igblast_comb_samp_subset = igblast_df_samp[
            igblast_df_samp["bc_wind_1"].isin(sample_wells)
        ]
        all_data_ids = igblast_df_samp.bc_wind_1.tolist()
        # Add check for all-well
        if igblast_comb_samp_subset.shape[0] == 0:
            print("No outputs for sample: " + sample_name)
        elif set(all_data_ids) == set(sample_wells):
            # If the wells are the same as all data, it is all-wells. So just copying files instead of re running.
            if filtered:
                mk_cmd = "mkdir -p " + main_odir + "/" + sample_name + "/TCR_filtered/"
                c = subprocess.getstatusoutput(mk_cmd)
                mk_cmd = "mkdir -p " + main_odir + "/" + sample_name + "/report/"
                c = subprocess.getstatusoutput(mk_cmd)
                cp_cmd = (
                    "cp "
                    + main_odir
                    + "/process/scratch/TCR_filtered/* "
                    + main_odir
                    + "/"
                    + sample_name
                    + "/TCR_filtered/"
                )
                c = subprocess.getstatusoutput(cp_cmd)
                mv_cmd = (
                    "mv "
                    + main_odir
                    + "/"
                    + sample_name
                    # v1.5 RTK
                    # + "/TCR_filtered/tcr_filt_metrics.txt "
                    + "/TCR_filtered/tcr_filt_metrics.csv "
                    + main_odir
                    + "/"
                    + sample_name
                    + "/report/"
                )
                c = subprocess.getstatusoutput(mv_cmd)

            else:
                mk_cmd = (
                    "mkdir -p " + main_odir + "/" + sample_name + "/TCR_unfiltered/"
                )
                c = subprocess.getstatusoutput(mk_cmd)
                mk_cmd = "mkdir -p " + main_odir + "/" + sample_name + "/report/"
                c = subprocess.getstatusoutput(mk_cmd)
                cp_cmd = (
                    "cp "
                    + main_odir
                    + "/process/scratch/TCR_unfiltered/* "
                    + main_odir
                    + "/"
                    + sample_name
                    + "/TCR_unfiltered/"
                )
                c = subprocess.getstatusoutput(cp_cmd)
                mv_cmd = (
                    "mv "
                    + main_odir
                    + "/"
                    + sample_name
                    # v1.5 RTK
                    # + "/TCR_unfiltered/tcr_unfilt_metrics.txt "
                    + "/TCR_unfiltered/tcr_unfilt_metrics.csv "
                    + main_odir
                    + "/"
                    + sample_name
                    + "/report/"
                )
                c = subprocess.getstatusoutput(mv_cmd)
        else:
            if filtered:
                mk_cmd = "mkdir -p " + main_odir + "/" + sample_name + "/TCR_filtered/"
                c = subprocess.getstatusoutput(mk_cmd)
                mk_cmd = "mkdir -p " + main_odir + "/" + sample_name + "/report/"
                c = subprocess.getstatusoutput(mk_cmd)
                prod_dub_info = generate_summary_files_threaded(
                    igblast_comb_samp_subset,
                    opfix=main_odir + "/" + sample_name + "/TCR_filtered/",
                    polyN_contig_flag=polyN_contig_flag,
                    nthreads=nthreads,
                    proc_dir=main_odir + "/process/tcr",
                )
                generate_metrics_file(
                    total_reads,
                    summary_df,
                    proc_dir=main_odir + "/process",
                    opfix=main_odir + "/" + sample_name + "/report/tcr_filt",
                    sample_wells=sample_wells,
                    pass_cells=pass_cells,
                    isSample=True,
                    filtered=True,
                    total_unfilt_rds=total_unfilt_rds,
                    isCombine=isCombine,
                )
            else:
                mk_cmd = (
                    "mkdir -p " + main_odir + "/" + sample_name + "/TCR_unfiltered/"
                )
                c = subprocess.getstatusoutput(mk_cmd)
                mk_cmd = "mkdir -p " + main_odir + "/" + sample_name + "/report/"
                c = subprocess.getstatusoutput(mk_cmd)
                prod_dub_info = generate_summary_files_threaded(
                    igblast_comb_samp_subset,
                    opfix=main_odir + "/" + sample_name + "/TCR_unfiltered/",
                    polyN_contig_flag=polyN_contig_flag,
                    nthreads=nthreads,
                    proc_dir=main_odir + "/process/tcr",
                )
                generate_metrics_file(
                    total_reads,
                    summary_df,
                    proc_dir=main_odir + "/process",
                    opfix=main_odir + "/" + sample_name + "/report/tcr_unfilt",
                    sample_wells=sample_wells,
                    isSample=True,
                    filtered=False,
                    total_unfilt_rds=total_unfilt_rds,
                    isCombine=isCombine,
                )


################################################
# COMBINE MODE


def combine_mode_summary_file(total_reads, summary_df, prod_dub_info, output_dir):
    # subset summary_df
    summary_df = summary_df[
        [
            "reads_per_cell",
            "isTCR",
            "VDJ_transcripts_per_cell",
            "reads_used",
            "reads_vj",
        ]
    ]
    if total_reads != summary_df.reads_per_cell.sum():
        print("Total reads not matching sum of reads per cell in make_summary_file.")

    # get productive info and doublet info
    prod_dub_df = pd.DataFrame.from_dict(
        prod_dub_info,
        orient="index",
        columns=["prod_a", "prod_b", "prod_pair", "inf_prod_pair", "isdub"],
    )
    summary_df = summary_df.join(prod_dub_df)

    # fill na
    summary_df = summary_df.fillna(0)

    #
    out_sum_file = output_dir + "/process/tcr_barcode_summary.txt"
    summary_df.to_csv(
        out_sum_file, sep="\t", header=True, index=True, index_label="barcode"
    )
    return summary_df


def run_combine(
    sub_file, output_dir, sample_flag, sampfile, sample_file_type, nthreads
):
    # Checks and see if there was filtered data generated.
    unfilt_list = []
    filt_list = []
    total_sublibraries = 0
    with open(sub_file, "r") as inFile:
        for aline in inFile.readlines():
            aline = aline.rstrip()
            if aline == "":
                continue
            total_sublibraries += 1
            # v1.5 RTK
            #ufile = aline + "/process/scratch/TCR_unfiltered/tcr_unfilt_metrics.txt"
            ufile = aline + "/process/scratch/TCR_unfiltered/tcr_unfilt_metrics.csv"
            uchk = os.path.isfile(ufile)
            if uchk:
                uschk = os.stat(ufile).st_size == 0
                if uschk:
                    print(ufile + " is empty. Skipping sublibrary.")
                else:
                    unfilt_list.append(aline)
            else:
                print(ufile + " does not exist. Skipping sublibrary.")

            # v1.5 RTK
            #ffile = aline + "/process/scratch/TCR_filtered/tcr_filt_metrics.txt"
            ffile = aline + "/process/scratch/TCR_filtered/tcr_filt_metrics.csv"
            fchk = os.path.isfile(ffile)
            if fchk:
                fschk = os.stat(ffile).st_size == 0
                if fschk:
                    print(ffile + " is empty. Skipping sublibrary.")
                else:
                    filt_list.append(aline)
            else:
                print(ffile + " does not exist. Skipping sublibrary.")

    if len(unfilt_list) == 0:
        print("No sublibraries to combine.")
        return True

    uflag = 0
    fflag = 0
    if len(unfilt_list) > 0:
        uflag = 1
    if len(filt_list) > 0:
        fflag = 1

    # Summary df and total reads
    sub_dict = {}
    sub_cnter = 0
    summ_all = pd.DataFrame()
    total_unfiltered_reads = 0
    total_vb = 0

    for aline in unfilt_list:
        sub_cnter += 1
        sub_suffix = "s" + str(sub_cnter)
        if aline not in sub_dict:
            sub_dict[aline] = sub_suffix
        summ_file = aline + "/process/tcr_barcode_summary.txt"
        if os.path.isfile(summ_file):
            summ_df = pd.read_csv(summ_file, sep="\t")
            summ_df["barcode"] = summ_df["barcode"].astype(str) + "__" + sub_suffix
            summ_all = pd.concat([summ_all, summ_df])
        else:
            print(summ_file, " not present")
        # v1.5 RTK
        #met_file = aline + "/process/scratch/TCR_unfiltered/tcr_unfilt_metrics.txt"
        met_file = aline + "/process/scratch/TCR_unfiltered/tcr_unfilt_metrics.csv"
        if os.path.isfile(met_file):
            with open(met_file, "r") as metFile:
                for bline in metFile.readlines():
                    bline = bline.rstrip()
                    if re.match("^tcr_n_unfiltered_reads.*", bline) is not None:
                        total_unfiltered_reads = total_unfiltered_reads + float(
                            bline.split(",")[1]
                        )
                    if re.match("^tcr_n_valid_bc_reads.*", bline) is not None:
                        total_vb = total_vb + float(bline.split(",")[1])
        else:
            print(met_file, " not present")
    # Get total filtered cells
    total_filtered_t_cells = 0
    for aline in filt_list:
        # v1.5 RTK
        #met_file = aline + "/process/scratch/TCR_filtered/tcr_filt_metrics.txt"
        met_file = aline + "/process/scratch/TCR_filtered/tcr_filt_metrics.csv"
        if os.path.isfile(met_file):
            with open(met_file, "r") as metFile:
                for bline in metFile.readlines():
                    bline = bline.rstrip()
                    if re.match("^tcr_number_of_cells.*", bline) is not None:
                        total_filtered_t_cells = total_filtered_t_cells + float(
                            bline.split(",")[1]
                        )
        else:
            print(met_file, " not present")

    summ_all = summ_all.set_index("barcode")
    # Combine
    polyN_contig_flag = 0
    summary_df = pd.DataFrame()
    if uflag == 1:
        print("Combining " + str(len(unfilt_list)) + " unfiltered sublibraries.")
        airr_all = pd.DataFrame()

        for aline in unfilt_list:
            aline = aline.rstrip()
            airr_file = (
                aline + "/process/scratch/TCR_unfiltered/tcr_annotation_airr.tsv"
            )
            if os.path.isfile(airr_file):
                airr_df = pd.read_csv(airr_file, sep="\t")
                airr_df["sequence_id"] = (
                    airr_df["sequence_id"].astype(str) + "__" + sub_dict[aline]
                )
                airr_df["cell_barcode"] = (
                    airr_df["cell_barcode"].astype(str) + "__" + sub_dict[aline]
                )
                airr_all = pd.concat([airr_all, airr_df])
            else:
                print(airr_file + " not present")

        if "transcript_count" in airr_all.columns:
            airr_all = airr_all.rename(
                columns={
                    "cdr3": "junction",
                    "cdr3_aa": "junction_aa",
                    "read_count": "read_support",
                    "cell_barcode": "bc_wind",
                    "transcript_count": "n10_supp",
                }
            )
            polyN_contig_flag = 1
        else:
            airr_all = airr_all.rename(
                columns={
                    "cdr3": "junction",
                    "cdr3_aa": "junction_aa",
                    "read_count": "read_support",
                    "cell_barcode": "bc_wind",
                }
            )

        mk_cmd = "mkdir -p " + output_dir + "/process/scratch/TCR_unfiltered/"
        c = subprocess.getstatusoutput(mk_cmd)
        prod_dub_info = generate_summary_files_threaded(
            airr_all,
            opfix=output_dir + "/process/scratch/TCR_unfiltered/",
            polyN_contig_flag=polyN_contig_flag,
            nthreads=nthreads,
            proc_dir=output_dir + "/process/tcr",
        )
        u_summary_df = combine_mode_summary_file(
            total_reads=total_vb,
            summary_df=summ_all,
            prod_dub_info=prod_dub_info,
            output_dir=output_dir,
        )
        summary_df = u_summary_df
        tmp = generate_metrics_file(
            total_reads=total_vb,
            summary_df=u_summary_df,
            opfix=output_dir + "/process/scratch/TCR_unfiltered/tcr_unfilt",
            isSample=False,
            filtered=False,
            total_unfilt_rds=total_unfiltered_reads,
            isCombine=True,
        )
        if sample_flag == 1:
            print("Generating combined unfiltered sample summary files")
            sample_dict = get_sample_wells(sampfile, intype=sample_file_type)
            tmp = get_sample_level_files(
                airr_all,
                sample_dict,
                main_odir=output_dir,
                polyN_contig_flag=polyN_contig_flag,
                summary_df=u_summary_df,
                total_reads=total_vb,
                filtered=False,
                total_unfilt_rds=total_unfiltered_reads,
                isCombine=True,
                nthreads=nthreads,
            )

    if fflag == 1:
        print("Combining " + str(len(filt_list)) + " filtered sublibraries.")
        airr_all = pd.DataFrame()
        for aline in filt_list:
            airr_file = aline + "/process/scratch/TCR_filtered/tcr_annotation_airr.tsv"
            if os.path.isfile(airr_file):
                airr_df = pd.read_csv(airr_file, sep="\t")
                airr_df["sequence_id"] = (
                    airr_df["sequence_id"].astype(str) + "__" + sub_dict[aline]
                )
                airr_df["cell_barcode"] = (
                    airr_df["cell_barcode"].astype(str) + "__" + sub_dict[aline]
                )
                airr_all = pd.concat([airr_all, airr_df])
            else:
                print(airr_file + " not present")

        if "transcript_count" in airr_all.columns:
            airr_all = airr_all.rename(
                columns={
                    "cdr3": "junction",
                    "cdr3_aa": "junction_aa",
                    "read_count": "read_support",
                    "cell_barcode": "bc_wind",
                    "transcript_count": "n10_supp",
                }
            )
            polyN_contig_flag = 1
        else:
            airr_all = airr_all.rename(
                columns={
                    "cdr3": "junction",
                    "cdr3_aa": "junction_aa",
                    "read_count": "read_support",
                    "cell_barcode": "bc_wind",
                }
            )

        mk_cmd = "mkdir -p " + output_dir + "/process/scratch/TCR_filtered/"
        c = subprocess.getstatusoutput(mk_cmd)
        prod_dub_info = generate_summary_files_threaded(
            airr_all,
            opfix=output_dir + "/process/scratch/TCR_filtered/",
            polyN_contig_flag=polyN_contig_flag,
            nthreads=nthreads,
            proc_dir=output_dir + "/process/tcr",
        )
        pass_cells = list(set(airr_all["bc_wind"].tolist()))
        if uflag == 1:
            u_summary_df = summary_df
        else:
            u_summary_df = combine_mode_summary_file(
                total_reads=total_vb,
                summary_df=summ_all,
                prod_dub_info=prod_dub_info,
                output_dir=output_dir,
            )
        tmp = generate_metrics_file(
            total_reads=total_vb,
            summary_df=u_summary_df,
            opfix=output_dir + "/process/scratch/TCR_filtered/tcr_filt",
            isSample=False,
            filtered=True,
            pass_cells=pass_cells,
            total_unfilt_rds=total_unfiltered_reads,
            isCombine=True,
            total_t_cells=total_filtered_t_cells,
        )
        if sample_flag == 1:
            print("Generating combined filtered sample summary files")
            sample_dict = get_sample_wells(sampfile, intype=sample_file_type)
            tmp = get_sample_level_files(
                airr_all,
                sample_dict,
                main_odir=output_dir,
                polyN_contig_flag=polyN_contig_flag,
                summary_df=u_summary_df,
                total_reads=total_vb,
                pass_cells=pass_cells,
                filtered=True,
                total_unfilt_rds=total_unfiltered_reads,
                isCombine=True,
                nthreads=nthreads,
            )


################################################
# # Transcriptome filtering


def process_cell_metadata(transcriptome_file):
    tscp_file = pd.read_csv(transcriptome_file, sep=",")
    return tscp_file.bc_wells.unique()


def get_tscp_cells_only(indf, transcriptome_file):
    keep_cells = process_cell_metadata(transcriptome_file)
    indf = indf[indf["bc_wind"].isin(keep_cells)]
    return indf, keep_cells


################################################
# Other functions


def explain_gfile_details():
    story = f"""
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
Description

    Analyze TCR seq data. This script analyzes once the "barcode_head.fastq.gz" file is generated. 

    Groups may be specified via a list file (--gfile <list>).
    
    For group specification via file (--gfile <list>), a simple, space delimited 
    text list file is used. Each file line should have:
        '<name1> <wells>'.
        '<name2> <other-wells>'.
        
    Wells are specified in blocks, ranges, or individually like this:
        'A1:C6' specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6.
        'A1-B6' specifies a range as [start]-[end]; A1-A12, B1-6.
        'C4' specifies a single well.
        Multiple selections are joined by commas (no space), e.g. 'A1-A6,B1:D3,C4'
"""
    print(story)


def check_infile(fname, verb=True, toxic=False):
    """Check file

    Return filename, possibly expanded if env var
    """
    r_fname = ""
    # Expand any env vars; User for ~/ and vars for $home
    fname = os.path.expanduser(fname)
    fname = os.path.expandvars(fname)
    # Read access means input is good
    if os.access(fname, os.R_OK):
        r_fname = fname
    else:
        story = f"Cannot read file: '{fname}'"
        if toxic:
            raise IOError(story)
        if verb:
            print(story)
    return r_fname


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ok = main()
