#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Make reference genome related stuff
#

import re
import gzip
import pandas as pd
import numpy as np
from collections import defaultdict

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import crispr
    import genes
    import utils
else:
    from splitpipe import crispr
    from splitpipe import genes
    from splitpipe import utils


# ---------------------------------------------------------------------------
def have_mkref_ins(spipe, verb=True):
    """Check if Mkref input files exist

    Return True/False, list of files
    """
    ok = True
    check_files = []

    # CRISPR case
    if spipe.is_focal_crispr():
        guides = spipe.get_par_val("crsp_guides", as_path=True)
        check_files.append(guides)
        if not utils.check_infile(check_files, verb=verb):
            ok = False
    else:
        # Normal genome inputs = three lists: genome_name + genes + fasta
        genome_names = spipe.get_par_val("genome_name", as_list=True)
        gtf_files = spipe.get_par_val("genes", as_list=True)
        fasta_files = spipe.get_par_val("fasta", as_list=True)
        nn = len(genome_names)
        ng = len(gtf_files)
        nf = len(fasta_files)
        if not (nn == ng == nf):
            # Add these (non-files) to check_files for possible reporting
            check_files.append(f"genome_name {genome_names}")
            check_files.append(
                f"PROBLEM: length of genome_name {nn}, genes {ng}, fasta {nf} differ"
            )
            ok = False

        # Add possible gfasta files to list then check files
        if ok:
            gfasta = spipe.get_par_val("gfasta", as_list=True)
            gfasta_files = []
            if gfasta:
                # List of [name, fname]
                for _, fname in gfasta:
                    gfasta_files.append(fname)
                # If list is good, input files all good
                check_files = gtf_files + fasta_files + gfasta_files
                if not utils.check_infile(check_files, verb=verb):
                    ok = False

    return ok, check_files


def run_mkref(spipe):
    """Run mkref step"""
    spipe.report_proc_step("Mkref")
    ok = True
    # Only input checked; Run regardless of output status
    i_ok, i_list = have_mkref_ins(spipe)
    # Run only if input
    if not i_ok:
        bad_list = utils.bad_infile_list(i_list)
        story = f"Don't have inputs for mkref: {bad_list}"
        spipe.set_problem(story)
        ok = False
    else:
        ok = make_ref_genome(spipe)

    ok = spipe.no_problems()
    spipe.report_proc_step("Mkref", status=ok)
    return ok


def make_ref_genome(spipe):
    """Make reference genome, preprocessing inputs then run indexing

    Return status
    """
    ok = True
    gfasta_info = gtf_info = None
    # Assume will be processing gtf and fasta for STAR
    gtf_fasta = True

    # fb/crispr guide processing case
    if spipe.is_focal_crispr():
        # Either guide direct matching or prep fasta inputs (e.g. for STAR)
        ok = crispr.mkref_prep_db(spipe)
        # Only prepare gtf and fasta if using STAR
        if spipe.is_focal_crispr(map_direct=True):
            gtf_fasta = False

    # Process any gfasta files (into fasta, gtf)
    if ok and gtf_fasta:
        ok, gfasta_info, gtf_info = prep_fasta_gtf_inputs(spipe)

    # Make a gene stats data table
    if ok:
        spipe.report_run_story("Generating gene data tables")
        if not make_gene_tables(spipe):
            ok = False
            story = "Failed to make gene data tables"
            spipe.set_problem(story)

    # Index genome; Use gtf_fasta flag for special crispr STAR case
    if ok and gtf_fasta:
        spipe.report_run_story("Creating genome index")
        if not generate_STAR_index(spipe):
            ok = False
            spipe.set_problem("Failed to generate STAR index")

    # Write info to file
    if ok:
        write_mkref_env_files(spipe, gtf_info, gfasta_info)
    else:
        spipe.report_run_story("Issues, so mkref def not saved")

    # Clean up
    if ok:
        spipe.report_run_story("Cleaning up")
        if not clean_up(spipe):
            # Don't need to fail, but still report problem
            spipe.set_problem("Failed mkref clean")

    return ok


def prep_fasta_gtf_inputs(spipe):
    """Prepare fasta and gtf (and gfasta) inputs for ref db buidling

    Auxiliary info from gfasta processing and gtf processing returned
    as pair of dicts (or None)

    Return (status, dict, dict)
    """
    ok = True
    gfasta_info = gtf_info = None

    # Process any gfasta files
    ok, gfasta_info = prep_gfasta_files(spipe)
    if not ok:
        story = "Failed to process gfasta record(s)"
        spipe.set_problem(story)
    # Generate genome sequence input file
    if ok:
        spipe.report_run_story("Processing fasta files")
        n_chr, n_base = make_combined_fasta(spipe, gfasta_info)
        if not n_chr:
            story = "Failed to make combined genome fasta"
            spipe.set_problem(story)
            ok = False
    spipe.report_run_story2(f"Total chromosomes {n_chr}, bases {n_base}")
    # Make a combine annotation from GTF
    if ok:
        spipe.report_run_story("Processing gene annotation files")
        gtf_info = make_gtf_annotations(spipe, gfasta_info)
        if not gtf_info:
            ok = False
            story = "Failed to make gtf annotations"
            spipe.set_problem(story)
    # Check fasta and gtf info is consistent
    if ok:
        spipe.report_run_story("Checking fasta gtf overlap")
        if not fasta_gtf_over_ok(spipe):
            ok = False
            story = "Issue(s) with fasta and gtf overlap"
            spipe.set_problem(story)

    return ok, gfasta_info, gtf_info


def prep_gfasta_files(spipe):
    """Prepare any gfasta-based input files

    Return tuple (status, dict) with info about processing
    """
    gfasta = spipe.get_par_val("gfasta", as_list=True)
    # Nothing to prepare
    if not gfasta:
        return True, None

    spipe.report_run_story(
        f"Preparing fasta and gtf files from {len(gfasta)} gfasta records"
    )
    # Collect filenames
    fasta_files = []
    gtf_files = []
    genome_list = []
    n_rec_list = []
    bname = spipe.filepath("GENO_GFASTA_PREF", None, genome_dir="output_dir")
    for i, (gname, fname) in enumerate(gfasta):
        fas_list = read_one_fasta(spipe, fname)
        if not fas_list:
            return False, None

        # Write out and save filenames
        fas_name = bname + f"_seqs{i+1}.fas"
        gtf_name = bname + f"_genes{i+1}.gtf"
        write_seq_list_fasta(fas_list, fas_name)
        n_rec = write_seq_list_gtf(fas_list, gtf_name)
        fasta_files.append(fas_name)
        gtf_files.append(gtf_name)
        genome_list.append(gname)
        n_rec_list.append(n_rec)

    new_dict = {
        "genome": genome_list,
        "fasta": fasta_files,
        "gtf": gtf_files,
        "n_records": n_rec_list,
    }

    return True, new_dict


def read_one_fasta(spipe, fname):
    """Load fasta file content

    Return list of (name, seq) tuples
    """
    try:
        seq_list = utils.fasta_name_seq_list(fname, one_name=False)
        return seq_list
    except Exception as e:
        story = f"Bad sample list file: {e} (exception)"
        spipe.set_problem(story)
        return None


def write_seq_list_fasta(seq_list, fname, line_len=100):
    """Write out (gfasta-based) fasta file

    seq_list = list of (name, seq) tuples

    Returns number of records
    """
    n_rec = 0
    with open(fname, "w") as OUTFILE:
        for name, seq in seq_list:
            n_rec += 1
            print(f">{name}", file=OUTFILE)
            spos = 0
            while spos < len(seq):
                print(seq[spos : spos + line_len], file=OUTFILE)
                spos += line_len
    return n_rec


def write_seq_list_gtf(seq_list, fname, biotype="gfasta"):
    """Write out (gfasta-based) gtf file

    seq_list = list of (name, seq) tuples

    Writes two rows per record; One with feature 'gene' and one with 'exon'

    Returns number of records
    """
    n_rec = 0
    with open(fname, "w") as OUTFILE:
        for name, seq in seq_list:
            n_rec += 1
            ats = f'gene_id "gene_id_{n_rec}"; gene_name "{name}"'
            if biotype:
                ats += f'; gene_biotype "{biotype}"'
            # Col 1: Seqname = fasta seq record name (i.e. chromosome)
            # Col 2: Source = gfasta processing
            # Col 3: Feature = 'gene' and 'exon'
            # Col 4-5: Coords = full length of fasta seq
            # Col 6-8: Score, strand, frame = generic
            # Col 9: Attributes string
            row = [
                name,
                "split-pipe-gfasta",
                "gene",
                "1",
                str(len(seq)),
                ".",
                "+",
                ".",
                ats,
            ]
            print("\t".join(row), file=OUTFILE)
            # Swap out feature from gene to exon and write again; All else identical
            row[2] = "exon"
            print("\t".join(row), file=OUTFILE)
    return n_rec


def make_combined_fasta(spipe, gfasta_info):
    """Process genome fasta input(s) into one file

    This is needed in case multi-species; Add species to chr names
    Also save chromosome names and lengths to file

    Return number of chromosomes, total length
    """
    n_chrom = tot_len = 0
    # Lists of fasta and names
    genome_names = spipe.get_par_val("genome_name", as_list=True)
    fasta_fnames = spipe.get_par_val("fasta", as_list=True)
    # If gfasta, add those files
    if gfasta_info:
        genome_names = genome_names + gfasta_info["genome"]
        fasta_fnames = fasta_fnames + gfasta_info["fasta"]

    # Output; For mkref, 'output_dir' is 'genome_dir'
    out_fname = spipe.filepath("GENO_PROC_FASTA_RAW", None, genome_dir="output_dir")
    OUTFILE = open(out_fname, "w")

    # Collect chromosome names for file
    ch_names = []
    ch_lens = []
    ch_length = 0

    # Loop over input fastas
    for fasta, gname in zip(fasta_fnames, genome_names):
        spipe.report_run_story2(f"Processing {gname} {fasta}")
        if fasta.endswith(".gz"):
            is_gz = True
            INFILE = gzip.open(fasta, "rb")
        else:
            INFILE = open(fasta)
            is_gz = False

        # All lines
        while True:
            line = INFILE.readline()
            if not line:
                break
            if is_gz:
                line = line.decode()
            # Ignore blank or comment lines
            line = line.split("#")[0].strip()
            if not line:
                continue
            # Header line gets reformat
            if line[0] == ">":
                #  input like '>X some more stories...'
                # output like '>genome_X some more stories...'
                parts = line[1:].split()
                rest_of = " ".join(parts[1:])
                chrom = f"{gname}_{parts[0]}"
                print(f">{chrom} {rest_of}", file=OUTFILE)
                # Counts, name and length
                n_chrom += 1
                ch_names.append(chrom)
                if ch_length:
                    ch_lens.append(ch_length)
                    ch_length = 0
            else:
                print(line, file=OUTFILE)
                tot_len += len(line)
                ch_length += len(line)
        INFILE.close()
    OUTFILE.close()

    # Save last length
    if ch_length:
        ch_lens.append(ch_length)
    # Save chromosome info for fasta part
    df = pd.DataFrame({"chrom": ch_names, "fasta": ch_lens})
    df.set_index("chrom", inplace=True)
    spipe.write_df(df, "GENO_CHROM_DATA", genome_dir="output_dir")

    return n_chrom, tot_len


def split_attributes(s):
    """Returns a dictionary of GTF/GFF file attributes

    s = string with gtf column 9 data, i.e. attributes. For example:
        'gene_id "ENSG00000188976"; gene_version "10"; transcript_id "ENST00000327044"; ...'

    Return dict
    """
    # Break into space-stripped key val strings
    att_list = [p.strip() for p in s.strip().split("; ")]
    # Keys = first space-delimited token
    att_keys = [a.split(" ")[0] for a in att_list]
    # Values = remaining text
    att_values = [a.split(" ", 1)[1] for a in att_list]
    # Attributes may or may not be quoted; Remove this and strip end space or semicolon
    att_values = [a.replace('"', "").strip().strip(";") for a in att_values]

    return dict(zip(att_keys, att_values))


def get_attribute(s, att):
    """Look for attribute value by key in attributes string

    Creates dict of attribute fields then looks in that for value
    (Not the most efficient way to process ... but works with pandas df apply)

    Return attribute value if found
    """
    att_value = ""
    try:
        att_value = split_attributes(s)[att]
    except Exception:
        pass
    return att_value


def get_one_att_value(att_string, key):
    """Attempt to get value for one attribute with given key

    att_string = "Attribues" data from gtf line; i.e. col 9
    Attributes are semicolon-separated key "value" pairs, like:
        'gene_id "Vmn2r116l-ps30"; transcript_id "";  ...  ; db_xref "GeneID:102552844"; '

    Return str
    """
    key += " "
    if att_string.startswith(key):
        start_index = len(key)
    else:
        index = att_string.find("; " + key)
        if index < 0:
            return ""
        start_index = index + len(key) + 2

    # End of key-value ... No error if not found (bogus format)
    end_index = att_string.find(";", start_index)
    if end_index < 0:
        return ""
    # Get rid of quotes and any trailing spaces
    return att_string[start_index:end_index].replace('"', "").strip()


def make_gtf_annotations(spipe, gfasta_info):
    """Make custom format gtf file

    Build gene_info structure and save to file
    Also save chrom info to file

    Return dict with processing info, or None
    """
    splicing = spipe.get_par_val("mkref_star_splicing", as_bool=True)
    min_exons = spipe.get_par_val("mkref_min_exons", as_int=True)
    min_genes = spipe.get_par_val("mkref_min_genes", as_int=True)
    # For format guessing
    gtf_guess_min = spipe.get_par_val("mkref_gtf_guess_min", as_int=True)

    # Genome names and gtf_filenames are lists
    genome_names = spipe.get_par_val("genome_name", as_list=True)
    gtf_filenames = spipe.get_par_val("genes", as_list=True)
    # Default min guess numbers for given gtf files (not gfasta)
    guess_mins = [gtf_guess_min] * len(gtf_filenames)

    # If no regular gtf files, need to reset min gene and exon for gfasta
    reset_min = True if not genome_names else False

    if gfasta_info:
        genome_names = genome_names + gfasta_info["genome"]
        gtf_filenames = gtf_filenames + gfasta_info["gtf"]
        # Set min guess, possibly shorten to record length
        for n_rec in gfasta_info["n_records"]:
            gtf_guess_min = min(n_rec, gtf_guess_min)
            guess_mins.append(gtf_guess_min)

    # Reset gene and exon min if no regular gtf; Only gfasta numbers are in list
    if reset_min:
        min_recs = min(guess_mins)
        min_exons = min(min_exons, min_recs)
        min_genes = min(min_genes, min_recs)

    spipe.report_run_story(f"Loading {len(genome_names)} listed gtfs")
    # Load gtf records; Merge together if same-name genome
    gtf_dict = {}
    for fname, genome, gtf_guess_min in zip(gtf_filenames, genome_names, guess_mins):
        spipe.report_run_story2(f"Reading GTF {fname}")

        # Load gtf dataframe; Doesn't have header or index col; Set cols after load
        gtf_df = utils.read_csv(
            fname, sep="\t", dtype="object", index_col=None, header=None
        )
        gtf_df.columns = genes.GTF_COLS

        # All rows need gene_id attribute
        spipe.report_run_story2(f"Checking minimal attributes")
        ok, story = check_gtf_att_data_df(gtf_df, key="gene_id")
        if not ok:
            spipe.set_problem(f"Unable to process gtf {fname}")
            spipe.set_problem(story)
            return None

        # Get info on format flavor
        gtf_info, gtf_story = genes.gtf_guess_flavor(gtf_df, gtf_guess_min)
        if not gtf_info or not gtf_info["trust"]:
            spipe.set_problem(f"Unable to process gtf {fname}")
            spipe.set_problem(gtf_story)
            # print("<<< NONE")
            return None
        grec = genes.GtfData(fname, name=genome, filename=fname, gtf_info=gtf_info)
        spipe.report_run_story2(
            f"{genome} GTF {fname}, format flavor {grec.get_flavor()}, {grec.num_records()} records"
        )
        # If already have this genome name, then merge
        if grec.get_name() in gtf_dict:
            gtf_dict[grec.get_name()].merge_df(grec)
        else:
            gtf_dict[grec.get_name()] = grec

    # print("GG4... mins", min_exons, min_genes)
    # Generate a combined GTF with only the exon annotations
    gtf_exon_combined = get_gtf_feature_df(gtf_dict, exon=True)
    spipe.report_run_story2(f"gtf_exon_combined {gtf_exon_combined.shape}")
    if len(gtf_exon_combined) < min_exons:
        story = f"Too few gtf exons {len(gtf_exon_combined)} (min {min_exons})"
        spipe.set_problem(story)
        return None
    write_geno_df(spipe, gtf_exon_combined, "GENO_PROC_EXON_GTF", sep="\t")

    # Generate a combined GTF with only the gene annotations
    if splicing:
        gtf_gene_combined = get_gtf_feature_df(gtf_dict, exon=False)
        spipe.report_run_story2(
            f"gtf_gene_combined, slicing True {gtf_gene_combined.shape}"
        )
    else:
        # If not splicing simply copy and rename feature exon >--> gene
        gtf_gene_combined = gtf_exon_combined.copy(deep=True)
        gtf_gene_combined["Feature"] = "gene"
        spipe.report_run_story2(
            f"gtf_gene_combined, slicing False {gtf_gene_combined.shape}"
        )
    if len(gtf_gene_combined) < min_genes:
        story = f"Too few genes {len(gtf_gene_combined)} (min {min_genes})"
        spipe.set_problem(story)
        return None
    write_geno_df(spipe, gtf_gene_combined, "GENO_PROC_GENE_GTF", sep="\t")

    # Get locations of genes. We are using the longest possible span of different transcripts here
    gtf_gene_combined.loc[:, "gene_id"] = gtf_gene_combined.Attributes.apply(
        lambda s: get_attribute(s, "gene_id")
    )

    # gene_starts = gtf_gene_combined.groupby("gene_id").Start.apply(min)
    # gene_ends = gtf_gene_combined.groupby("gene_id").End.apply(max)
    gene_starts = gtf_gene_combined.groupby("gene_id").Start.apply(np.minimum.reduce)
    gene_ends = gtf_gene_combined.groupby("gene_id").End.apply(np.maximum.reduce)

    chroms = gtf_gene_combined.groupby("gene_id").Chromosome.apply(lambda s: list(s)[0])
    strands = gtf_gene_combined.groupby("gene_id").Strand.apply(lambda s: list(s)[0])

    # Note; Only shift dicts *after* gtf has been written out (above)
    # Coord shift from 1-based GTF coords to 0-based pysam
    gene_starts = gene_starts - 1
    gene_ends = gene_ends - 1

    # Get stepsize for processing new genome; Saved with gene_info
    gtf_dict_stepsize = spipe.get_par_val("mkref_gtf_dict_stepsize", as_int=True)

    # Create a dictionary for each "bin" of the genome, that maps to a list of genes within or overlapping
    # that bin. The bin size is determined by gtf_dict_stepsize.
    starts_rounded = gene_starts.apply(
        lambda s: np.floor(s / gtf_dict_stepsize) * gtf_dict_stepsize
    ).values
    ends_rounded = gene_ends.apply(
        lambda s: np.ceil(s / gtf_dict_stepsize) * gtf_dict_stepsize
    ).values
    gene_ids = gene_starts.index
    start_dict = gene_starts.to_dict()
    end_dict = gene_ends.to_dict()
    gene_dict = defaultdict(list)
    for i in range(len(gene_starts)):
        # 2024-01-06 Pandas future warning says need iloc for int address into series
        # cur_chrom = chroms[i]
        # cur_strand = strands[i]
        cur_chrom = chroms.iloc[i]
        cur_strand = strands.iloc[i]
        cur_start = int(starts_rounded[i])
        cur_end = int(ends_rounded[i])
        cur_gene_id = gene_ids[i]
        for coord in range(cur_start, cur_end + 1, gtf_dict_stepsize):
            if not (cur_gene_id in gene_dict[cur_chrom + ":" + str(coord)]):
                gene_dict[cur_chrom + ":" + str(coord) + ":" + cur_strand].append(
                    cur_gene_id
                )

    # Create a dictionary from genes to exons
    exon_gene_ids = gtf_exon_combined.Attributes.apply(
        lambda s: get_attribute(s, "gene_id")
    ).values

    # Note; Only shift dicts *after* gtf has been written out (above)
    # Coord shift from 1-based GTF coords to 0-based pysam
    exon_starts = gtf_exon_combined.Start.values - 1
    exon_ends = gtf_exon_combined.End.values - 1

    exon_gene_start_end_dict = defaultdict(dict)
    for i in range(len(exon_gene_ids)):
        cur_gene_id = exon_gene_ids[i]
        cur_exon_start = exon_starts[i]
        cur_exon_ends = exon_ends[i]
        exon_gene_start_end_dict[cur_gene_id][cur_exon_start] = cur_exon_ends

    # Add gene_id col then make mappings to that
    gtf_gene_combined["gene_id"] = gtf_gene_combined["Attributes"].apply(
        lambda s: get_attribute(s, "gene_id")
    )
    gene_id_to_gene_names = dict(
        zip(
            gtf_gene_combined["gene_id"].values,
            gtf_gene_combined["Attributes"].apply(
                lambda s: get_attribute(s, "gene_name")
            ),
        )
    )
    gene_id_to_genome = dict(
        zip(
            gtf_gene_combined["gene_id"].values,
            gtf_gene_combined["Chromosome"].apply(lambda s: s.split("_")[0]),
        )
    )
    gene_id_to_strand = dict(
        zip(gtf_gene_combined["gene_id"], gtf_gene_combined["Strand"])
    )
    gene_id_to_chrom = dict(
        zip(gtf_gene_combined["gene_id"], gtf_gene_combined["Chromosome"])
    )
    bt_key = "gene_biotype"
    spipe.report_run_story(f"Collecting biotype info with key='{bt_key}'")
    gene_id_to_biotype = dict(
        zip(
            gtf_gene_combined["gene_id"],
            gtf_gene_combined["Attributes"].apply(lambda s: get_attribute(s, bt_key)),
        )
    )

    # Depending on annotation, 'biotype' may not exist; Set something
    for k, v in gene_id_to_biotype.items():
        if not v:
            gene_id_to_biotype[k] = "NA"

    # Save dictionary with gene info and settings
    gene_info = {
        "gene_bins": gene_dict,
        "genes_to_exons": exon_gene_start_end_dict,
        "gene_starts": start_dict,
        "gene_ends": end_dict,
        "gene_id_to_name": gene_id_to_gene_names,
        "gene_id_to_genome": gene_id_to_genome,
        "gene_id_to_chrom": gene_id_to_chrom,
        "gene_id_to_strand": gene_id_to_strand,
        "gene_id_to_biotype": gene_id_to_biotype,
        "gtf_dict_stepsize": gtf_dict_stepsize,
    }

    spipe.report_run_story2("Saving gene_info file")
    spipe.report_run_story2(f"Number of gene_bins {len(gene_info['gene_bins'])}")

    fname = spipe.filepath("GENO_GENE_INFO", None, genome_dir="output_dir")
    story = genes.write_gene_info(gene_info, fname, verb=False)
    spipe.report_run_story2(f"Write_gene_info: {story}")

    # Save chromosome info
    update_chrom_info(spipe, gtf_exon_combined, gtf_gene_combined)

    # Also save gene_info struct in spipe
    spipe._set_gene_info(gene_info)

    # Return info for json record
    new_dict = {
        "gtf_dict": gtf_dict,
        "gtf_filenames": gtf_filenames,
    }
    return new_dict


def check_gtf_att_data_df(gtf_df, key="gene_id", max_missing=0, max_rows=10):
    """Check gtf attribute data exists for key or list of keys

    max_missing = number of rows that can have missing key
    max_rows = max number of row (indexes) to report if issue

    Return status, story
    """
    assert isinstance(gtf_df, pd.DataFrame), f"Expected DataFrame got {type(gtf_df)}"
    ok = True
    story = ""

    if isinstance(key, list):
        for one_k in key:
            ok, story = check_gtf_att_data_df(
                gtf_df, key=one_k, max_missing=max_missing, max_rows=max_rows
            )
            if not ok:
                break
    else:
        # Missing values have zero len string
        # missing = gtf_df["Attributes"].apply(lambda s: get_one_att_value(s, key)).str.len() < 1
        missing = (
            gtf_df["Attributes"].apply(lambda s: get_attribute(s, key)).str.len() < 1
        )
        if sum(missing) > max_missing:
            ok = False
            story = f"Problem parsing gtf attributes for {sum(missing)} rows (> {max_missing} allowed)\n"
            n_row = 0
            row_list = []
            for index, atts in zip(
                gtf_df[missing].index, gtf_df[missing]["Attributes"].to_list()
            ):
                n_row += 1
                if max_rows and (n_row > max_rows):
                    break
                # Zero based index; report 1-based
                row_list.append(f"Row {index + 1}\t{atts}")
            # One more story line then rows
            if max_rows and (len(missing) > max_rows):
                max_story = f"first {max_rows}"
            else:
                max_story = "all"
            story += f"Listing {max_story} bad data rows\n"
            story += "\n".join(row_list)

    return ok, story


def get_gtf_feature_df(gtf_dict, exon=False):
    """Get single gtf with specific feature type for all species

    gtf_dict = dict with GtfData objs
    exon = flag to annotate as exon, else gene

    Prepends 'Chromosome' name with species name (e.g.'21' >--> 'hg38_21')

    Return gtf dataframe
    """
    feature = "exon" if exon else "gene"

    species = sorted(list(gtf_dict.keys()))
    # Start with gtf dataframe for first species
    new_df = gtf_dict[species[0]].get_df()
    new_df = new_df[new_df["Feature"] == feature].copy()
    # Cast Chromosome to str as it might be int
    new_df["Chromosome"] = species[0] + "_" + new_df["Chromosome"].astype(str)

    # Add any additional species
    for i in range(1, len(species)):
        temp_df = gtf_dict[species[i]].get_df()
        temp_df = temp_df[temp_df["Feature"] == feature].copy()
        temp_df["Chromosome"] = species[i] + "_" + temp_df["Chromosome"].astype(str)
        new_df = pd.concat([new_df, temp_df])
    new_df.index = range(len(new_df))
    return new_df


def write_geno_df(spipe, df, fkey, sep=",", index=False):
    """Write genome-dir dataframe (csv file)

    Return string
    """
    # For mkref, 'output_dir' is 'genome_dir'
    story = spipe.write_df(df, fkey, genome_dir="output_dir", sep=sep, index=index)
    return story


def update_chrom_info(spipe, gtf_exon_combined, gtf_gene_combined):
    """Combine gtf chrom info (given) with fasta (from file)

    Updated file written

    Return string
    """
    f_df = spipe.read_csv("GENO_CHROM_DATA", genome_dir="output_dir")
    # gtf exon col with values = count of records
    e_df = pd.DataFrame(gtf_exon_combined.groupby("Chromosome").size())
    e_df.index.name = "chrom"
    e_df.columns = ["exon"]
    # gtf gene
    g_df = pd.DataFrame(gtf_gene_combined.groupby("Chromosome").size())
    g_df.index.name = "chrom"
    g_df.columns = ["gene"]
    # Combine, fill with zero and save as int
    new_df = pd.concat([f_df, g_df, e_df], axis=1).fillna(0).astype(int)
    story = spipe.write_df(new_df, "GENO_CHROM_DATA", genome_dir="output_dir")
    return story


def fasta_gtf_over_ok(spipe):
    """Check overlap of chrom records between gtf and fasta

    Return boolean
    """
    ok = True
    df = spipe.read_csv("GENO_CHROM_DATA", genome_dir="output_dir")
    # Get explicit genome col from index (as <genome>_<chrom>)
    df["genome"] = df.index.map(lambda x: x.split("_")[0])

    # Add column with count of "unplaced features"
    # This assumes genes and exons are at least 1 base, so fasta seq has to
    #   be at least long enough for all those assigned to each chromosome.
    # If names don't match (e.g. "hg38_1" -vs- "hg38_chr1" this will flag
    #   the names in the gtf that aren't in the fasta. The count of not covered
    #   features is then compared to max allowed to flag issue
    # Limit unplaced to zero in case have to report
    df["unplaced"] = df[["exon", "gene"]].max(axis=1) - df["fasta"]
    df["unplaced"] = df["unplaced"].apply(lambda x: max(0, x))

    # max allowed gtf features not in fasta
    gtf_not_fas = spipe.get_par_val("mkref_max_gtf_not_fas", as_int=True)
    # Uncovered chrom and feature counts
    unplac_ch = len(df[df["unplaced"] > 0])
    unplac_feat = df[df["unplaced"] > 0]["unplaced"].sum()

    # Any unplaced = report or fail
    if unplac_feat:
        # Get str for worst offenders in dataframe and slice extremes
        df.sort_values(by=["unplaced", "fasta"], ascending=[False, True], inplace=True)
        df = pd.concat([df.head(5), df.tail(5)])
        df_story = utils.report_df_str(df, pref="", index_name="", row_num=False)

        if unplac_feat > gtf_not_fas:
            ok = False
            spipe.set_problem(
                "Chromosome feature overlap issue between fasta(s) and gtf(s)"
            )
            spipe.set_problem(
                f"{unplac_feat} features on {unplac_ch} chrom not in fasta; max {gtf_not_fas}"
            )
            spipe.set_problem(df_story)

        else:
            spipe.report_warning(
                "Chromosome feature overlap issue between fasta(s) and gtf(s)"
            )
            spipe.report_warning(
                f"{unplac_feat} features on {unplac_ch} chrom not in fasta; max {gtf_not_fas}"
            )
            spipe.report_warning(df_story)
    else:
        spipe.report_run_story("All gtf genome features overlap with seq chromosomes")

    return ok


def make_gene_tables(spipe):
    """Create and save gene data tables for current genome

    Return status
    """
    if spipe.is_focal_crispr(map_direct=True):
        spipe.report_run_story("Direct mapping fb/crispr guides; No gene table")
        return True

    gene_info = spipe.get_gene_info()
    # Gene data (e.g. what goes with DGE files)
    spipe.report_run_story2("Building gene data table")
    gene_id_to_name = gene_info["gene_id_to_name"]
    gene_id_to_genome = gene_info["gene_id_to_genome"]
    gene_df = pd.DataFrame()
    gene_df["gene_id"] = sorted(list(gene_info["gene_id_to_name"].keys()))
    gene_df["gene_name"] = gene_df["gene_id"].apply(lambda s: gene_id_to_name[s])
    gene_df["genome"] = gene_df["gene_id"].apply(lambda s: gene_id_to_genome[s])

    # Missing names?
    # Convert empty names to nan, then count and maybe fill in with ID
    # 2024-01-26 Pandas future warning says to replace inplace=True is going away
    # gene_df["gene_name"].replace(r"^\s*$", np.nan, regex=True, inplace=True)
    gene_df["gene_name"] = gene_df["gene_name"].replace(r"^\s*$", np.nan, regex=True)

    n_missing = gene_df["gene_name"].isna().sum()
    spipe.report_run_story2(f"Missing gene names: {n_missing}")
    if n_missing and spipe.get_par_val("mkref_id_for_no_name", as_bool=True):
        # gene_df["gene_name"].fillna(gene_df["gene_id"], inplace=True)
        gene_df["gene_name"] = gene_df["gene_name"].fillna(gene_df["gene_id"])

        spipe.report_run_story(f"Replaced {n_missing} missing gene names with IDs")

    write_geno_df(spipe, gene_df, "GENO_GENE_DATA", index=False)

    # Gene stats
    fasta = spipe.filepath("GENO_PROC_FASTA_RAW", None, genome_dir="output_dir")
    spipe.report_run_story2(f"Loading sequences {fasta}")
    fasdic = utils.fasta_name_seq_dict(fasta)
    spipe.report_run_story2("Building stats table")
    stats_df = geno_stats_df(gene_info, fasdic)
    write_geno_df(spipe, stats_df, "GENO_GENE_STATS", index=False)

    return True


def gc_am_perc(seq, percent=2):
    """Calculate GC and ambiguous base percents in sequence

    Return tuple (gc, am) as fractions or percentages
    """
    gc = seq.count("C") + seq.count("G")
    am = seq.count("N")
    denom = len(seq) - am
    if denom == 0:
        gc = -1
    else:
        if percent:
            gc = round(gc * 100 / denom, percent)
            am = round(am * 100 / denom, percent)
        else:
            gc = gc / denom
            am = am / denom
    return (gc, am)


def geno_stats_df(gene_info, fasdic):
    """Generate dataframe with gene stat info per id

    Return dataframe
    """
    rlis = []
    for gene_id in gene_info["gene_id_to_name"].keys():
        name = gene_info["gene_id_to_name"][gene_id]
        ch = gene_info["gene_id_to_chrom"][gene_id]
        st = gene_info["gene_starts"][gene_id]
        en = gene_info["gene_ends"][gene_id]
        strand = gene_info["gene_id_to_strand"][gene_id]
        biotype = gene_info["gene_id_to_biotype"][gene_id]

        # Whole gene seq and GC
        g_seq = fasdic[ch][st:en]
        g_gc, _ = gc_am_perc(g_seq)

        # Exon-only seq and GC; Exons as dict with start = key, end = val
        xdic = gene_info["genes_to_exons"][gene_id]
        ex_starts = sorted(xdic.keys())
        x_seq = ""
        for x_st in ex_starts:
            x_seq += fasdic[ch][x_st - 1 : xdic[x_st]]
        x_gc, _ = gc_am_perc(x_seq)

        rlis.append(
            [
                gene_id,
                name,
                ch,
                st,
                len(g_seq),
                strand,
                len(xdic),
                len(x_seq),
                g_gc,
                x_gc,
                biotype,
            ]
        )

    cols = "gene_id,name,chrom,g_start,g_len,strand,ex_num,ex_len,g_gc,ex_gc,biotype".split(
        ","
    )
    df = pd.DataFrame(rlis, columns=cols)
    df["g_gc"] = df["g_gc"].apply(lambda x: f"{x:0.2f}")
    df["ex_gc"] = df["ex_gc"].apply(lambda x: f"{x:0.2f}")

    return df


def get_STAR_mkref_call(spipe, clean_ss=True):
    """Get command to call STAR for mkref

    Return string
    """
    output_dir = spipe.get_par_val("output_dir", as_path=True)
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    # Get adjustable parameters
    genomeSAindexNbases = spipe.get_par_val(
        "mkref_star_genomeSAindexNbases", as_int=True
    )
    genomeSAsparseD = spipe.get_par_val("mkref_star_genomeSAsparseD", as_int=True)
    lggram = spipe.get_par_val("mkref_star_limitGenomeGenerateRAM", as_int=True)
    # Output; For mkref, 'output_dir' is 'genome_dir'
    in_fasta = spipe.filepath("GENO_PROC_FASTA_RAW", None, genome_dir=output_dir)
    exon_gtf = spipe.filepath("GENO_PROC_EXON_GTF", None, genome_dir=output_dir)

    # Build command line
    xpath = spipe.get_exec("STAR")
    star_command = f"{xpath} \
                    --runMode genomeGenerate \
                    --genomeDir {output_dir} \
                    --genomeFastaFiles {in_fasta} \
                    --runThreadN {nthreads} \
                    --limitGenomeGenerateRAM {lggram} \
                    --genomeSAindexNbases {genomeSAindexNbases} \
                    --genomeSAsparseD {genomeSAsparseD} \
                    "
    # Splicing has additional arg
    if spipe.get_par_val("mkref_star_splicing", as_bool=True):
        star_command += f"--sjdbGTFfile {exon_gtf}"

    # Any extra args?
    exargs = spipe.get_par_val("mkref_star_extra_args", as_str=True)
    if exargs:
        star_command += " "
        star_command += exargs.replace(",", " ")

    # Clean to single space?
    if clean_ss:
        star_command = re.sub(" +", " ", star_command)

    return star_command


def generate_STAR_index(spipe):
    """Call STAR to generate index

    Return status
    """
    command = get_STAR_mkref_call(spipe)
    ok, callout, ex_story = utils.call_exec(command, verb=True)
    # Story to log only (call verb=True already prints outcome)
    spipe.report_run_story(ex_story, to_log=True)

    # Output saved to file
    ver = spipe.get_version()
    fpath = spipe.filepath("GENO_STAR_MKREF_OUT", None, genome_dir="output_dir")
    utils.file_header(ofile=fpath, story="STAR mkref call", ver=ver)
    with open(fpath, "a") as OUTFILE:
        com_str = utils.wrap_exe_com_args(command)
        print("\n# Call:", file=OUTFILE)
        print(com_str, file=OUTFILE)
        print("\n# Output:", file=OUTFILE)
        print(callout, file=OUTFILE)
        spipe.report_run_story(f"STAR output saved to {fpath}")

    if not ok:
        spipe.set_problem("Subprocess failed: " + command)
        spipe.set_problem(ex_story)

    return ok


def write_mkref_env_files(spipe, gtf_info, gfasta_info):
    """Write details of mkref processing to file

    spipe should have gene_info initialized (i.e. as above)
    """
    o_dict = {}
    # Output; For mkref, 'output_dir' is 'genome_dir'
    fkey = "GENO_MKREF_DEF"
    fpath = spipe.filepath(fkey, None, genome_dir="output_dir")
    # Possibly bumping file versions?
    b = utils.bump_logs(fpath)
    story = f"(bumped {b} prior versions)" if b else ""
    db_format = "guide_mm" if spipe.is_focal_crispr(map_direct=True) else "STAR"

    o_dict["header"] = spipe.get_header_dict(desc="mkref process info")
    o_dict["genome_dir"] = spipe.get_genome_dir()
    o_dict["db_format"] = db_format
    o_dict["mkref"] = get_mkref_specs(spipe, gtf_info, gfasta_info)
    o_dict["genome"] = spipe._get_genome_specs()

    # Save file
    utils.write_json(o_dict, fpath)
    spipe.report_run_story(f"Saved env info {fpath} {story}")


def get_mkref_specs(spipe, gtf_info, gfasta_info):
    """Get define settings for mkref (e.g. for json)

    Return dict
    """
    new_dict = {}

    fasta_gtf = True
    if spipe.is_focal_crispr(map_direct=True):
        new_dict["db_format"] = "guide_mm"
        mm_dist = spipe.get_par_val("crsp_match_hamming_dist", as_int=True)
        new_dict["match_hamming_dist"] = f"{mm_dist}"
        fasta_gtf = False

    # fasta_gtf and STAR
    if fasta_gtf:
        new_dict["db_format"] = "STAR"
        new_dict["db_call"] = get_STAR_mkref_call(spipe)
        # Python < 3.9 cannot use this syntax
        # new_dict = new_dict | get_mkref_fasta_gtf_specs(spipe, gtf_info, gfasta_info)
        mkref_dict = get_mkref_fasta_gtf_specs(spipe, gtf_info, gfasta_info)
        new_dict.update(mkref_dict)

    return new_dict


def get_mkref_fasta_gtf_specs(spipe, gtf_info, gfasta_info):
    """Get define settings for fasta/gtf (e.g. STAR) mkref

    Return dict
    """
    new_dict = {}

    # List of things per genome
    genomes = spipe.get_gene_info()["genome_list"]
    new_dict["genome_num"] = len(genomes)
    new_dict["genome_names"] = ",".join(genomes)

    # GTF story
    info_list = []
    for genome, gdata in gtf_info["gtf_dict"].items():
        info_list.append(
            {
                "genome": genome,
                "filename": gdata.get_filename(),
                "flavor": gdata.get_flavor(),
                "num_records": gdata.num_records(),
            }
        )
    new_dict["gtf_data"] = info_list

    # Raw files story
    genome_list = spipe.get_par_val("genome_name", as_list=True)
    fasta_list = spipe.get_par_val("fasta", as_list=True)
    gtf_list = spipe.get_par_val("genes", as_list=True)
    if gfasta_info:
        genome_list = genome_list + gfasta_info["genome"]
        fasta_list = fasta_list + gfasta_info["fasta"]
        gtf_list = gtf_list + gfasta_info["gtf"]
    info_list = []
    for genome, fasta, gtf in zip(genome_list, fasta_list, gtf_list):
        f_info = utils.report_file_info_str(fasta)
        g_info = utils.report_file_info_str(gtf)
        info_list.append({"genome": genome, "fasta": f_info, "gtf": g_info})
    new_dict["input_files"] = info_list

    return new_dict


def clean_up(spipe):
    """Clean up for mkref step"""
    # Files to compress
    # Output; For mkref, 'output_dir' is 'genome_dir'
    gzip_files = []
    gzip_files.append(
        spipe.filepath("GENO_PROC_FASTA_RAW", None, genome_dir="output_dir")
    )
    gzip_files.append(
        spipe.filepath("GENO_PROC_GENE_GTF", None, genome_dir="output_dir")
    )
    gzip_files.append(
        spipe.filepath("GENO_PROC_EXON_GTF", None, genome_dir="output_dir")
    )
    comp_lev = spipe.get_par_val("gzip_compresslevel", as_int=True)
    if not utils.check_and_gzip_file(gzip_files, comp_lev=comp_lev):
        spipe.set_problem(f"Compressing geno files failed {gzip_files}")
        return False

    return True
