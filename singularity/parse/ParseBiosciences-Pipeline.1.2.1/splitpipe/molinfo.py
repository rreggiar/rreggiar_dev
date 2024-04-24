#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Converting bam to read information
# Cell barcode, transcript, gene alignment, and count
#

from multiprocessing import Process
import pandas as pd
from collections import defaultdict
import numpy as np
import pysam


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


# ---------------------------------------------------------------------------
# Columns for transcript assignment csv file
TSCP_CSV_COLS = [
    "bc_wells",
    "genome",
    "gene",
    "gene_name",
    "count",
    "exonic",
    "rt_type",
    "bc_bcis",
    "cell_barcode",
    "polyN",
]
# Column index for 'count' in above list
TSCP_COUNT_COL_IDX = TSCP_CSV_COLS.index("count")


def have_molinfo_ins(spipe, verb=True):
    """Check if Molinfo input files exist

    Return True/False, list of files
    """
    ok = False
    # sorted bam file
    fname = spipe.filepath("PF_BAM_SORTED", None)
    check_files = [fname]
    # If list is good, we have ins
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def have_molinfo_outs(spipe, verb=False):
    """Check if Molinfo output files exist

    Return True/False, list of files
    """
    ok = False
    # Transcript assignment file
    fname = spipe.filepath("PF_TSCP_ASSIGN", None)
    check_files = [fname]
    # If list is good, we have outs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def run_molinfo(spipe):
    """Run molinfo step"""
    spipe.report_proc_step("Molinfo")
    # Pysam generates a warning because bam has no index:
    #   "[E::idx_find_and_load] Could not retrieve index file for <file>"
    # Story: https://github.com/pysam-developers/pysam/issues/939
    # Save and set verbosity zero
    save_pyvb = pysam.set_verbosity(0)

    ok = True
    # Input and output status and list
    i_ok, i_list = have_molinfo_ins(spipe)
    o_ok, o_list = have_molinfo_outs(spipe)
    # Run only if don't have output or fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for molinfo: {bad_list}"
            spipe.set_problem(story)
            ok = False
        else:
            ok = molecule_info(spipe)
    else:
        spipe.report_run_story("Skipping Molinfo and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")

    # Reset pysam verbosty (needed?)
    pysam.set_verbosity(save_pyvb)

    ok = spipe.no_problems()
    spipe.report_proc_step("Molinfo", status=ok)
    return ok


def molecule_info(spipe):
    """Gets molecular info for a bam file. Splits the bamfile into
    nthread chunks and runs in parallel

    Returns status
    """
    # Break bam into parts for threaded processing
    split_bam(spipe)

    Pros = []
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    for i in range(nthreads):
        spipe.report_run_story2(f"Starting thread {i+1}")
        p = Process(target=molecule_info_chunk, args=(spipe, i + 1))
        Pros.append(p)
        p.start()

    for t in Pros:
        t.join()

    # Combine parts from multi-threading
    # transcripts (lines) and reads; Output to stats file
    ok = join_anno_bam_files(spipe)
    if ok:
        n_reads, n_tscp = join_tscp_assignment_files(spipe)
        if n_reads:
            update_pipeline_stats(spipe, n_reads, n_tscp)
        else:
            ok = False
            spipe.set_problem("Combinging tscp_assignment parts failed")

    if ok:
        molinfo_clean_up(spipe)

    return ok


# -------------------------- multithread stuff  ------------------------------


def chunk_filepath(spipe, ibam=False, obam=False, ocsv=False, chunk=None):
    """Get filepath for chunk file

    ibam    = flag for bam input (sorted) file
    obam    = flag for bam output (annotated) file
    ocsv    = flag for csv output (transcript assignment) file
    chunk   = chunk number for multithreading.

    If chunk, filename get this inserted before extension:
        <file>.<ext>    >--> <file>.<chunk>.<ext>

    Return filepath string
    """
    if ocsv:
        ras_fname = spipe.filepath("PF_TSCP_ASSIGN_RAW", None)
        if chunk is None:
            fname = ras_fname
        else:
            dpart, fname = utils.fname_dir_base(ras_fname)
            fname = f"{fname.split('.csv')[0]}.chunk{chunk}.csv"
            fpath = dpart + fname
    elif obam:
        obam_fname = spipe.filepath("PF_BAM_ANNO", None)
        if chunk is None:
            fname = obam_fname
        else:
            dpart, fname = utils.fname_dir_base(obam_fname)
            fname = f"{fname.split('.bam')[0]}.chunk{chunk}.bam"
            fpath = dpart + fname
    elif ibam:
        ibam_fname = spipe.filepath("PF_BAM_SORTED", None)
        if chunk is None:
            fname = ibam_fname
        else:
            dpart, fname = utils.fname_dir_base(ibam_fname)
            fname = f"{fname.split('.bam')[0]}.chunk{chunk}.bam"
            fpath = dpart + fname
    else:
        story = "chunk_filepath called without specifying file type"
        raise ValueError(story)

    return fpath


def split_bam(spipe):
    """Split bam file for multi-thread processing

    Bam expected to be sorted by name (i.e. encoded well indexes, etc)

    Note that partitioned bam chunks often start with many-read cells. This
    is because reads for any cell are not split across different bam chunks
    and chunks are only written (in full) after the *next* different cell
    read has been encountered; Cells with many reads have a higher chance of
    crossing a chunk-size boundary, and these will then get written at the
    start of the next chunk bam file.

    Returns nothing
    """
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    spipe.report_run_story(f"Splitting bam file ({nthreads} threads)")

    bam_fname = spipe.filepath("PF_BAM_SORTED", None)
    samfile = pysam.Samfile(bam_fname)

    # Init (input) bam outputs for each thread
    samfile_chunks = {}
    for i in range(nthreads):
        chunk_fname = chunk_filepath(spipe, ibam=True, chunk=(i + 1))
        samfile_chunks[i] = pysam.Samfile(chunk_fname, "wb", template=samfile)

    stat_dict = spipe.get_seq_pipe_stats_dict()
    aligned_reads = int(
        stat_dict["reads_align_unique"] + stat_dict["reads_align_multimap"]
    )
    # If keeping unmapped reads (i.e. in STAR bam file), total = all inputs
    if spipe.get_par_val("align_star_keep_no_map", as_bool=True):
        total_reads = int(stat_dict["reads_align_input"])
    else:
        total_reads = aligned_reads

    spipe.report_run_story(
        f"Aligned reads = {aligned_reads}, total reads to process = {total_reads}"
    )

    # Number of reads per file. Round up to ensure we don't write (nthreads +1) files.
    reads_per_chunk = int(np.ceil(total_reads / nthreads))
    spipe.report_run_story2(f"Reads per chunk {reads_per_chunk}")

    c = 0
    prev_cell_barcode = ""
    reads = []
    for read in samfile:
        if not read.is_secondary:
            reads.append(read)
            # Encoding of read header line; v0.9.4:
            #   w1_w2_w3__pt__i1_i2_i3__bc1_bc2_bc3__polyN
            #
            # First part w1_w2_w3 is common wells == common cells
            cell_barcode = read.qname.split("__")[0]
            # Only write reads after all reads for a cell have been loaded to avoid
            #   splitting one transcriptome into multiple files:
            if cell_barcode != prev_cell_barcode:
                chunk = int(np.floor(c / reads_per_chunk))
                for r in reads[:-1]:
                    samfile_chunks[chunk].write(r)
                reads = reads[-1:]
            prev_cell_barcode = cell_barcode
            c += 1

    for i in range(nthreads):
        samfile_chunks[i].close()
    samfile.close()


# ------------------------ Main work function (each thread)


def molecule_info_chunk(spipe, chunk=None, fname=None):
    """Gets the molecular info from bamfile.

    Outputs transcript assignment csv file
    Outputs annotated bam file
    """
    min_mapq = spipe.get_par_val("mol_min_mapq_score", as_int=True)

    if chunk is None:
        thread_str = ""
    else:
        thread_str = f" thread {chunk}"

    # In / out file paths
    if fname:
        ibamfile = fname
    else:
        ibamfile = chunk_filepath(spipe, ibam=True, chunk=chunk)
    # print("Input", ibamfile)

    obamfile = chunk_filepath(spipe, obam=True, chunk=chunk)
    ocsvfile = chunk_filepath(spipe, ocsv=True, chunk=chunk)

    # Gene / guide info
    if spipe.is_focal_crispr(map_direct=True):
        guide_dm = True
        guide_info = spipe.get_guide_info()
        gpx_guide_name_to_gene_id = guide_info["gpx_guide_name_to_gene_id"]
        gpx_guide_name_to_guide_name = guide_info["gpx_guide_name_to_guide_name"]
    else:
        guide_dm = False
        gene_info = spipe.get_gene_info()
        gene_dict = gene_info["gene_bins"]
        gtf_dict_stepsize = gene_info["gtf_dict_stepsize"]
        gene_starts = gene_info["gene_starts"]
        gene_ends = gene_info["gene_ends"]
        gene_id_to_name = gene_info["gene_id_to_name"]
        gene_id_to_biotype = gene_info["gene_id_to_biotype"]
        exon_gene_start_end_dict = gene_info["genes_to_exons"]

    # To annotate sample name; dic[well_index] = name
    # Well index may be multi-plate (i.e. "06_23") with one or two barcodes
    wind_to_sample = spipe.get_wind_to_sample_dict()
    multi_plate = True if spipe.num_sample_plates() > 1 else False

    # Collapse similar polyN for the same gene-cell_barcode combination.
    # Write output info for each bc to tscp_assign file
    # Only single header line here; Data added later via "to_csv" in append mode
    with open(ocsvfile, "w") as f:
        print(",".join(TSCP_CSV_COLS), file=f)

    # Get bam input and output
    samfile = pysam.Samfile(ibamfile)
    o_samfile = pysam.Samfile(obamfile, "wb", template=samfile)

    # Collect chromosomes from bam header
    chrom_dict = pd.DataFrame(samfile.header["SQ"]).SN.to_dict()

    bc_wind_list = []
    rt_type_list = []
    bc_bci_list = []
    bc_seq_list = []
    polyN_list = []
    gene_list = []
    exon_asns = []
    cell = 0
    dump = 0
    next_cell = False
    for read in samfile:
        # Header encoded info; v0.9.7:
        #   w1_w2_w3__pt__i1_i2_i3__bc1_bc2_bc3__polyN__dstamp_sindex
        # Split out parts (all strings)
        #   bc_wind = w1_w2_w3 like '02_14_64' = triple of well indexes
        #   rt_type = pt like 'R' 'T'
        #   bc_bcis = i1_i2_i3 like '52_4_64' = triple of bc indexes
        #   bc_seqs = bc1_bc2_bc3 = bc seqs (corrected) from read
        #   polyN = polyN from read
        #   dstamp_sindex = date stamp and sequencing index
        parts = read.query_name.split("__")
        bc_wind, rt_type, bc_bcis, bc_seqs, polyN_seq = parts[:5]
        # Sample name from first well index triple
        wind_parts = bc_wind.split("_")
        mp_wind = wind_parts[0]
        if multi_plate:
            mp_wind += "_" + wind_parts[1]
        samp_name = wind_to_sample[mp_wind]

        # Init gene info in case no gene; For output annotated bam tags
        gene_id = ""
        gene_name = ""
        geno_region = False
        geno_reg_code = "N"

        # Only single, (perfect) quality alignment
        if (not read.is_secondary) and (read.mapping_quality >= min_mapq):
            one_match = False

            # fb/crispr direct map case
            if guide_dm:
                one_match = True
                geno_region = True
                geno_reg_code = "E"
                gene_name = gpx_guide_name_to_guide_name[read.reference_name]
                gene_id = gpx_guide_name_to_gene_id[read.reference_name]

            # Normal genome mapping; gene info, including possible overlap picking
            else:
                gene_match = get_gene(
                    read,
                    chrom_dict,
                    gene_dict,
                    gtf_dict_stepsize,
                    gene_starts,
                    gene_ends,
                    gene_id_to_name,
                    gene_id_to_biotype,
                    exon_gene_start_end_dict,
                )
                # Nothing or ambig
                if len(gene_match) == 1:
                    one_match = True
                    # "E" = exon, "I" = intron
                    geno_region = bool(
                        check_exon_alignment(
                            read, gene_match[0], exon_gene_start_end_dict
                        )
                        > 0.5
                    )
                    geno_reg_code = "E" if geno_region else "I"
                    gene_id = gene_match[0]
                    gene_name = gene_id_to_name[gene_id]

            # Only do something if we've got a gene
            if one_match:
                bc_wind_list.append(bc_wind)
                rt_type_list.append(rt_type)
                bc_bci_list.append(bc_bcis)
                bc_seq_list.append(bc_seqs)
                polyN_list.append(polyN_seq)
                gene_list.append(gene_id)
                exon_asns.append(geno_region)

                # Only want to collapse polyN and write reads to file once, after
                # all the reads from a cell have been loaded. As soon as we see
                # a read from a different cell than the previous read, we process
                # all the reads from the first cell.
                cell += 1
                if cell > 1:
                    # Compare well indexes
                    #  (Not bc seqs, which are uncorrected so may differ by chance)
                    #  (Not bci barcode indexes, which may map to same well==cell )
                    if bc_wind_list[-1] != bc_wind_list[-2]:
                        next_cell = True

                # Process reads from one cell
                if next_cell:
                    output_one_cell_csv(
                        spipe,
                        bc_wind_list,
                        rt_type_list,
                        bc_bci_list,
                        bc_seq_list,
                        polyN_list,
                        gene_list,
                        exon_asns,
                        ocsvfile,
                    )

                    # Reset lists, keeping the first read from the new cell
                    bc_wind_list = [bc_wind]
                    rt_type_list = [rt_type]
                    bc_bci_list = [bc_bcis]
                    bc_seq_list = [bc_seqs]
                    polyN_list = [polyN_seq]
                    gene_list = [gene_id]
                    exon_asns = [geno_region]

                    cell = 1
                    next_cell = False

        # Annotate read and write to output bam
        set_read_tags(
            read,
            gene_id,
            gene_name,
            geno_reg_code,
            bc_wind,
            bc_bcis,
            bc_seqs,
            polyN_seq,
            samp_name,
        )
        o_samfile.write(read)

        dump += 1
        if dump % 100000 == 0:
            spipe.report_run_story3(
                f"Molinfo processed {dump} aligned reads... {thread_str}"
            )

    # Catch any last cell and write this
    if cell > 1:
        output_one_cell_csv(
            spipe,
            bc_wind_list,
            rt_type_list,
            bc_bci_list,
            bc_seq_list,
            polyN_list,
            gene_list,
            exon_asns,
            ocsvfile,
        )

    samfile.close()
    o_samfile.close()


def get_gene(
    read,
    chrom_dict,
    gene_dict,
    gtf_dict_stepsize,
    gene_starts,
    gene_ends,
    gene_id_to_name,
    gene_id_to_biotype,
    exon_gene_start_end_dict,
):
    """Returns a list of genes that read alignment maps to.
    The current implementation is stranded and will not return a gene on the
    opposite strand.
    Input: a read object from pysam Alignment file

    If overlapping multiple genes for this read, pick single 'best' one
    Rules favor coding > 'simple' name > exonic

    Output: a list with one or zero gene_id
    """
    if read.is_reverse:
        strand = "-"
    else:
        strand = "+"

    # Assign each read into a "bin" by rounding the position (<chrom>:<pos>:<strand>) to the nearest
    #   gtf_dict_stepsize
    p = int(np.floor(read.positions[0] / gtf_dict_stepsize) * gtf_dict_stepsize)
    read_pos = chrom_dict[read.tid] + ":" + str(p) + ":" + strand

    # The gene_dict contains all genes that overlap a given "bin". This
    # does not guarantee that the read actually overlaps these genes,
    # but reduces the list of possible genes to a few from >20k.
    potential_genes = gene_dict[read_pos]

    # Now we explicity check to see if the read maps to gene(s) from potential_genes.
    # The current implementation is strict such the full read must be contained in the gene and
    # no part of the read can occur outside the gene.
    matching_genes = []
    for gene_id in potential_genes:
        gene_start, gene_end = gene_starts[gene_id], gene_ends[gene_id]
        if (gene_start <= read.positions[0]) and (read.positions[-1] <= gene_end):
            matching_genes.append(gene_id)

    # Score and pick best if more than one match
    if len(matching_genes) < 2:
        one_gene_list = matching_genes
    else:
        one_gene_list = []
        # List of suspect chars for names
        susp_chars = [".", "-"]
        best_score = -1
        for gene_id in matching_genes:
            gene_name = gene_id_to_name[gene_id]
            gene_biotype = gene_id_to_biotype[gene_id]
            frac_exonic = check_exon_alignment(read, gene_id, exon_gene_start_end_dict)
            # If name has 'supect' chars, set this score 0
            susp_char_score = 1
            for susp_c in susp_chars:
                if susp_c in gene_name:
                    susp_char_score = 0
            # Final score = sum of weighted component scores
            gene_score = 0
            gene_score += 4 * (gene_biotype == "protein_coding")
            gene_score += 2 * susp_char_score
            gene_score += 1 * bool(frac_exonic > 0.5)
            # Best score so far?
            if gene_score > best_score:
                best_score = gene_score
                best_gene = gene_id
        # Save best
        one_gene_list.append(best_gene)
    return one_gene_list


def check_exon_alignment(read, gene_id, exon_gene_start_end_dict):
    """Returns the fraction of bases in a read that align to at least one exonic base"""
    possible_exons = exon_gene_start_end_dict[gene_id]

    # pysam is zero based and gtf is one based
    aligned_positions = np.array(read.positions)  # + 1
    align_array = np.zeros(len(aligned_positions), dtype=bool)

    for start, end in possible_exons.items():
        # print(read.positions)
        # print(start,end)
        # Bitwise or to set boolean region overlap or not
        align_array = align_array | (
            (aligned_positions >= (start)) & (aligned_positions <= (end))
        )
    return np.mean(align_array)


def set_read_tags(
    read,
    gene_id,
    gene_name,
    geno_region,
    bc_wind,
    bc_bcis,
    bc_seqs,
    polyN_seq,
    samp_name,
):
    """Set tags in (pysam) read structure for bam file output
    See: https://samtools.github.io/hts-specs/SAMtags.pdf

    Tag definitions:
        GX:Z: Gene ID
        GN:Z: Gene name
        # OX:Z: UMI                     <-- Official tag; Not used
        pN:Z: PolyN                     <-- Custom tag; Used
        CR:Z: Barcode (uncorrected = seqs)
        CB:Z: Barcode (corrected = 1-based well/cell indexes)
        pB:Z: Barcode index string (1-based, all bc)        <-- Custom tag
        pS:Z: Sample name (from round 1 well indexes)       <-- Custom tag
        RE:A:<E|N|I> Region indicator: Exon, iNtron, Intergenic
    """
    read.set_tag(tag="GX", value=gene_id)
    read.set_tag(tag="GN", value=gene_name)
    read.set_tag(tag="pN", value=polyN_seq)
    read.set_tag(tag="CR", value=bc_seqs)
    read.set_tag(tag="CB", value=bc_wind)
    read.set_tag(tag="pB", value=bc_bcis)
    read.set_tag(tag="pS", value=samp_name)
    read.set_tag(tag="RE", value=geno_region)


def output_one_cell_csv(
    spipe,
    bc_wind_list,
    rt_type_list,
    bc_bci_list,
    bc_seq_list,
    polyN_list,
    gene_list,
    exon_asns,
    ofname,
):
    """Handle one cell worth of transcript assignment csv file

    bc_wind_list = list of well index strings (e.g. '02_14_56') == identifier for this cell
    rt_type_list = list of primer types (e.g. 'R' or 'T' etc)
    bc_bci_list = list of barcode index string (e.g. '9_84_26')
    bc_seq_list = list of (uncorrected) barcode srings (e.g. AAACGATA_AACCGAGA_GTACGCAA)
    polyN = list of seqs
    gene_list = lists of gene ids
    exon_asns = list of exon assignment values

    Appends to output file, ofname
    """
    if spipe.is_focal_crispr(map_direct=True):
        guide_info = spipe.get_guide_info()
        guide_id_to_name = guide_info["gene_id_to_name"]
        genome = guide_info["genome"]
    else:
        gene_info = spipe.get_gene_info()
        gene_id_to_name = gene_info["gene_id_to_name"]
        gene_id_to_genome = gene_info["gene_id_to_genome"]

    # Don't use last item in each list; These belongs to next cell
    df = pd.DataFrame(
        {
            "bc_wells": bc_wind_list[:-1],
            "rt_type": rt_type_list[:-1],
            "bc_bcis": bc_bci_list[:-1],
            "cell_barcode": bc_seq_list[:-1],
            "polyN": polyN_list[:-1],
            "gene": pd.Series(gene_list[:-1]),
            "exonic": exon_asns[:-1],
        }
    )

    # Collapsed data frame combines same + similar gene and polyN together, yielding count
    # In other words, if same gene and polyN, assume one tscp molecule; Count = duplicates
    df_collapsed = collapse_polyN_dataframe(df)
    df_collapsed.columns = ["gene", "bc_bcis", "polyN", "count"]

    # Exon = True if more than half of tscp are classified as such
    exonic_fraction = df.groupby(["gene", "bc_bcis", "polyN"]).exonic.mean()
    d_lis = list(
        zip(
            df_collapsed["gene"].values,
            df_collapsed["bc_bcis"].values,
            df_collapsed["polyN"].values,
        )
    )
    df_collapsed.loc[:, "exonic"] = exonic_fraction.loc[d_lis].values > 0.5

    # Shrink read dataframe to tscp dataframe
    # Group and sort as collapsed, then reindex to align values
    df = df.drop_duplicates(["gene", "bc_bcis", "polyN"]).sort_values(
        ["gene", "bc_bcis", "polyN"]
    )
    df.reset_index(drop=True, inplace=True)
    # Copy over
    df_collapsed["cell_barcode"] = df["cell_barcode"]
    df_collapsed["bc_wells"] = df["bc_wells"]
    df_collapsed["rt_type"] = df["rt_type"]

    if spipe.is_focal_crispr(map_direct=True):
        df_collapsed.loc[:, "gene_name"] = df_collapsed.gene.apply(
            lambda s: guide_id_to_name[s]
        )
        df_collapsed.loc[:, "genome"] = genome
    else:
        # Gene info derived from gene_id
        df_collapsed.loc[:, "gene_name"] = df_collapsed.gene.apply(
            lambda s: gene_id_to_name[s]
        )
        df_collapsed.loc[:, "genome"] = df_collapsed.gene.apply(
            lambda s: gene_id_to_genome[s]
        )

    # Data to file; Column header line written outside loop (so header=False here)
    df_collapsed[TSCP_CSV_COLS].to_csv(ofname, header=False, index=False, mode="a")


def collapse_polyN_dataframe(df):
    """Collapses similar (<= 1 hamming dist) polyN seqs with the same cell_barcode-gene combination

    Return dataframe
    """
    # If records share gene, barcode indexes, (identical) polyN, collapse these and get count
    counts = df.groupby(["gene", "bc_bcis", "polyN"]).size().reset_index()
    counts.columns = ["gene", "bc_bcis", "polyN", "count"]
    # Now see if groups of polyN may be further combined (including mismatch)
    counts_df = counts.groupby(["gene", "bc_bcis"]).apply(
        lambda x: pd.Series(collapse_polyN(list(x["polyN"]), list(x["count"])))
    )

    # Pandas returns different content depending on groupby results
    # Need to do this in case counts_df only has one
    if not isinstance(counts_df, pd.Series):
        polyN = counts_df.columns[0]
        polyN_counts = counts_df.iloc[:, 0]
        counts_df["polyN"] = polyN
        counts_df["count"] = polyN_counts
        counts_df = counts_df[["polyN", "count"]]

    counts_df = counts_df.reset_index()
    return counts_df


bases = list("ACGT")


def single_mut_seqs(seq, n):
    """Return a list of all sequences with a 1 nt mutation (no ins/del)."""
    mut_seqs = []
    for i in range(n):
        for b in bases:
            if b != seq[i]:
                mut_seqs.append(seq[:i] + b + seq[i + 1 :])
    return mut_seqs


def hamdist2(str1, str2):
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


def collapse_polyN(kmer_list, counts):
    """Collapses polyN that are within 1 hamming dist of each other.
    Collapses by adding together counts of both similar seqs.

    Returns dict with d[polyN] = count
    """
    kmer_len = len(kmer_list[0])
    kmer_seqs_dict = dict(zip(kmer_list, np.zeros(len(kmer_list))))
    ham_dict = defaultdict(list)
    seq_len = len(kmer_list)
    # Collect lists of within-Hamming-dist seqs for each seq
    # If there are many seqs, create a dict with all one-offs
    #    Total one-off for 8-mer = 24 (8 x 3)
    if seq_len > 50:
        for seq in kmer_list:
            mut_seqs = single_mut_seqs(seq, kmer_len)
            for ms in mut_seqs:
                try:
                    kmer_seqs_dict[ms]
                except Exception:
                    pass
                else:
                    ham_dict[ms].append(seq)
    # Else direct seq-to-seq comparison
    # (if only single seq, seq_len==1, this does nothing)
    else:
        for i in range(seq_len):
            for j in range(i + 1, seq_len):
                s1, s2 = kmer_list[i], kmer_list[j]
                ham_dist = hamdist2(s1, s2)
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
                    break
            except KeyError:
                pass

    return kmer_counts


def join_anno_bam_files(spipe):
    """Join annotated bam file chunks from multi-threading into single file

    Return status
    """
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    spipe.report_run_story2(f"Combining {nthreads} annotated bam parts")

    # Write out all chunk filenames to temp file and collect name list
    t_fname = spipe.filepath("TMP_TEMP1_FILE", None)
    tfile = open(t_fname, "w")
    chunk_names = []
    for i in range(nthreads):
        fname = chunk_filepath(spipe, obam=True, chunk=(i + 1))
        chunk_names.append(fname)
        print(fname, file=tfile)
        spipe.report_run_story3(f"bam part file {fname}")
    tfile.close()

    # Concat parts together with samtools
    xpath = spipe.get_exec("samtools")
    bam_fname = spipe.filepath("PF_BAM_ANNO", None)
    command = f"{xpath} cat -b {t_fname} -o {bam_fname}"

    # Any extra args?
    exargs = spipe.get_par_val("mol_samtools_extra_args", as_str=True)
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

    # Clean up; Remove temp file and chunk files
    if not spipe.keep_temp_files():
        chunk_names.append(t_fname)
        utils.check_and_rm_file(chunk_names)
    else:
        story = "Not cleaning mol(info) temp anno bam files"
        spipe.report_run_story(story)

    return True


def join_tscp_assignment_files(spipe):
    """Combine (multi-thread) tscp assignment parts into single file

    Return tuple (total number of reads, transcripts)
    """
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    spipe.report_run_story2(f"Combining {nthreads} transcript assignment parts")

    # Collect filename list
    filenames = []
    for i in range(nthreads):
        filenames.append(chunk_filepath(spipe, ocsv=True, chunk=(i + 1)))

    # New combo output file; 'RAW' isn't compressed
    oname = spipe.filepath("PF_TSCP_ASSIGN_RAW", None)
    nlines = 0
    ncounts = 0
    with open(oname, "w") as outfile:
        # Write header
        print(",".join(TSCP_CSV_COLS), file=outfile)
        for fname in filenames:
            with open(fname) as infile:
                # Skip (don't copy) header line from each file:
                infile.readline()
                for line in infile:
                    outfile.write(line)
                    nlines += 1
                    ncounts += int(line.split(",")[TSCP_COUNT_COL_IDX])
    return ncounts, nlines


def update_pipeline_stats(spipe, n_reads, n_tscp):
    """Save stats; Append so no header"""
    # Update stats file; Load data frame, update rows, write back out
    df = spipe.read_csv("PF_STAT_PIPE", names="statistic,value", verb=False)
    df.loc["reads_map_transcriptome", "value"] = n_reads
    df.loc["number_of_tscp", "value"] = n_tscp
    # Dump all lines to log; n_rows=0
    spipe.write_df(df, "PF_STAT_PIPE", n_rows=0)


def molinfo_clean_up(spipe):
    """Clean up molinfo step files"""
    # print(">> molinfo_clean_up")
    # Remove temporary split files
    if not spipe.keep_temp_files():
        nthreads = spipe.get_par_val("nthreads", as_int=True)
        spipe.report_run_story("Cleaning molinfo temp files")
        for i in range(nthreads):
            chunk_csv = chunk_filepath(spipe, ocsv=True, chunk=(i + 1))
            chunk_bam = chunk_filepath(spipe, ibam=True, chunk=(i + 1))
            # print(f"++ removing '{chunk_csv}' '{chunk_bam}'")
            utils.check_and_rm_file(chunk_csv)
            utils.check_and_rm_file(chunk_bam)
        # Non-annotated bam file
        fname = spipe.filepath("PF_BAM_SORTED", None)
        utils.check_and_rm_file(fname)
    else:
        spipe.report_run_story("Not cleaning mol(info) temp csv, bam files")

    # Compress tscp assignments
    spipe.report_run_story("Compressing transcript assignments csv")
    fname = spipe.filepath("PF_TSCP_ASSIGN_RAW", None)
    utils.check_and_gzip_file(fname)
    # print("<< molinfo_clean_up")
