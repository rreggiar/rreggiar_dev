#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    # import utils
    pass
else:
    # from splitpipe import utils
    pass


def explain_story(version):
    story = f"""
--------------------------------------------------------------------------------
Parse Biosciences data analysis pipeline.
{version}
--------------------------------------------------------------------------------
Modes

    --------------------
    mkref   Make reference genome; Need genome fasta and gtf file(s), name(s).
            Also need to specify output directory for processed genome files.

    --------------------
    all     >> All below modes, run in order, for 1 sublibrary <<

    inqc    Input QC (check fastq inputs)
    pre     Pre-process fastq files (barcode correction + cDNA)
    align   Align barcode-corrected reads (STAR).
    post    Post-process read alignment (sort bam, samtools).
    mol     Molecule info extraction (transcripts and genes from reads).
    [dge, TCR...]     
            Count matrix creation and stats. Specific step depends on use case. 
    ana     Analysis of cells and samples (clustering, etc).
    rep     Generate output reports for each sample.
    clean   Clean up intermediate files (zip / compress / delete).

    --------------------
    comb    Combine outputs from multiple (>= 2) sublibraries.

            Inputs must be processed up to at least 'dge' step; Need to supply
                list of 'output_dir' paths via '--sublibraries X/ Y/ Z/'.
            This step is ~ equivalent to the 'dge' step but new files / stats are
                derived by combining previously calculated values.
            Samples, chemistry, kit, genome, etc are taken from sublib run info;
                These parameters are not needed (and are ignored).

    To run a subset of steps, start at <mode1> and use '--until_step' <mode2>
        to specify the end mode. See list of full-run modes with --dryrun.

--------------------------------------------------------------------------------
Specifying chemistry and kit

    The chemistry version for the data must be specified (via --chemistry)

    Kit is normally determined automatically from a subsample of data.
    However, this may be specified (via --kit) if needed.

    To see valid chemistry and kit names, use --kit_list

    To bypass kit score checking, use --kit_score_skip

--------------------------------------------------------------------------------
Specifying sample well format

    --sample <name> <wells>

    Samples are specified by <name> and <wells>

    Wells are specified in blocks, ranges, or individually like this:
        'A1:C6' specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6.
        'A1-B6' specifies a range as [start]-[end]; A1-A12, B1-6.
        'C4' specifies a single well.
        Multiple selections are joined by commas (NO SPACES), e.g. 'A1-A6,B1:D3,C4'

    Sample specification may be supplied via --sample, --samp_list, or --samp_sltab.

    Sample names must be unique and cannot be ambiguous with respect to wells; If a
        sample is named as a well, it is prepended like 'd4' >--> 'sample_d4'.

--------------------------------------------------------------------------------
TCR single sublibrary analysis or combining TCR outputs

    The --trc_analysis option is needed for single sublib analysis (i.e. --mode all)

    The script install_TCR.sh must be run (and pipeline reinstalled) for TCR analysis.

    Use --tcr_check to check on installed TCR databases"

--------------------------------------------------------------------------------
CRISPR single sublibrary analysis or combining CRISPR outputs

    The --crispr option also requires specification of guides (--crsp_guides)
        and a parent output path (--parent_dir)

--------------------------------------------------------------------------------
Output files

    Html reports are named as <sample_name>_summary_analysis.html. Aggregated
        metrics for all samples are in the agg_samp_ana_summary.csv file.

    Digital gene expression data is in <sample_name>/DGE_* subdirs. Per cell
        counts are in 'count_matrix.mtx' sparse matrix files. Cell and gene data
        are in accompanying '.csv' files. Matrix rows are cells, columns are genes.


    Other sample-specific data files are in <sample_name>/<report>.

    Global (all sample) data files are the process subdir.

    Parameter options determine which working files are kept, and new files saved.
"""
    print(story)
