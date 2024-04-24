#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

from collections import OrderedDict, namedtuple

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import kits
    import utils
else:
    from splitpipe import kits
    from splitpipe import utils


# ----------------- Default parameter values ---------------------------------


def get_spipe_default_pars(par_dict=None):
    """Get dict with default pipeline parameters

    NOTE: If changing default parameters, update parameter file 'spipe.par' too

    Return dict
    """
    if not par_dict:
        par_dict = OrderedDict()

    # Keys used below can have '.' and UpperCase; Cleaned prior to use
    # Run mode and kit parameters empty
    par_dict["mode"] = ""
    par_dict["chemistry"] = ""
    par_dict["kit"] = ""
    par_dict["kit_source"] = ""
    par_dict["use_case"] = ""
    par_dict["run_name"] = ""
    # paths / files
    par_dict["parfile"] = None
    par_dict["input_dir"] = None
    par_dict["output_dir"] = None
    par_dict["parent_dir"] = None  # Path to output dir
    par_dict["fq1"] = None
    par_dict["fq2"] = None
    par_dict["sublib_list"] = None  # File listing sublibrary paths
    par_dict["sublib_pref"] = ""  # Sublibrary list path prefix
    par_dict["sublib_suff"] = ""  # Sublibrary list path suffix
    par_dict["sublibraries"] = []  # list of <path>
    par_dict["samp_list"] = None  # File listing sample defs, <name> <wells>
    par_dict["samp_sltab"] = None  # SampleLoadingTable sample defs (spreadsheet)
    par_dict["sample"] = []  # list of <name> <wells>
    par_dict["genome_dir"] = None

    # Mkref genome stuff
    par_dict["genome_name"] = []  # Names for genomes
    par_dict["genes"] = []  # Gene definitions; gtf files
    par_dict["fasta"] = []  # Genome sequence; fasta files
    par_dict["gfasta"] = []  # Combined genome name and seq (e.g. for focal / crispr)

    # Target enrichment / focal barcoding
    par_dict["targeted_list"] = None  # List of target genes

    # Run-time settings
    par_dict["nthreads"] = 0
    par_dict["max_threads"] = 32
    par_dict["rseed"] = 42
    par_dict["reuse"] = False
    par_dict["make_allwell"] = False
    par_dict["make_allsample"] = True
    par_dict["samp_wells_unique"] = True
    par_dict["keep_going"] = True
    par_dict["keep_temps"] = False
    par_dict["no_compress"] = False
    par_dict["gzip_compresslevel"] = 4
    par_dict["log_memory_use"] = True
    par_dict["one_step"] = False
    par_dict["until_step"] = ""
    par_dict["start_timeout"] = 30
    par_dict["clear_runproc"] = False
    par_dict["dryrun"] = False
    par_dict["kit_list"] = False  # Flag to list kits, chem
    par_dict["chem_list"] = False  # Flag to list kits, chem
    par_dict["bc_list"] = False  # Flag to list barcode sets

    # Barcode
    par_dict["bc_save_data_clean"] = True
    par_dict["bc_amp_seq"] = ""  # Chem-specific default set later if needed
    par_dict["bc_round_set"] = []  # Kit-specific default set later if needed
    par_dict["bc_edit_dist"] = 2
    # Kit barcode score
    par_dict["kit_score_warn"] = 0.7
    par_dict["kit_score_min"] = 0.1
    par_dict["kit_score_remain_frac"] = 0.2
    par_dict["kit_score_skip"] = False

    # TSO primer
    par_dict["tso_seq"] = "AACGCAGAGTGAATGGG"
    par_dict["tso_read_limit"] = 15

    # Fastq read sampling (stats, bc correct, kit guess)
    par_dict["fastq_samp_slice"] = [100000, 1100000]

    # Parameters for different pipeline steps
    par_dict["mkref_id_for_no_name"] = True
    par_dict["mkref_gtf_guess_min"] = 1
    par_dict["mkref_min_genes"] = 1
    par_dict["mkref_min_exons"] = 1
    par_dict["mkref_max_gtf_not_fas"] = 0
    par_dict["mkref_gtf_dict_stepsize"] = 10000
    par_dict["mkref_star_splicing"] = True
    par_dict["mkref_star_genomeSAindexNbases"] = 14
    par_dict["mkref_star_genomeSAsparseD"] = 2
    par_dict["mkref_star_limitGenomeGenerateRAM"] = 24000000000
    par_dict["mkref_star_extra_args"] = None

    par_dict["inqc_fastqc_run"] = True
    par_dict["inqc_fastqc_show_warn"] = False
    par_dict["inqc_min_reads"] = 100000
    par_dict["inqc_min_good_R1len_frac"] = 0.5
    par_dict["inqc_min_good_R2len_frac"] = 0.9
    par_dict["inqc_min_perf_frac"] = 0.3
    par_dict["inqc_warn_perf_frac"] = 0.5

    par_dict["pre_min_R1_len"] = 45
    par_dict["pre_check_name_match"] = True
    par_dict["pre_check_name_R1R2"] = True
    par_dict["pre_head_keep_seq_index"] = True
    # Regex to capture read <name> <read1|2> <seq_index>
    # Expecting something like this:
    #   @A01366:302:HGFYTDMXY:1:1101:1090:1000 1:N:0:CAGATCAT
    par_dict["pre_head_regex"] = r"@(\S+)\s+([12])\S*:(\S*)$"
    par_dict["pre_perfect_bc_cell_thresh"] = 0.92
    par_dict["pre_keep_no_bc_fq"] = False

    # TCR specific
    par_dict["tcr_analysis"] = False  # Use case flag
    par_dict["tcr_genome"] = "human"
    par_dict["tcr_use_imgt_db"] = False
    par_dict["tcr_kmer"] = 21
    par_dict["tcr_read_threshold"] = 2
    par_dict["tcr_trim"] = 30

    # Focal/CRISPR specific
    par_dict["focal_barcoding"] = False
    par_dict["crsp_analysis"] = False
    par_dict["crsp_guides"] = ""
    par_dict["crsp_use_star"] = False
    par_dict["crsp_guide_offset"] = 0
    par_dict["crsp_pref_samp_len"] = 10
    par_dict["crsp_match_hamming_dist"] = 1

    # Normal case align
    par_dict["align_min_genome_version"] = "v1.1.3"
    par_dict["align_min_start_reads"] = 50000
    par_dict["align_star_outFilterMultimapNmax"] = 3
    par_dict["align_star_keep_no_map"] = False
    par_dict["align_star_extra_args"] = None
    # Post-(alignment)-processing
    par_dict["post_warn_map_frac"] = 0.33
    par_dict["post_min_map_frac"] = 0.1
    par_dict["post_samtools_extra_args"] = None
    # Molinfo (bam to tscp assignment)
    par_dict["mol_min_mapq_score"] = 255
    par_dict["mol_samtools_extra_args"] = None

    # DGE step
    par_dict["dge_min_start_tscp"] = 10000
    par_dict["dge_copy_dge_genes_file"] = True
    par_dict["dge_save_express_gene_file"] = True
    par_dict["dge_ssc_xstep_size"] = 1000
    par_dict["dge_ssc_xstep_ndup"] = 5
    par_dict["dge_species_purity_thresh"] = 0.9
    # Filtered DGE; transcript cutoff controls
    par_dict["dge_tscp_cut_given"] = None
    par_dict["dge_tscp_cut_min"] = 30
    par_dict["dge_tscp_cut_max"] = None
    par_dict["dge_tscp_cut_scale_fac"] = 1.0
    # Filtered DGE; cell def controls
    par_dict["dge_cell_est"] = None
    par_dict["dge_cell_given"] = None
    par_dict["dge_cell_gxfold"] = 3
    par_dict["dge_cell_min"] = 10
    par_dict["dge_cell_auc_min"] = 0.4
    par_dict["dge_cell_max"] = 100000  # NOTE; Kit-specific default set later
    par_dict["dge_cell_scale_by_frac"] = True
    par_dict["dge_cell_hard_min"] = 5
    par_dict["dge_cell_list"] = None  # List of barcodes to use
    # Filtered DGE; genes
    par_dict["dge_gene_cell_min"] = 3
    # Unfiltered
    par_dict["dge_tscp_cut_umin"] = 10
    par_dict["dge_gene_cell_umin"] = 1
    # Target enrichment
    par_dict["dge_tscp_cut_targeted_min"] = 15
    par_dict["dge_tscp_cut_targeted_umin"] = 5
    par_dict["dge_gene_cell_targeted_min"] = 2
    # Transcript min read / tscp support; General and focal barcoding / crispr
    par_dict["dge_tscp_read_min"] = 1
    par_dict["dge_gene_tscp_min"] = 1
    # par_dict["dge_tscp_read_fb_min"] = 3
    # par_dict["dge_gene_tscp_fb_min"] = 5
    par_dict["dge_tscp_read_fb_min"] = 1
    par_dict["dge_gene_tscp_fb_min"] = 1
    par_dict["dge_gene_tscp_fb_umin"] = 1  # Used as a default filter in spclass
    par_dict["dge_tscp_read_fb_umin"] = 1  # Used as a default filter in spclass

    # Combine mode
    par_dict["comb_new_tscp"] = False
    par_dict["comb_unfilt_dge"] = True
    par_dict["comb_max_sublibs_mega"] = 16
    par_dict["comb_max_sublibs_wt"] = 8
    par_dict["comb_max_sublibs_mini"] = 2

    # Analysis
    par_dict["ana_scanpy_verbosity"] = ""
    par_dict["ana_save_anndata"] = True
    par_dict["ana_mem_eff_cell_thresh"] = 100000
    par_dict["ana_norm_target_sum"] = 10000
    par_dict["ana_scale_max_value"] = 10
    par_dict["ana_scale_zero_center"] = True
    par_dict["ana_pca_min_hvar"] = 100
    par_dict["ana_pca_dim_max"] = 50
    par_dict["ana_pca_dim_min"] = 10
    par_dict["ana_pca_dim_targeted_max"] = 10
    par_dict["ana_pca_irlbpy"] = True
    par_dict["ana_cluster_min_bc"] = 30
    par_dict["ana_cluster_1base"] = True
    par_dict["ana_save_umap"] = True
    par_dict["ana_umap_n_neighbors"] = 50
    par_dict["ana_umap_min_dist"] = 0.5
    par_dict["ana_cluster_alg"] = "leiden"
    par_dict["ana_ctab_top_n"] = 20
    par_dict["ana_ctab_sco_min_den"] = 0.01
    par_dict["ana_gene_cell_min"] = 3

    # Report output stuff
    par_dict["rep_save_figs"] = True
    par_dict["rep_plate_min_zero"] = True
    par_dict["rep_umap_subsamp"] = 15000
    par_dict["rep_umap_ss_clus_min"] = 10
    par_dict["rep_umap_norm_fac"] = 1e4
    par_dict["rep_umap_fb_zero_nc"] = True
    # Urls to get report on specific gene (by id)
    # NCBI is not organism specific
    ncbi = "https://www.ncbi.nlm.nih.gov/gene/?term={gene_id}"
    par_dict["rep_gene_def_url"] = ncbi
    # embl = "https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={gene_id}"
    # par_dict['rep_gene_def_url'] = embl
    par_dict["rep_gene_url"] = []  # list of <genome> <url>
    par_dict["rep_alt_sty_fmt"] = False

    # Clean step stuff
    par_dict["clean_per_step"] = False
    par_dict["clean_make_summaries_zip"] = True

    # Map target stuff
    par_dict["map_targ_genome_dir_path"] = "map_targs"
    par_dict["map_targ"] = []  # list of <name> <fasta>
    par_dict["map_targ_kmer_len"] = 30
    par_dict["map_targ_kmer_step"] = 1
    par_dict["map_targ_kmer_search_max"] = 41
    par_dict["map_targ_kmer_search_step"] = 10

    # Standardized keys
    par_dict = clean_par_dict_keys(par_dict)

    return par_dict


def clean_param_key(key):
    """Return cleaned up version of parameter key"""
    return key.lower().replace(".", "_")


def param_num_vals(key):
    """Return the number of values for parameter key item

    For example
        'sample' = <name> <wells = 2
        'output_dir' = <path> = 1

    Return int
    """
    twos = [
        "sample",
        "map_targ",
        "rep_gene_url",
        "bc_round_set",
        "gfasta",
        "fastq_samp_slice",
    ]

    key = clean_param_key(key)
    if key in twos:
        num = 2
    else:
        num = 1
    return num


def get_spipe_temp_par_dict():
    """Get dict with command-line-only (temp) parameters mapped to full length

    Return dict of tuples
    """
    TEMP_PAR_MAPPING = """
        # Short command line and internal parameter versions with possibile logic
        tscp_use            dge_tscp_cut_given
        tscp_min            dge_tscp_cut_min
        tscp_max            dge_tscp_cut_max
        cell_use            dge_cell_given
        cell_est            dge_cell_est
        cell_xf             dge_cell_gxfold
        cell_min            dge_cell_min
        cell_max            dge_cell_max
        cell_list           dge_cell_list
        save_anndata        ana_save_anndata
        no_save_anndata     ana_save_anndata        not
        save_figs           rep_save_figs
        no_keep_going       keep_going              not
        yes_allwell         make_allwell
        no_allwell          make_allwell            not
        yes_allsample       make_allsample
        no_allsample        make_allsample          not
        crispr              crsp_analysis
        crsp_read_thresh    dge_tscp_read_fb_min
        crsp_tscp_thresh    dge_gene_tscp_fb_min
        crsp_max_mm         rsp_match_hamming_dist
        use_imgt_db         tcr_use_imgt_db

        # These are command line only; There is no long parameter
        clear_runproc       none
        kit_list            none
        chem_list           none
        bc_list             none
        tcr_check           none
        explain             none
    """
    par_dict = {}
    for line in TEMP_PAR_MAPPING.split("\n"):
        # Split off anything after comment, then split on words
        parts = line.split("#")[0].split()
        if len(parts) < 2:
            continue

        temp_par = parts[0]
        long_par = parts[1]
        if long_par.lower() == "none":
            long_par = ""
        logic = "" if len(parts) < 3 else parts[2]

        par_dict[temp_par] = (long_par, logic)

    return par_dict


def clean_par_dict_keys(par_dict):
    """Clean up parameter keys to lower case, underscore

    Return dict
    """
    # If not a dict, assume argparse.ArgumentParser and get dict
    if not isinstance(par_dict, dict):
        if par_dict is None:
            par_dict = {}
        else:
            par_dict = vars(par_dict)

    # Clean up keys and copy values
    new_pars = {}
    for k, v in par_dict.items():
        new_k = k.lower().replace(".", "_")
        new_pars[new_k] = v
    return new_pars


def get_kit_default_pars(kit, chem):
    """Get kit-specific default parameters

    Return dict
    """
    par_dict = {}
    # Parse first
    chem = kits.parse_chemistry(chem)
    kit_s, _ = kits.parse_kit(kit, chem=chem)

    # Barcodes for given kit and chem
    bc_rounds = []
    if not kits.is_custom_kit(kit_s):
        bc_rounds = kits.kit_bc_set_list(kit_s, chem, as_tup=True)
    par_dict["bc_round_set"] = bc_rounds

    # Max cells for mega is 150K, else 30K
    cell_max = 30000
    if kit == "WT_mega":
        cell_max = 150000
    par_dict["dge_cell_max"] = cell_max

    return par_dict


def get_chem_default_pars(chem):
    """Get chem-specific default parameters

    Return dict
    """
    par_dict = {}

    chem = kits.parse_chemistry(chem)
    # Read structure
    # v3 chemistry
    if kits.parse_chemistry(chem, as_num=True) == 3:
        seq = "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111"
    else:
        seq = "NNNNNNNNNN33333333GTGGCCGATGTTTCGCATCGGCGTACGACT22222222ATCCACGTGCTTGAGACTGTGG11111111"
    par_dict["bc_amp_seq"] = seq

    return par_dict


def get_crispr_default_pars():
    """Get default parameters for crispr use case

    Return dict
    """
    PARAM_STR = """
        mkref_star_genomeSAindexNbases 3
        pre_min_R1_len      10
        dge_min_start_tscp  5000
        dge_tscp_cut_min    1
        dge_tscp_cut_umin   1
        dge_gene_cell_min   1
        dge_gene_cell_umin  1
        dge_cell_min        1
        dge_cell_hard_min   1
    """
    par_dict = {}
    for line in PARAM_STR.split("\n"):
        # Split off anything after comment, then split on words
        parts = line.split("#")[0].split()
        if len(parts) == 2:
            par_dict[parts[0]] = parts[1]
    return par_dict


# ----------------- Dir and file names --------------------------------------
# These formats are used to get actual paths via 'filepath(key,sample_name)'
FILE_DEFS = """
    # Output dirs
    DIR_TOP             {top_dir}
    DIR_PROC            {top_dir}/process
    DIR_FASTQC          {top_dir}/process/fastQC
    DIR_SAMP            {top_dir}/{sample_name}
    DIR_S_REP           {top_dir}/{sample_name}/report
    DIR_S_FIG           {top_dir}/{sample_name}/figures
    DIR_S_DGE_N_F       {top_dir}/{sample_name}/DGE_filtered
    DIR_S_DGE_N_UF      {top_dir}/{sample_name}/DGE_unfiltered
    DIR_S_DGE_T_F       {top_dir}/{sample_name}/DGE_targeted_filtered
    DIR_S_DGE_T_UF      {top_dir}/{sample_name}/DGE_targeted_unfiltered
    DIR_S_DGE_FB_F      {top_dir}/{sample_name}/focal_filtered
    DIR_S_DGE_FB_UF     {top_dir}/{sample_name}/focal_unfiltered
    DIR_S_DGE_CRSP_F    {top_dir}/{sample_name}/guide_RNAs_filtered
    DIR_S_DGE_CRSP_UF   {top_dir}/{sample_name}/guide_RNAs_unfiltered
    DIR_S_TCR_F         {top_dir}/{sample_name}/TCR_filtered
    DIR_S_TCR_UF        {top_dir}/{sample_name}/TCR_unfiltered
    DIR_SCRATCH         {top_dir}/process/scratch/
    DIR_PROC_FB_GENO    {top_dir}/process/scratch/focal_genome/
    DIR_PROC_TCR_F      {top_dir}/process/scratch/TCR_filtered/
    DIR_PROC_TCR_UF     {top_dir}/process/scratch/TCR_unfiltered/
    DIR_PROC_CRSP       {top_dir}/process/scratch/CRISPR/

    # Genome dir files
    DIR_GENO            {genome_dir}
    GENO_GENE_INFO      {genome_dir}/gene_info.json
    GENO_GUIDE_INFO     {genome_dir}/guide_info.json
    GENO_GENE_DATA      {genome_dir}/all_genes.csv
    GENO_GUIDE_DATA     {genome_dir}/all_guides.csv
    GENO_MKREF_DEF      {genome_dir}/process/mkref_def.json
    GENO_PROC_FASTA_RAW {genome_dir}/process/genome.fas
    GENO_PROC_GENE_GTF  {genome_dir}/process/genes.gtf
    GENO_PROC_EXON_GTF  {genome_dir}/process/exons.gtf
    GENO_GENE_STATS     {genome_dir}/process/all_gene_stats.csv
    GENO_CHROM_DATA     {genome_dir}/process/all_chrom.csv
    GENO_STAR_MKREF_OUT {genome_dir}/process/STAR_mkref_sdtout.txt
    GENO_GFASTA_PREF    {genome_dir}/process/gfasta             # For gfasta-derived

    # Process files
    PF_RUNPROC_DEF      {top_dir}/process/run_proc_def.json
    PF_STAT_PIPE        {top_dir}/process/pre_align_stats.csv
    PF_STAT_SEQ         {top_dir}/process/sequencing_stats.csv
    PF_ALL_GENES        {top_dir}/process/all_genes.csv
    PF_TARGET_GENES     {top_dir}/process/target_genes.csv
    PF_ALL_GUIDES       {top_dir}/process/all_guides.csv
    PF_TSCP_ASSIGN_RAW  {top_dir}/process/tscp_assignment.csv
    PF_TSCP_ASSIGN      {top_dir}/process/tscp_assignment.csv.gz
    PF_FASTQ_BC_RAW     {top_dir}/process/barcode_head.fastq
    PF_FASTQ_BC_PREF    {top_dir}/process/barcode_head          # Prefix (align STAR)
    PF_FASTQ_BC         {top_dir}/process/barcode_head.fastq.gz
    PF_FASTQ_FIN_OUT    {top_dir}/process/barcode_headLog.final.out
    PF_FASTQ_SJ_TAB     {top_dir}/process/barcode_headSJ.out.tab
    PF_BAM_RAW          {top_dir}/process/barcode_headAligned.out.bam
    PF_BAM_SORTED       {top_dir}/process/barcode_headAligned_sorted.bam
    PF_BAM_SORTED_PREF  {top_dir}/process/barcode_headAligned_sorted    # Prefix (samtools); not filename
    PF_BAM_ANNO         {top_dir}/process/barcode_headAligned_anno.bam
    PF_BAM_NOMAP        {top_dir}/process/barcode_headAligned_nomap.bam
    PF_FASTQ_NOBC_R1    {top_dir}/process/no_barcode_R1.fastq
    PF_FASTQ_NOBC_R2    {top_dir}/process/no_barcode_R2.fastq
    PF_BARCODE_DATA     {top_dir}/process/barcode_data.csv
    PF_VAL_BC_DATA      {top_dir}/process/well_vbc_counts.csv
    PF_VAL_2BC_DATA     {top_dir}/process/well_2vbc_counts.csv
    PF_MAP_TARG_CT      {top_dir}/process/map_targ_bc_counts.csv
    PF_BCC_READ_LENS    {top_dir}/process/bcc_read_length_counts.csv
    PF_KIT_SCORE        {top_dir}/process/kit_bc_scores.csv
    PF_CONDA_INFO       {top_dir}/process/env_conda_info.txt
    PF_VERSION_INFO     {top_dir}/process/env_version_info.txt
    PF_STAR_ALIGN_OUT   {top_dir}/process/STAR_align_sdtout.txt
    PF_FASTQC_SUM_OUT   {top_dir}/process/fastQC_pass_fail_sum.txt
    PF_TCR_SUBPROC_OUT  {top_dir}/process/TCR_subproc_out.txt
    PF_MEM_PROFILE      {top_dir}/process/mem_profile.txt
    PF_CRSP_FASTA       {top_dir}/process/scratch/CRISPR/gRNA.fasta
    PF_CRSP_MKREF_PAR   {top_dir}/process/scratch/CRISPR/mkref.par


    # Temp files
    TMP_SUM_ZIP_DIR     {top_dir}/process/all_sum_temp
    TMP_TEMP1_FILE      {top_dir}/process/temp1.temp
    TMP_TEMP2_FILE      {top_dir}/process/temp2.temp
    TMP_TEMP3_FILE      {top_dir}/process/temp3.temp
    TMP_SF_ASUM_HTML    {top_dir}/process/{sample_name}_tmp_ana_sum.html

    # Top-level, multi-sample
    SF_SUMMARIES_ZIP    {top_dir}/all_summaries.zip
    SF_AGG_ASUM_CSV     {top_dir}/agg_samp_ana_summary.csv

    # Sample files
    SF_ASUM_HTML        {top_dir}/{sample_name}_analysis_summary.html
    # DGE
    SF_NF_MTX           {top_dir}/{sample_name}/DGE_filtered/count_matrix.mtx
    SF_NF_CELL          {top_dir}/{sample_name}/DGE_filtered/cell_metadata.csv
    SF_NF_GENE          {top_dir}/{sample_name}/DGE_filtered/all_genes.csv
    SF_NF_ANA_H5AD      {top_dir}/{sample_name}/DGE_filtered/anndata.h5ad
    # DGE normal unfiltered
    SF_NU_MTX           {top_dir}/{sample_name}/DGE_unfiltered/count_matrix.mtx
    SF_NU_CELL          {top_dir}/{sample_name}/DGE_unfiltered/cell_metadata.csv
    SF_NU_GENE          {top_dir}/{sample_name}/DGE_unfiltered/all_genes.csv
    # DGE targeted filtered
    SF_TF_MTX           {top_dir}/{sample_name}/DGE_targeted_filtered/count_matrix.mtx
    SF_TF_CELL          {top_dir}/{sample_name}/DGE_targeted_filtered/cell_metadata.csv
    SF_TF_GENE          {top_dir}/{sample_name}/DGE_targeted_filtered/target_genes.csv
    SF_TF_ANA_H5AD      {top_dir}/{sample_name}/DGE_targeted_filtered/anndata.h5ad
    # DGE targeted unfiltered
    SF_TU_MTX           {top_dir}/{sample_name}/DGE_targeted_unfiltered/count_matrix.mtx
    SF_TU_CELL          {top_dir}/{sample_name}/DGE_targeted_unfiltered/cell_metadata.csv
    SF_TU_GENE          {top_dir}/{sample_name}/DGE_targeted_unfiltered/target_genes.csv
    # DGE focal barcoding filtered
    SF_FF_MTX           {top_dir}/{sample_name}/focal_filtered/count_matrix.mtx
    SF_FF_CELL          {top_dir}/{sample_name}/focal_filtered/cell_metadata.csv
    SF_FF_GENE          {top_dir}/{sample_name}/focal_filtered/target_genes.csv
    SF_FF_ANA_H5AD      {top_dir}/{sample_name}/focal_filtered/anndata.h5ad
    # DGE focal barcoding unfiltered
    SF_FU_MTX           {top_dir}/{sample_name}/focal_unfiltered/count_matrix.mtx
    SF_FU_CELL          {top_dir}/{sample_name}/focal_unfiltered/cell_metadata.csv
    SF_FU_GENE          {top_dir}/{sample_name}/focal_unfiltered/target_genes.csv
    # CRISPR DGE barcoding filtered
    SF_CF_MTX           {top_dir}/{sample_name}/guide_RNAs_filtered/count_matrix.mtx
    SF_CF_CELL          {top_dir}/{sample_name}/guide_RNAs_filtered/cell_metadata.csv
    SF_CF_GENE          {top_dir}/{sample_name}/guide_RNAs_filtered/all_guides.csv
    SF_CF_ANA_H5AD      {top_dir}/{sample_name}/guide_RNAs_filtered/anndata.h5ad
    # CRISPR DGE barcoding unfiltered
    SF_CU_MTX           {top_dir}/{sample_name}/guide_RNAs_unfiltered/count_matrix.mtx
    SF_CU_CELL          {top_dir}/{sample_name}/guide_RNAs_unfiltered/cell_metadata.csv
    SF_CU_GENE          {top_dir}/{sample_name}/guide_RNAs_unfiltered/all_guides.csv
    # TCR filtered
    SF_TCR_F_STAT       {top_dir}/{sample_name}/report/tcr_filt_metrics.csv
    SF_TCR_F_AIRR       {top_dir}/{sample_name}/TCR_filtered/tcr_annotation_airr.tsv
    SF_TCR_F_CONTIGS    {top_dir}/{sample_name}/TCR_filtered/tcr_contigs.fa
    SF_TCR_F_BC_REPORT  {top_dir}/{sample_name}/TCR_filtered/barcode_report.tsv
    SF_TCR_F_CLON_FREQ  {top_dir}/{sample_name}/TCR_filtered/clonotype_frequency.tsv
    # TCR unfiltered
    SF_TCR_UF_STAT      {top_dir}/{sample_name}/report/tcr_unfilt_metrics.csv
    SF_TCR_UF_AIRR      {top_dir}/{sample_name}/TCR_unfiltered/tcr_annotation_airr.tsv
    SF_TCR_UF_CONTIGS   {top_dir}/{sample_name}/TCR_unfiltered/tcr_contigs.fa
    SF_TCR_UF_BC_REPORT {top_dir}/{sample_name}/TCR_unfiltered/barcode_report.tsv
    SF_TCR_UF_CLON_FREQ {top_dir}/{sample_name}/TCR_unfiltered/clonotype_frequency.tsv

    # Sample files; Analysis
    SFR_ASUM_CSV        {top_dir}/{sample_name}/report/analysis_summary.csv
    SFR_ANA_PROC_DEF    {top_dir}/{sample_name}/report/analysis_process.json
    SFR_CLUST_DIFF_EXP  {top_dir}/{sample_name}/report/cluster_diff_exp.csv
    SFR_CLUST_ASSIGN    {top_dir}/{sample_name}/report/cluster_assignment.csv
    SFR_CLUST_UMAP      {top_dir}/{sample_name}/report/cluster_umap.csv

    # Sample files; Report context
    SFR_ALLSTATS        {top_dir}/{sample_name}/report/sample_all_stats.csv
    SFR_SAMP_DEF        {top_dir}/{sample_name}/report/sample_status.json
    SFR_TSCP_CUTOFF     {top_dir}/{sample_name}/report/tscp_cutoff_calc.csv
    SFR_TSCP_CT         {top_dir}/{sample_name}/report/tscp_counts.csv
    SFR_EXPRESS_GENE    {top_dir}/{sample_name}/report/expressed_genes.csv
    SFR_ENRICHMENT      {top_dir}/{sample_name}/report/target_enrichment.csv
    SFR_MAP_TARG_CT     {top_dir}/{sample_name}/report/map_targ_counts.csv
    SFR_SPEC_READ_CT    {top_dir}/{sample_name}/report/species_read_counts.csv
    SFR_SPEC_TSCP_CT    {top_dir}/{sample_name}/report/species_tscp_counts.csv
    SFR_SS_TSCP_CT      {top_dir}/{sample_name}/report/tscp_counts_subsampled.csv
    SFR_SS_GENE_CT      {top_dir}/{sample_name}/report/gene_counts_subsampled.csv
    SFR_SAMP_MRTG_CT    {top_dir}/{sample_name}/report/sample_mrtg_counts.csv
    SFR_TSCP_R1W        {top_dir}/{sample_name}/report/tscp_median_by_rnd1_well.csv
    SFR_TSCP_R2W        {top_dir}/{sample_name}/report/tscp_median_by_rnd2_well.csv
    SFR_TSCP_R3W        {top_dir}/{sample_name}/report/tscp_median_by_rnd3_well.csv
    SFR_CELL_R1W        {top_dir}/{sample_name}/report/cell_counts_by_rnd1_well.csv
    SFR_CELL_R2W        {top_dir}/{sample_name}/report/cell_counts_by_rnd2_well.csv
    SFR_CELL_R3W        {top_dir}/{sample_name}/report/cell_counts_by_rnd3_well.csv
    # Focal and crispr files
    SFR_FB_CELL_CT      {top_dir}/{sample_name}/report/guide_cell_counts.csv
    SFR_FB_CELL_MAP_TAB {top_dir}/{sample_name}/focal_filtered/guide_assignment.csv
    SFR_FB_CELL_MAP_DIC {top_dir}/{sample_name}/focal_filtered/guide_assignment.json
    SFR_CRSP_CELL_CT    {top_dir}/{sample_name}/report/guide_cell_counts.csv
    SFR_CRSP_CELL_MAP_TAB {top_dir}/{sample_name}/guide_RNAs_filtered/guide_assignment.csv
    SFR_CRSP_CELL_MAP_DIC {top_dir}/{sample_name}/guide_RNAs_filtered/guide_assignment.json

    # Individual figures
    SFR_FIG_CELL_CUTOFF {top_dir}/{sample_name}/figures/fig_tscp_cell_cutoff.png
    SFR_FIG_TSCP_R1W    {top_dir}/{sample_name}/figures/fig_tscp_by_rnd1_well.png
    SFR_FIG_TSCP_R2W    {top_dir}/{sample_name}/figures/fig_tscp_by_rnd2_well.png
    SFR_FIG_TSCP_R3W    {top_dir}/{sample_name}/figures/fig_tscp_by_rnd3_well.png
    SFR_FIG_CELL_R1W    {top_dir}/{sample_name}/figures/fig_cell_by_rnd1_well.png
    SFR_FIG_CELL_R2W    {top_dir}/{sample_name}/figures/fig_cell_by_rnd2_well.png
    SFR_FIG_CELL_R3W    {top_dir}/{sample_name}/figures/fig_cell_by_rnd3_well.png
    SFR_FIG_SS_TSCP     {top_dir}/{sample_name}/figures/fig_ss_tscp_per_cell.png
    SFR_FIG_SS_GENE     {top_dir}/{sample_name}/figures/fig_ss_gene_per_cell.png
    SFR_FIG_UMAP_CLUS   {top_dir}/{sample_name}/figures/fig_umap_cluster.png
    SFR_FIG_UMAP_SAMP   {top_dir}/{sample_name}/figures/fig_umap_sample.png
    SFR_FIG_UMAP_SPEC   {top_dir}/{sample_name}/figures/fig_umap_species.png
    SFR_FIG_UMAP_GENE   {top_dir}/{sample_name}/figures/fig_umap_gene_count.png
    SFR_FIG_UMAP_TSCP   {top_dir}/{sample_name}/figures/fig_umap_tscp_count.png
    SFR_FIG_UMAP_READ   {top_dir}/{sample_name}/figures/fig_umap_read_count.png
    SFR_FIG_BARNYARD    {top_dir}/{sample_name}/figures/fig_barnyard.png
    """


def get_spipe_filepath_dict():
    """Get dict with filepath format templates

    Strings are used by 'spipe.filepath(key, sample)' function

    Return dict
    """
    fp_dict = {}
    for line in FILE_DEFS.split("\n"):
        # Split off anything after comment, then split on words
        parts = line.split("#")[0].split()
        if len(parts) == 2:
            fp_dict[parts[0]] = parts[1]

    return fp_dict


# ----------------- Executables ---------------------------------------------
# Initial names to call; Paths get expanded prior to use
EXEC_DEFS = """
    # Data lines: <exe> <arg> <num> [<opt>]
    # Executable, args for version, line for version output, optional status
    STAR        --version   1
    samtools    --version   1
    fastqc      --version   1   optional
    gzip        --version   1
    pigz        --version   1
    unzip       -v          1
    """


# To hold exec info
ExecInfo = namedtuple(
    "ExecInfo", ["name", "status", "optional", "path", "ver", "story"]
)


def get_spipe_exec_dict(call=True, ignore=None, timeout=None):
    """Get dict with exec paths

    call = flag to actually call executables (i.e. to test out)
    ignore = possible list of names to ignore
    timeout = passed to execute call

    Return dict [name] = ExecInfo
    """
    exe_dict = {}
    for line in EXEC_DEFS.split("\n"):
        # Split off anything after comment, then split on words
        parts = line.split("#")[0].split()
        if len(parts) < 2:
            continue

        # name, version call arg, version line, optionality flag
        name = parts[0]
        if ignore and name in ignore:
            continue
        ver_call = parts[1]
        ver_line = 1 if len(parts) < 3 else int(parts[2])

        # Optional?
        is_opt = False
        if (len(parts) > 3) and (parts[3].lower().startswith("opt")):
            is_opt = True

        # Calling to get version?
        if call:
            # Get path and version
            ok, path, version, story = utils.get_exec_path_and_ver(
                name, ver_call, ver_line=ver_line, timeout=timeout
            )
        else:
            ok, path, story = utils.check_exec(name, timeout=timeout)
            version = "Version not checked"

        exe_dict[name] = ExecInfo(name, ok, is_opt, path, version, story)

    return exe_dict


def use_case_out_dir_key(spipe, use_case="", top_dir=None, as_list=False):
    """
    Return file key or list of keys for use case
    """
    if not use_case:
        if top_dir:
            use_case, _ = top_dir.get_runproc_info(key="run.use_case")
        else:
            use_case = spipe.get_use_case()

    # Make lists of DGE keys; Primary (analysis + report) first, then...
    if use_case == "normal":
        dir_key = "DIR_S_DGE_N_F,DIR_S_DGE_N_UF"
    elif use_case == "target_only":
        dir_key = "DIR_S_DGE_T_F,DIR_S_DGE_T_UF"
    elif use_case == "target_parent":
        dir_key = "DIR_S_DGE_T_F,DIR_S_DGE_T_UF"
    elif use_case == "focal_bc":
        dir_key = "DIR_S_DGE_FB_F,DIR_S_DGE_FB_UF"
    elif use_case == "crispr":
        dir_key = "DIR_S_DGE_CRSP_F,DIR_S_DGE_CRSP_UF"
    elif use_case == "tcr_only":
        dir_key = "DIR_S_TCR_UF"
    elif use_case == "tcr_parent":
        dir_key = "DIR_S_TCR_F,DIR_S_TCR_UF"
    elif use_case == "mkref":
        dir_key = ""
    else:
        raise ValueError(f"Unknown use_case {use_case} for dir keys")
    # As list or just primary key
    dir_key = dir_key.split(",")
    if as_list:
        dir_key = [k for k in dir_key if k]
    else:
        dir_key = dir_key[0]
    return dir_key


def out_filepaths(
    spipe, samp, dir_key="", use_case=True, targeted=False, top_dir="", as_key=False
):
    """Get filepath names or keys for DGE given sample and type

    dir_key = DGE top level key
    use_case = Flag to base answer on use_case if no dir_key is given
    targeted = Flag for targeted vs not
    top_dir = top dir (obj or str) to get use_case if using that
    as_key  = Flag to return file keys, not full paths

    Return dict with file paths/keys
    """
    c_use_case = ""
    # If no given key, assign via use case
    if (not dir_key) and use_case:
        c_use_case = spipe.get_use_case(top_dir=top_dir)
        dir_key = use_case_out_dir_key(spipe, use_case=c_use_case)

    # Source gene (i.e. not output; Input from genome or targeted list)
    src_gene = "PF_TARGET_GENES" if targeted else "GENO_GENE_DATA"

    # Default base case (filtered use_case 'normal')
    file_keys = {
        "dir": dir_key,
        "mtx": "SF_NF_MTX",
        "gene": "SF_NF_GENE",
        "cell": "SF_NF_CELL",
        "anndata": "SF_NF_ANA_H5AD",
        "src_gene": src_gene,
    }
    # Update to subdir-specific file keys
    for k, v in file_keys.items():
        if dir_key == "DIR_S_DGE_N_F":
            pass
        elif dir_key == "DIR_S_DGE_N_UF":
            v = v.replace("SF_NF_", "SF_NU_")
        # target
        elif dir_key == "DIR_S_DGE_T_F":
            v = v.replace("SF_NF_", "SF_TF_")
        elif dir_key == "DIR_S_DGE_T_UF":
            v = v.replace("SF_NF_", "SF_TU_")
        # focal barcoding
        elif dir_key == "DIR_S_DGE_FB_F":
            v = v.replace("SF_NF_", "SF_FF_")
        elif dir_key == "DIR_S_DGE_FB_UF":
            v = v.replace("SF_NF_", "SF_FU_")
        # crispr
        elif dir_key == "DIR_S_DGE_CRSP_F":
            v = v.replace("SF_NF_", "SF_CF_")
        elif dir_key == "DIR_S_DGE_CRSP_UF":
            v = v.replace("SF_NF_", "SF_CU_")
        else:
            raise ValueError(
                f"Unknown DGE output dir key {dir_key}, use_case {c_use_case}"
            )
        file_keys[k] = v
    # Expand to actual paths?
    if not as_key:
        for k, v in file_keys.items():
            # Check if key exists
            if v in spipe._filepath_dict:
                fpath = spipe.filepath(v, samp, top_dir=top_dir)
            else:
                fpath = None
            file_keys[k] = fpath
    return file_keys


def geno_filepaths(
    spipe, samp, dir_key="", use_case=True, targeted=False, top_dir="", as_key=False
):
    """ """
    new_dict = {}
    fkey = "PF_TARGET_GENES" if targeted else "PF_ALL_GENES"

    new_dict["genes"] = fkey
    if not as_key:
        new_dict["genes"] = spipe.filepath(new_dict["genes"], samp, top_dir=top_dir)

    return new_dict
