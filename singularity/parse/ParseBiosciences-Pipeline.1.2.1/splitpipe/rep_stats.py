#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Report (HTML) stats prep / formatting stuff
#

from collections import OrderedDict


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


### ------------------- Stats key-value list for summary -------------------


def get_report_stats_kv_list(spipe, samp, stat_df, targeted=False, focal=False):
    """Get (html) report stats with pretty print label and value ready for jinja

    Order of stats (name, val) in returned list = display order on page

    Return list of (key,val) stats in order for reporting on html
    """
    # Make dict from dataframe first col
    stat_dict = stat_df.iloc[:, 0].to_dict()
    # List of key:values to be used
    stat_list = []

    stat_list += get_cell_num_kv_list(spipe, stat_dict, targeted=targeted)

    if focal:
        stat_list += get_focal_tscp_gene_kv_list(spipe, stat_dict, targeted=targeted)
    else:
        stat_list += get_tscp_gene_kv_list(spipe, stat_dict, targeted=targeted)

    stat_list += get_reads_kv_list(spipe, stat_dict, targeted=targeted)
    if focal:
        stat_list += get_focal_kv_list(spipe, stat_dict, targeted=targeted)
    if targeted:
        stat_list += get_target_kv_list(spipe, stat_dict)
    stat_list += get_seq_kv_list(spipe, stat_dict, targeted=targeted)

    return stat_list


def get_cell_num_kv_list(spipe, stat_dict, targeted=False):
    """Get the cell number parts of report into key:val dict list

    stat_dict = dict with all stats
    targeted = flag for targeted case

    Return list of k:v dicts
    """
    # Base key for number of cells
    xkey = "targeted_" if targeted else ""
    b_key = f"{xkey}number_of_cells"

    # Total first
    new_dict = OrderedDict()
    new_dict["Estimated Number of Cells"] = stat_dict[b_key]

    # Only if multiple genomes
    if spipe.num_genomes() > 1:
        # Collect species specific keys; Multiplet is last
        k_list = []
        mult_key = ""
        for k in sorted(stat_dict.keys()):
            if k == b_key:
                continue
            if k.startswith("multiplet_"):
                mult_key = k
                continue
            if k.endswith(b_key):
                k_list.append(k)
        k_list.append(mult_key)
        # Final format
        for k in k_list:
            species = k.split("_")[0]
            new_dict[f"{species} Number of Cells Detected"] = stat_dict[k]

    stat_list = utils.dict_to_kv_dict_list(new_dict)
    return stat_list


def get_tscp_gene_kv_list(spipe, stat_dict, targeted=False):
    """Get the transcript and gene parts of report into key:val dict list

    stat_dict = dict with all stats
    targeted = flag for targeted case

    Return list of k:v dicts
    """
    # Pretty print extra keyword
    p_xkey = "Targeted " if targeted else ""

    new_dict = OrderedDict()
    # Two pass, tscp first
    for k in sorted(stat_dict.keys()):
        if k.endswith("_median_tscp_per_cell"):
            # Only preface with species if more than one
            if spipe.num_genomes() > 1:
                species = k.split("_")[0]
                new_dict[f"{species} Median {p_xkey}Transcripts/Cell"] = stat_dict[k]
            else:
                new_dict[f"Median {p_xkey}Transcripts/Cell"] = stat_dict[k]

    for k in sorted(stat_dict.keys()):
        if k.endswith("_median_genes_per_cell"):
            if spipe.num_genomes() > 1:
                species = k.split("_")[0]
                new_dict[f"{species} Median {p_xkey}Genes/Cell"] = stat_dict[k]
            else:
                new_dict[f"Median {p_xkey}Genes/Cell"] = stat_dict[k]

    stat_list = utils.dict_to_kv_dict_list(new_dict)
    return stat_list


def get_focal_tscp_gene_kv_list(spipe, stat_dict, targeted=False):
    """Focal median tscp per cell

    Return list of k:v dicts
    """
    xkey = "targeted_" if targeted else ""

    species = spipe.get_genome_list()[0]
    stat_key = f"{species}_{xkey}median_tscp_per_cell"

    new_dict = OrderedDict()
    new_dict["Median Guide Transcripts/Cell"] = stat_dict[stat_key]

    stat_list = utils.dict_to_kv_dict_list(new_dict)
    return stat_list


def get_reads_kv_list(spipe, stat_dict, targeted=False):
    """Get the read parts of report into key:val dict list

    stat_dict = dict with all stats
    targeted = flag for targeted case

    Return list of k:v dicts
    """
    # Extra key if targeted
    xkey = "targeted_" if targeted else ""

    new_dict = OrderedDict()

    if spipe.is_focal_crispr(crispr=True):
        new_dict["Guide Mean Reads/Cell"] = (
            stat_dict["number_of_reads"] / stat_dict[f"{xkey}number_of_cells"]
        )
        new_dict["Guide Number of Reads"] = stat_dict["number_of_reads"]
    else:
        new_dict["Mean Reads/Cell"] = (
            stat_dict["number_of_reads"] / stat_dict[f"{xkey}number_of_cells"]
        )
        new_dict["Number of Reads"] = stat_dict["number_of_reads"]

    stat_list = utils.dict_to_kv_dict_list(new_dict)
    return stat_list


def get_focal_kv_list(spipe, stat_dict, targeted=False):
    """Get the focal / crispr parts of report into key:val dict list

    stat_dict = dict with all stats

    Return list of k:v dicts
    """
    new_dict = OrderedDict()
    xkey = "targeted_" if targeted else ""

    n_cells = stat_dict[f"{xkey}number_of_cells"]
    n_any = n_cells - stat_dict["num_cells_with_0_crispr_guide"]
    # new_dict['Fraction Cells with Focal Genes'] = n_any / n_cells
    # new_dict['Usable Guide Fraction'] = n_any / n_cells
    new_dict["Fraction of Cells with Guides"] = n_any / n_cells

    new_dict["Number of Cells with 1 Guide"] = stat_dict[
        "num_cells_with_1_crispr_guide"
    ]
    new_dict["Number of Cells with >= 2 Guides"] = stat_dict[
        "num_cells_with_2plus_crispr_guide"
    ]

    stat_list = utils.dict_to_kv_dict_list(new_dict)
    return stat_list


def get_target_kv_list(spipe, stat_dict):
    """Get the targeted parts of report into key:val dict list

    stat_dict = dict with all stats

    Return list of k:v dicts
    """
    new_dict = OrderedDict()
    new_dict["Fraction Reads on Target"] = stat_dict["targeted_read_fraction"]
    new_dict["Targeted Number of Genes"] = stat_dict["target_number"]
    new_dict["Fraction of Target Genes Found"] = stat_dict["targeted_gene_fraction"]

    stat_list = utils.dict_to_kv_dict_list(new_dict)
    return stat_list


def get_seq_kv_list(spipe, stat_dict, targeted=False):
    """Get the sequencing parts of report into key:val dict list

    stat_dict = dict with all stats

    Return list of k:v dicts
    """
    # Extra key and pretty keys if targeted or focal
    xkey = "targeted_" if targeted else ""
    p_xkey = "Targeted " if targeted else ""
    f_xkey = "Guide " if spipe.is_focal_crispr(crispr=True) else ""

    new_dict = OrderedDict()
    new_dict[f"{p_xkey}{f_xkey}Sequencing Saturation"] = stat_dict[
        f"{xkey}sequencing_saturation"
    ]
    new_dict[f"{f_xkey}BC1 (RT) >Q30"] = stat_dict["bc1_Q30"]
    new_dict[f"{f_xkey}BC2 >Q30"] = stat_dict["bc2_Q30"]
    new_dict[f"{f_xkey}BC3 >Q30"] = stat_dict["bc3_Q30"]
    new_dict[f"{f_xkey}cDNA >Q30"] = stat_dict["cDNA_Q30"]

    stat_list = utils.dict_to_kv_dict_list(new_dict)
    return stat_list
