#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Report (HTML) generation stuff
#

import sys
import re
import traceback
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as tickr
import matplotlib.colors as mc
import colorsys
import seaborn as sns
import base64
import jinja2
import multiprocessing as mp
from collections import OrderedDict
from scipy.sparse import csc_matrix
from natsort import natsort_keygen


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import pb_const as const
    import rep_pp
    import rep_stats
    import utils
else:
    from splitpipe import pb_const as const
    from splitpipe import rep_pp
    from splitpipe import rep_stats
    from splitpipe import utils


# ---------------------------------------------------------------------------
def have_reports_ins(spipe, verb=True, samp=None):
    """Check if reports inputs exist

    Return True/False, list of files
    """
    ok = False
    check_files = []
    # Files for given or all samples
    if samp:
        # All stats
        check_files.append(spipe.filepath("SFR_ALLSTATS", samp))
        # Analysis (clustering) info; Not for TCR
        if not spipe.is_tcr():
            check_files.append(spipe.filepath("SFR_ANA_PROC_DEF", samp))
    else:
        for samp in spipe.get_samples():
            check_files.append(spipe.filepath("SFR_ALLSTATS", samp))
            if not spipe.is_tcr():
                check_files.append(spipe.filepath("SFR_ANA_PROC_DEF", samp))
    # If list is good, we have inputs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def have_reports_outs(spipe, verb=False, samp=None):
    """Check if reports output files exist

    Return True/False, list of files
    """
    ok = False
    check_files = []
    # HTML files
    if samp:
        check_files.append(spipe.filepath("SF_ASUM_HTML", samp))
    else:
        for samp in spipe.get_samples():
            check_files.append(spipe.filepath("SF_ASUM_HTML", samp))
    # If list is good, we have outs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def run_reports(spipe):
    """Run reports step

    Return status
    """
    spipe.report_proc_step("Reports")
    # Input and output status and list
    # If keep_going, input checks are per sample (ignore here)
    if spipe.keep_going():
        i_ok = True
    else:
        i_ok, i_list = have_reports_ins(spipe)
    o_ok, o_list = have_reports_outs(spipe)
    # Run only if don't have output or fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for report: {bad_list}"
            spipe.set_problem(story)
        else:
            generate_all_html_reports(spipe)
    else:
        spipe.report_run_story("Skipping Reports and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")

    clean_up(spipe)
    ok = spipe.no_problems()
    spipe.report_proc_step("Reports", status=ok)
    return ok


def generate_all_html_reports(spipe):
    """Function to loop through every sample and process

    Return nothing; If any problems, these are reported / remembered
    """
    # Each sample
    samp_list = spipe.get_samples()
    for i, samp in enumerate(samp_list):
        ok = True
        spipe.report_run_story("")
        story = (
            f"Report processing sample ({i+1} of {len(samp_list)}) {samp.get_name()}"
        )
        spipe.report_run_story(story)

        # If already not good, skip
        if not samp.get_status(is_good=True):
            stat_step = samp.get_status(step=True)
            story = f"Status of sample {samp.get_name()} bad ({stat_step}); Skipping"
            spipe.report_run_story(story)
            continue

        spipe.report_run_mem(story=f"Report samp {i+1}")

        # Have required inputs?
        i_ok, i_list = have_reports_ins(spipe, samp=samp)
        if not i_ok:
            ok = False
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for report: {bad_list}"
            spipe.set_problem(story)
        else:
            # Catch possible problems (e.g. plotting, scanpy ...)
            try:
                if spipe.is_tcr():
                    ok = generate_single_tcr_html_report(spipe, samp)
                else:
                    ok = generate_single_html_report(spipe, samp)
            except Exception as e:
                ok = False
                # Dump traceback, stdout and log
                traceback.print_exc(file=sys.stdout)
                sys.stdout.flush()
                traceback.print_exc(file=spipe.get_log_fh())
                story = f"Failed generate_single_html_report: {e}"
                spipe.set_problem(story)

        # Handle possible post processing (e.g. merging new with parent html)
        if ok:
            ok = rep_pp.post_proc_html_report(spipe, samp)

        # Problem? report / possibly bail
        if not ok:
            story = f"Report problem with sample {samp.get_name()}"
            spipe.set_problem(story)
            samp.set_status(False, "report")
            if spipe.keep_going():
                spipe.report_run_story("Keep-going True... So, keep on going!")
                spipe.clear_problems()
                print()
            else:
                break


def generate_single_html_report(spipe, samp):
    """Generate html report for one sample; Content depends on loaded files

    Returns status; Outputs saved to files
    """
    # Various things depend on targeted and/or focal barcoding
    targeted = spipe.is_targeted()
    focal_crispr = spipe.is_focal_crispr()

    # Load info and anndata with clusters, etc
    ana_info, adata = load_samp_analysis(spipe, samp)
    # Have to at least have a status (i.e. file from analysis step)
    if (not ana_info["status"]) or (adata is None):
        spipe.set_problem(f"No analysis done for sample {samp.get_name()}")
        return False

    min_bc = ana_info["cluster_min_bc"]
    n_cells, n_genes = adata.shape

    # Stats; All stats dataframe
    stat_df = spipe.get_stats_data(samp)
    if spipe.have_parent():
        # Parent dir data
        par_dir = spipe.get_parent_info(key="path")
        par_stat_df = spipe.get_stats_data(samp, top_dir=par_dir)
    else:
        par_stat_df = None

    # Package jinja vars for pipeline
    pipeline_json = get_pipeline_json(spipe, samp)

    # statistics collection for report display
    spipe.report_run_story2("Collecting stats to report")
    stat_list = rep_stats.get_report_stats_kv_list(
        spipe, samp, stat_df, targeted=targeted, focal=focal_crispr
    )
    sum_json = {
        "stat_list": stat_list,
        "cur_samp_name": samp.get_name(),
        "cur_samp_wells": samp.get_wspec(),
        "cur_samp_nsublibs": samp.num_sublibs(),
    }

    # Init jinja vars; Some to dummy values
    umap_json = "umap_json"
    diff_exp_table_html = "diff_exp_table_html"
    clust_names = "clust_names"

    # ============ Independent of clusering ==========
    spipe.report_run_story2("Generating plate figures")
    plates_json = plot_plates(spipe, samp, targeted=targeted, focal=focal_crispr)

    # Sequencing saturation figures
    spipe.report_run_story2("Generating subsample plots")
    tscp_json = plot_gt_subsamp(spipe, samp, tscp=True, targeted=targeted)
    gene_json = plot_gt_subsamp(spipe, samp, tscp=False, targeted=targeted)
    sat_json = {"tscp_sat": tscp_json, "gene_sat": gene_json}

    # Get values for ploty figures on pg1 of report
    # "Identified cells" Transcript -vs- sorted BC fig
    spipe.report_run_story2("Generating tscp cell count figure")
    bc_rank_json = plot_bc_rank(
        spipe, samp, stat_df, targeted=targeted, par_stat_df=par_stat_df
    )

    # Barnyard plot
    if spipe.num_genomes() > 1:
        if spipe.get_par_val("rep_save_figs", as_bool=True):
            spipe.report_run_story2("Generating barnyard plot")
            plot_barnyard(spipe, samp)

    # This communicates issues to change HTML; too few cells / genes for plot, etc
    warning_json = {
        "status": "good",
        "n_cells": n_cells,
        "n_genes": n_genes,
        "min_bc": min_bc,
    }

    # Use cluster number to determine if we have cluster data
    # Also focal / crispr, getting data from parent
    if focal_crispr or (int(ana_info.get("cluster_num", 0)) > 0):
        # Doesn't return anything; Table written to temporary file
        spipe.report_run_story2("Rendering gene enrichment table (pb)")
        diff_exp_table_html = render_clus_gene_tab(spipe, samp)

        # Cell subsample for UMAP plots; If no subsampling, adata is just returned
        if focal_crispr:
            # Cell subsamp based on parent
            adata = get_umap_pl_par_cell_ss(spipe, samp, adata)
        else:
            adata = get_umap_pl_cell_subsamp(spipe, samp, adata)

        # Get UMAP plotting info (excep gex) into dataframe
        if focal_crispr:
            par_ana_info, par_adata = load_samp_analysis(spipe, samp, parent=True)
            umap_df = get_umap_df(
                spipe, samp, par_adata, par_ana_info, targeted=targeted
            )
        else:
            umap_df = get_umap_df(spipe, samp, adata, ana_info, targeted=targeted)

        xkey = "targeted_" if targeted else ""
        count_cols = [f"{xkey}gene_count", f"{xkey}tscp_count", f"{xkey}mread_count"]
        umap_df[count_cols] = adata.obs[count_cols]

        umap_json = plot_umap(spipe, samp, umap_df, adata, ana_info, targeted=targeted)

        # Cluster 'names' for selection legend; Unique (short, sorted list)
        clust_names = get_cell_cluster_names(spipe, samp, adata, ana_info, unique=True)
    else:
        if "too_few" in ana_info["status"]:
            warning_json["status"] = "cellwarning"
            story = f"Too few cells; No clustering: {adata.shape}"
            spipe.report_run_story(story)
        else:
            warning_json["status"] = "clustwarning"
            story = f"No clustering data: {ana_info['status']}"
            spipe.report_run_story(story)

    j_args = {
        "pipeline_json": pipeline_json,
        "diff_exp_table_html": diff_exp_table_html,
        "bc_rank_json": bc_rank_json,
        "warning_json": warning_json,
        "clust_names": clust_names,
        "plates_json": plates_json,
        "umap_json": umap_json,
        "sum_json": sum_json,
        "sat_json": sat_json,
    }
    output_jinja_html(spipe, samp, j_args)

    return True


def output_jinja_html(spipe, samp, j_args):
    """Call Jinja2 to process then save HTML

    j_args = dict with variables to pass to Jinja2

    Return status
    """
    # output filename and pkg path
    rendered_file_path = spipe.filepath("SF_ASUM_HTML", samp)
    script_path = spipe.get_pkg_path()

    # Input HTML templates for jinja2
    if spipe.is_tcr():
        template_filename = "/templates/template_tcr.html.j2"
    else:
        template_filename = "templates/template.html.j2"

    # Global flag var for alt styling; Lowercase javascript
    alt_style_fmt = (
        "true" if spipe.get_par_val("rep_alt_sty_fmt", as_bool=True) else "false"
    )
    # Update jinja arg collection
    j_args["PATH"] = script_path
    j_args["alt_style_fmt"] = alt_style_fmt
    j_args["logo"] = get_logo_img(spipe)
    j_args["checkbox"] = get_chk_img(spipe)

    # Jinja2 calls
    environment = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath=script_path)
    )
    output_text = environment.get_template(template_filename).render(j_args)

    with open(rendered_file_path, "w") as result_file:
        result_file.write(output_text)
        spipe.report_run_story(f"New report html: {rendered_file_path}")

    return True


def get_pipeline_json(spipe, samp):
    """Get pipeline / run setting / version variables into json stucture

    Return dict
    """
    # Only get the number part of version (i.e. 'split-pipe v123' >--> v123)
    pipeline_version = spipe.get_version().split()[-1]
    chem = spipe.get_chemistry()
    kit_s, kit_n, _ = spipe.get_kit()
    use_case = spipe.get_use_case()
    # footer = f"Parse Biosciences - {kit_s} {chem}, Pipeline {pipeline_version}"
    footer = f"Parse Biosciences - Pipeline {pipeline_version}"
    # Format flag
    alt_style_fmt = (
        "true" if spipe.get_par_val("rep_alt_sty_fmt", as_bool=True) else "false"
    )

    # Important: Javascript booleans are lowercase
    is_targeted = "true" if spipe.is_targeted() else "false"
    is_comb = "true" if spipe.is_combine() else "false"
    is_tcr_only = "true" if spipe.is_tcr(only=True) else "false"
    is_tcr_parent = "true" if spipe.is_tcr(parent=True) else "false"
    is_focal = "true" if spipe.is_focal_crispr(focal=True) else "false"
    is_crispr = "true" if spipe.is_focal_crispr(crispr=True) else "false"
    is_meta_sample = "true" if samp.is_meta() else "false"
    have_parent = "true" if spipe.have_parent() else "false"

    pipeline_json = {
        "version_footer": footer,
        "pipeline_version": pipeline_version,
        "chemistry": chem,
        "kit": kit_s,
        "use_case": use_case,
        "is_targeted": is_targeted,
        "is_comb": is_comb,
        "is_tcr_only": is_tcr_only,
        "is_tcr_parent": is_tcr_parent,
        "is_focal": is_focal,
        "is_crispr": is_crispr,
        "is_meta_sample": is_meta_sample,
        "have_parent": have_parent,
        "alt_style_fmt": alt_style_fmt,
    }
    return pipeline_json


### ------------------------- TCR specific ------------------------------


def generate_single_tcr_html_report(spipe, samp):
    """Generate TCR html report for one sample

    Returns status; Outputs saved to files
    """

    # Stats; All stats dataframe
    if spipe.have_parent():
        filtered = True
    else:
        filtered = False

    # Package jinja vars for pipeline
    pipeline_json = get_pipeline_json(spipe, samp)

    # statistics collection for report display
    spipe.report_run_story2("Collecting stats to report")
    sum_json = get_tcr_sum_stat_json(spipe, samp)

    spipe.report_run_story2(f"Collecting clonotype data (filtered {filtered})")
    clonotype_dict, clonotype_table_html = get_tcr_clonotype_table(
        spipe, samp, filtered=filtered
    )

    # This communicates issues to change HTML; too few clonotypes for plot, etc
    n_clonotypes = clonotype_dict["n_clonotypes"]
    status = "good" if n_clonotypes > 0 else "clonowarning"

    warning_json = {"status": status, "n_clonotypes": n_clonotypes}

    j_args = {
        "sum_json": sum_json,
        "pipeline_json": pipeline_json,
        "warning_json": warning_json,
        "clonotype_dict": clonotype_dict,
        "clonotype_table_html": clonotype_table_html,
    }
    output_jinja_html(spipe, samp, j_args)

    return True


def get_tcr_sum_stat_json(spipe, samp, round_to=3):
    """Get TCR summary stat collection

    Return json dict
    """
    stat_df = spipe.get_stats_data(samp)
    raw_dict = stat_df.iloc[:, 0].to_dict()

    stat_dict = OrderedDict()
    stat_dict["Estimated Number of T Cells"] = raw_dict["tcr_number_of_cells"]
    stat_dict["Mean Reads/Cell"] = raw_dict["tcr_mean_reads_per_cell"]
    # stat_dict['Median VDJ Transcripts/Cell'] = raw_dict['tcr_median_VDJ_tscp_per_cell']

    # As fraction of all cells
    nc = raw_dict["tcr_number_of_cells"]
    stat_dict["T Cells with Productive TRA"] = round(
        raw_dict["tcr_tcells_with_productive_TRA"] / nc, round_to
    )
    stat_dict["T Cells with Productive TRB"] = round(
        raw_dict["tcr_tcells_with_productive_TRB"] / nc, round_to
    )
    stat_dict["T Cells with Detected Productive Pair (TRA, TRB)"] = round(
        raw_dict["tcr_tcells_with_productive_pair"] / nc, round_to
    )

    # stat_dict['T Cells with Inferred Productive Pair (TRA, TRB)'] = round(raw_dict['tcr_tcells_with_productive_inferred_pair'] / nc, round_to)

    stat_dict["BC1 >Q30"] = raw_dict["bc1_Q30"]
    stat_dict["BC2 >Q30"] = raw_dict["bc2_Q30"]
    stat_dict["BC3 >Q30"] = raw_dict["bc3_Q30"]
    stat_dict["cDNA >Q30"] = raw_dict["cDNA_Q30"]

    # Get list of kev-val dicts for json to javascript
    stat_list = utils.dict_to_kv_dict_list(stat_dict)

    sum_json = {
        "stat_list": stat_list,
        "cur_samp_name": samp.get_name(),
        "cur_samp_wells": samp.get_wspec(),
        "cur_samp_nsublibs": samp.num_sublibs(),
    }
    return sum_json


def get_tcr_clonotype_table(spipe, samp, filtered=False):
    """Get clonotype data

    Return data dict, html content string
    """
    # File key depends on filtered or not
    fkey = "SF_TCR_F_CLON_FREQ" if filtered else "SF_TCR_UF_CLON_FREQ"
    # File should be present even if zero data
    clonotype_df = spipe.read_csv(fkey, samp, sep="\t", index_col=None)
    n_clonotypes = len(clonotype_df)

    n_clonotype_to_plot = min(10, n_clonotypes)
    n_clonotype_table = min(100, n_clonotypes)

    # Possibly has no clonotypes
    clonotype_df = spipe.read_csv(fkey, samp, sep="\t", index_col=None)

    # Have data
    if n_clonotypes > 0:
        # Reduce decimals and
        clonotype_df["frequency"] = clonotype_df["frequency"].round(4)
        clonotype_df["TRA"] = clonotype_df["TRA"].fillna("-")
        clonotype_df["TRB"] = clonotype_df["TRB"].fillna("-")

        clonotype_table = clonotype_df[
            ["clonotype_id", "TRA", "TRB", "count", "frequency"]
        ]
        clonotype_table.columns = [["Clonotype ID", "TRA", "TRB", "Count", "Frequency"]]
        clonotype_table_html = clonotype_table[:n_clonotype_table].to_html(
            border=0, index=False, classes=["diff-table"]
        )
        clonotype_table_html = clonotype_table_html.replace("dataframe", "table")
    else:
        clonotype_table_html = ""

    # Reverse order so highest clonotype is on top
    clonotype_df = clonotype_df[:n_clonotype_to_plot][::-1]

    hover_text = []
    for i in range(n_clonotype_to_plot):
        clonotype_text = str(
            clonotype_df.iloc[i].loc[["TRA", "TRB", "count"]].to_dict()
        )[1:-1]
        # Plotly uses breaks not \n
        clonotype_text = clonotype_text.replace(", ", "<br />")
        # Remove excessive single quotes
        clonotype_text = clonotype_text.replace("'", "")
        hover_text.append(clonotype_text)

    clonotype_dict = {
        "n_clonotypes": n_clonotypes,
        "clonotype_id": list(clonotype_df.clonotype_id.values),
        "frequency": list(clonotype_df.frequency.values),
        "hover_text": hover_text,
    }

    return clonotype_dict, clonotype_table_html


### --------------------------- Cluster analysis ----------------------------


def load_samp_analysis(spipe, samp, parent=False):
    """Load data from analysis of one sample

    parent = flag to load from parent_dir

    Return dict with info and possibly andata
    """
    # If parent_dir (e.g. for focal bc)
    top_dir = None
    if parent:
        top_dir = spipe.get_parent_info(key="path")

    # Info file
    fname = spipe.filepath("SFR_ANA_PROC_DEF", samp, top_dir=top_dir)
    if utils.check_infile(fname, verb=False, toxic=False):
        pdict = utils.read_json(fname)
        info = pdict["analysis"]
    else:
        info = {"status": None}

    dge_files = const.out_filepaths(spipe, samp, use_case=True, top_dir=top_dir)
    # anndata file; If status good and file exists
    if utils.check_infile(dge_files["anndata"], verb=False, toxic=False):
        adata = utils.load_parsebio_dge(dge_files["dir"])
        spipe.report_run_story2(f"Loaded anndata {adata.shape} from {dge_files['dir']}")
    else:
        adata = None
        spipe.report_run_story2(f"No analysis anndata found: {dge_files['anndata']}")

    return info, adata


def get_cell_cluster_names(spipe, samp, adata, ana_info, unique=False):
    """Get list of cluster names

    unique = Flag to make unique (and sorted) list (e.g. for legend)

    Return list
    """
    if spipe.is_focal_crispr():
        # Read cluster dataframe
        par_dir = spipe.get_parent_info(key="path")
        fname = spipe.filepath("SFR_CLUST_ASSIGN", samp, top_dir=par_dir)
        df = utils.read_csv(fname, index_col=0)
        clist = df["cluster"]
    else:
        clust_alg = spipe.get_par_val("ana_cluster_alg", as_str=True)
        clist = adata.obs[clust_alg].astype(int)

    # If 1-based and min is zero, shift
    if spipe.get_par_val("ana_cluster_1base", as_bool=True) and (min(clist) == 0):
        clist += 1

    # Simple list of int; If unique, also sort
    if unique:
        clist = sorted(clist.unique(), key=natsort_keygen())
    else:
        clist = clist.values
    # Cast as str
    clist = [str(c) for c in clist]
    return clist


### ------------------------ Page data and plotting -----------------------


def clear_plot_data():
    """Close matplotlib stuff (should release memory)"""
    plt.close()
    plt.clf()
    plt.cla()


def get_logo_img(spipe):
    """Get encoded image of company logo"""
    # Convert company logo to bytecode
    if spipe.get_par_val("rep_alt_sty_fmt", as_bool=True):
        logo_path = spipe.get_pkg_path() + "/templates/" + "Parse_logo-alt.png"
        spipe.report_run_story2(f"Using alt logo {logo_path}")
    else:
        logo_path = spipe.get_pkg_path() + "/templates/" + "Parse_logo.png"
        spipe.report_run_story2(f"Using norm logo {logo_path}")

    with open(logo_path, "rb") as image_file:
        logo = (
            '"data:image/png;base64,'
            + base64.b64encode(image_file.read()).decode()
            + '"'
        )
        return logo


def get_chk_img(spipe):
    """Get encoded image checkboxes for diff table"""
    chk_path = spipe.get_pkg_path() + "/templates/" + "checkbox.png"
    with open(chk_path, "rb") as image_file:
        chk = (
            '"data:image/png;base64,'
            + base64.b64encode(image_file.read()).decode()
            + '"'
        )
        return chk


def adjust_color(color, light=0.5, d_hue=0):
    """Change lightness of given color"""
    try:
        c = mc.cnames[color]
    except Exception:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0] + d_hue, max(0, min(1, light * c[1])), c[2])


def subsamp_plot_color(spec, sublib):
    """Get color for sublibary

    spec = species counter var (0, 1, etc)
    sublib = sublibrary counter var (0, 1, 2, etc)
    """
    if (spec % 2) == 0:
        color = (0.3, 0.75, 0.95)
        light = 0.65
        d_light = 0.08
        d_hue = 0.05
    else:
        color = (0.95, 0.5, 0.1)
        light = 0.85
        d_light = 0.05
        d_hue = 0.02
    color = adjust_color(color, light=light + d_light * sublib, d_hue=d_hue * sublib)
    return color


def plot_plates(spipe, samp, targeted=False, focal=False):
    """Handle plotting of all plate data heatmaps (cell and tscp for each bc round)

    Return json formatted list of heatmap values and layout for each plate
    """
    margin = {"l": 40, "t": 40, "r": 40, "b": 40}
    xaxis = {
        "tickmode": "array",
        "tickvals": list(range(0, 12)),
        "ticktext": list(range(1, 13)),
    }
    # Original color (pre v1.2); Reds
    # colorscale = [[0, "#ffeded"], [1, "#d90000"]]
    # Purple
    # colorscale = [[0, "#f6e6ff"], [1, "#7100b3"]]
    # "Mint-1"
    # colorscale = [[0, "#e6fff9"], [1, "#008062"]]
    # Magenta to Bluish grey; Low saturation
    # colorscale = [[0, "#f5eff5"], [1, "#9f609f"]]
    # colorscale = [[0, "#f4eff5"], [1, "#8d619e"]]
    # colorscale = [[0, "#f2eff6"], [1, "#735194"]]
    # colorscale = [[0, "#f0eff6"], [1, "#675ba4"]]
    # Magenta; old Parse?
    colorscale = [[0, "#f8edf8"], [1, "#943894"]]
    # Muted cyan like disc color pallet
    # colorscale = [[0, "#ecf8f9"], [1, "#37a6ae"]]

    all_plates_json = []
    for rnd in [1, 2, 3]:
        for tscp in [True, False]:
            if tscp:
                p_type = "tscp"
                if targeted:
                    p_title = f"Round {rnd} - Median Targeted Transcripts Per Well"
                    hovertemplate = (
                        "<b>%{y}%{x}</b> Targeted Transcripts: %{z} <extra></extra> "
                    )
                else:
                    if focal:
                        p_title = f"Round {rnd} - Median Guides Per Well"
                        hovertemplate = "<b>%{y}%{x}</b> Guides: %{z} <extra></extra> "
                    else:
                        p_title = f"Round {rnd} - Median Transcripts Per Well"
                        hovertemplate = (
                            "<b>%{y}%{x}</b> Transcripts: %{z} <extra></extra> "
                        )
            else:
                p_type = "cell"
                p_title = f"Round {rnd} - Cells Per Well"
                hovertemplate = "<b>%{y}%{x}</b> Cells: %{z} <extra></extra> "

            # Plate rows depend on well count for round
            # rows, cols = spipe.get_bc_plate_wells(rnd, as_list=True)
            rows, cols = spipe.get_bc_rows_cols(bc_round=rnd, as_list=True)
            # Rows are reversed for plotting, e.g. [D,C,B,A]
            rows = sorted(rows, reverse=True)

            # plot_df = plate_df["count"].values.reshape(len(rows), len(cols))

            # Valid bc / kit have rows = 1, 4, 8
            if len(rows) == 1:
                height = 100
                # Only single row heatmap; Need to extend colorbar, limit ticks
                # Heatmap colorbar parameters:
                #    https://plotly.com/python/reference/#heatmap-colorbar
                # To set colorbar starting at zero, need function on colorscale:
                #   https://community.plotly.com/t/heatmap-color-scale/46282
                #
                nticks = 3
                colorbar_len = 3.5
                yaxis_list = list(range(0, 1))
            elif len(rows) == 4:
                height = 170
                nticks = 4
                colorbar_len = 1
                yaxis_list = list(range(0, 4))
            elif len(rows) == 8:
                height = 245
                nticks = 4
                colorbar_len = 1
                yaxis_list = list(range(0, 8))
            else:
                # print(f"PPLATES nrows = {len(rows)}")
                # print(f"PPLATES {rows}")
                height = 75 + len(rows) * 25
                colorbar_len = max(1, (4.5 - len(rows)))
                nticks = 3 if len(rows) < 3 else 4
                yaxis_list = list(range(0, len(rows)))

            # Get plate data; This also possibly plots figure
            well_z, min_z, max_z = plot_one_plate(spipe, samp, rnd, tscp=tscp)

            # Force zero for color scale?
            if spipe.get_par_val("rep_plate_min_zero", as_bool=True):
                min_z = 0

            # Package for plotly
            yaxis = {
                "tickmode": "array",
                "tickvals": yaxis_list,
                "ticktext": rows,
            }
            well_json = [
                {
                    "type": "heatmap",
                    "z": well_z,
                    "zmin": min_z,
                    "zmax": max_z,
                    "colorscale": colorscale,
                    "colorbar": {
                        "len": colorbar_len,
                        "nticks": nticks,
                    },
                    "hovertemplate": hovertemplate,
                }
            ]
            lyt_json = {
                "height": height,
                "width": 544,
                "xaxis": xaxis,
                "yaxis": yaxis,
                "margin": margin,
                "title": p_title,
            }

            # well_data = {p_type + str(rnd): well_json}
            # well_layout = {"lyt_" + p_type + str(rnd): lyt_json}
            # all_plates_json.extend([well_data, well_layout])
            # well_data = {p_type + str(rnd): well_json}
            # well_layout = {"lyt_" + p_type + str(rnd): lyt_json}
            all_plates_json.extend([well_json, lyt_json])

    return all_plates_json


def plot_one_plate(spipe, samp, rnd, tscp=False):
    """Plot transcript / cell counts per well (i.e. plates)

    rnd = barcoding round (1, 2, 3)
    tscp = flag for transcripts; else cells

    return tuple ([list of well values], min, max)
    """
    # Hardcoded
    colormap = plt.cm.Reds

    if tscp:
        title = f"Round {rnd} - median transcripts per well"
        # file name keys like 'SFR_TSCP_R3W'
        fkey = f"SFR_TSCP_R{rnd}W"
        figkey = f"SFR_FIG_TSCP_R{rnd}W"
    else:
        title = f"Round {rnd} - cells per well"
        # file name keys like 'SFR_CELL_R3W'
        fkey = f"SFR_CELL_R{rnd}W"
        figkey = f"SFR_FIG_CELL_R{rnd}W"

    # Wells named, possibly with plate (e.g. A8, D10, ... p2_B4)
    count_df = spipe.read_csv(fkey, samp=samp)
    count_dict = count_df["count"].to_dict()

    # List of rows and columns (both lists as strings)
    row_list, col_list = spipe.get_bc_rows_cols(bc_round=rnd, as_list=True)
    # Init empty plate then fill with any counts
    plate_df = pd.DataFrame(0, index=row_list, columns=col_list)
    for row in row_list:
        for col in col_list:
            well = row + col
            if well in count_dict:
                plate_df.at[row, col] = count_dict[well]

    yticks = range(0, len(row_list))
    yticklabels = row_list
    xticks = range(0, len(col_list))
    xticklabels = col_list

    # Save figure to file?
    if spipe.get_par_val("rep_save_figs", as_bool=True):
        fig = plt.figure(figsize=(6, 5), dpi=200)
        ax = fig.add_subplot(111)
        cm = plt.imshow(plate_df, cmap=colormap, vmin=0)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_title(title)
        fig.colorbar(cm, pad=0.02, aspect=10, shrink=0.7)
        fname = spipe.filepath(figkey, samp)
        fig.savefig(fname, format="png", bbox_inches="tight", dpi=200)

    clear_plot_data()

    # Data for plotly / json
    plot_mat = []
    for row in row_list:
        plot_mat.append(list(plate_df.loc[row]))

    # Need to reverse rows to display plates correctly
    plot_mat.reverse()

    return (plot_mat, plate_df.values.min(), plate_df.values.max())


def plot_gt_subsamp(spipe, samp, tscp=False, targeted=False):
    """Plot subsample curves of gene / transcript count with reads

    tscp = flag for transcripts; else cells
    targeted = flag for target enrichment

    return dictionary with coordinates for plotly
    """
    # Extra keys and word part of keys if targeted
    p_xkey = "Targeted " if targeted else ""

    # Transcript or cell
    if tscp:
        fkey = "SFR_SS_TSCP_CT"
        xlab = "Sequencing Reads per Cell"
        ylab = f"Median {p_xkey}Transcripts per Cell"
        title = f"{p_xkey}Transcripts per Cell"
        figkey = "SFR_FIG_SS_TSCP"
    else:
        fkey = "SFR_SS_GENE_CT"
        xlab = "Sequencing Reads per Cell"
        ylab = f"Median {p_xkey}Genes per Cell"
        title = f"{p_xkey}Genes per Cell"
        figkey = "SFR_FIG_SS_GENE"
    ss_df = spipe.read_csv(fkey, samp=samp)

    # Shouldn't ever have more than two species
    # colors = [[0.00, 0.25, 0.8], [0.75, 0.25, 0.0], ["#2CA02C"]]

    fig = plt.figure(figsize=(6, 5), dpi=200)
    ax = fig.add_subplot(111)

    species = spipe.get_genome_list()
    # is_comb = spipe.is_combine()
    is_comb = samp.is_combine()
    # Combine any sublibraries into single lines
    for i, spec in enumerate(species):
        if is_comb:
            s_cols = [x for x in ss_df.columns if x.startswith(f"{spec}__")]
        else:
            s_cols = [x for x in ss_df.columns if x.startswith(spec)]
        # Each sublib (maybe just one)
        for j, col in enumerate(s_cols):
            lab = f"{spec}_{j+1}" if is_comb else spec
            color = subsamp_plot_color(i, j)
            spec_df = ss_df[col].dropna()
            spec_df.plot(label=lab, ax=ax, marker="o", color=color)

    # Save figure to file
    if spipe.get_par_val("rep_save_figs", as_bool=True):
        ax.legend()
        ax.set_ylabel(ylab)
        ax.set_xlabel(xlab)
        ax.set_title(title)
        fname = spipe.filepath(figkey, samp)
        fig.savefig(fname, format="png", bbox_inches="tight", dpi=200)

    clear_plot_data()

    subsamp_json = []
    for column in ss_df:
        tmp = ss_df[column].dropna(axis=0)
        subsamp_json.append(
            {
                "name": tmp.name,
                "mode": "lines",
                "x": tmp.index.to_list(),
                "y": tmp.values.tolist(),
            }
        )

    return subsamp_json


def plot_bc_rank(spipe, samp, stat_df, targeted=False, par_stat_df=None):
    """Plot tscp-per-barcode, cell cutoff, tscp-rank "knee plot"

    stat_df = dataframe with stats (i.e. analysis_summary.csv; Maybe combined)

    return json object for plotly
    """
    focal = spipe.is_focal_crispr()
    if focal:
        par_stat_dict = par_stat_df.iloc[:, 0].to_dict()

    # Get series with counts, number of these over cutoffs
    tscp_counts = spipe.read_csv("SFR_TSCP_CT", samp=samp)

    # Focal uses parent info
    if focal:
        cut_min = 1
        # cut_max = int(tscp_counts.max())
        cut_max = int(max(tscp_counts.iloc[0]))

        key = "number_of_cells"
        if key not in par_stat_dict:
            xkey = "targeted_" if targeted else ""
            key = f"{xkey}number_of_cells"

        print("focal key", key)
        n_cells = int(par_stat_dict[key])

    else:
        # Cutoff; Maybe one or combined
        cut_list = list(stat_df.loc["cell_tscp_cutoff"].values)
        cut_min = min(cut_list)
        cut_max = max(cut_list)

        # Number of cells; Previously determined (not via thresh cutoff)
        xkey = "targeted_" if targeted else ""
        row = f"{xkey}number_of_cells"
        # print("non-focal row", row)
        n_cells = int(list(stat_df.loc[row].values)[0])

    counts_df = tscp_counts["count"].sort_values(ascending=False)
    n_top = len(counts_df[counts_df >= cut_max])
    n_bot = len(counts_df[counts_df >= cut_min])

    # Static plot
    if spipe.get_par_val("rep_save_figs", as_bool=True):
        fig = plt.figure(figsize=(6, 5), dpi=200)
        ax = fig.add_subplot(111)

        # Full set (all counts_df) as grey
        ax.plot(range(len(counts_df)), counts_df.values, color="lightgray", linewidth=2)
        # If more than one (i.e. combined) thresh, middle as light green dots
        if n_top != n_bot:
            ax.plot(
                range(n_bot),
                counts_df.values[:n_bot],
                color="#A0D3AA",
                linewidth=0,
                marker=".",
            )
        # Passing cell set as green dots
        ax.plot(
            range(n_top), counts_df.values[:n_top], color="g", linewidth=0, marker="."
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        _ = ax.set_xlabel("# Barcodes (logscale)")
        if focal:
            _ = ax.set_ylabel("# Guides (logscale)")
        else:
            _ = ax.set_ylabel("# Transcripts (logscale)")

        # Annotation text
        cell_str = utils.report_num_str(n_cells)
        if cut_min != cut_max:
            tscp_str = (
                f"({utils.report_num_str(cut_min)} - {utils.report_num_str(cut_max)})"
            )
        else:
            tscp_str = f"{utils.report_num_str(cut_min)}"
        median_tscp = utils.report_num_str(counts_df[:n_cells].median())
        tx_str = f" n_cells: {cell_str}\n tscp_cutoff: {tscp_str}\n median_tscp: {median_tscp}"
        ax.text(1, 10, tx_str)
        ax.set_title("Identified Cells")

        fname = spipe.filepath("SFR_FIG_CELL_CUTOFF", samp)
        fig.savefig(fname, format="png", bbox_inches="tight", dpi=200)

    clear_plot_data()

    # Calc and collect values for plotly
    bc_x_raw = list(range(len(counts_df)))
    bc_y_raw = list(counts_df.values)

    # Subsample values in bc rank plot
    max_x = len(bc_x_raw)
    num_points = 5000
    x_step = max_x ** (1 / num_points)
    x_int = np.arange(num_points)
    bc_sub_x = np.unique((x_step**x_int).astype(int)).tolist()
    bc_sub_y = [bc_y_raw[i] for i in bc_sub_x]

    # print(f"BBB lens {len(bc_x_raw)} {len(bc_y_raw)}")
    if n_bot >= len(bc_y_raw):
        print(f"# Trimming plot bc_y_raw... (n_bot {n_bot} to {len(bc_y_raw) - 1})")
        n_bot = len(bc_y_raw) - 1

    if n_top >= len(bc_y_raw):
        print(f"# Trimming plot bc_y_raw... (n_top {n_top} to {len(bc_y_raw) - 1})")
        n_top = len(bc_y_raw) - 1

    # Insert n_top and n_bot back into bc list
    bc_sub_x.extend([n_top, n_bot])
    bc_sub_x.sort()

    bc_sub_y.extend([bc_y_raw[n_top], bc_y_raw[n_bot]])
    bc_sub_y.sort(reverse=True)

    # Convert to string for jinja (can't load int)
    bc_rank_x = list(map(str, bc_sub_x))
    bc_rank_y = list(map(str, bc_sub_y))

    bc_rank_json = [bc_rank_x, bc_rank_y, str(n_top), str(n_bot)]

    return bc_rank_json


def render_clus_gene_tab(spipe, samp):
    """Generate table with contents of table for per-cluster gene info

    Returns string content of HTML table
    """
    num_genomes = spipe.num_genomes()

    focal = spipe.is_focal_crispr()
    top_dir = None
    if focal:
        # Parent path to dataframe
        top_dir = spipe.get_parent_info(key="path")

    # Load data table
    # Default is to use first col as index, but want 'cluster' as column
    diff_exp_df = spipe.read_csv(
        "SFR_CLUST_DIFF_EXP", samp=samp, index_col=None, top_dir=top_dir
    )

    # Replace simple gene_name with url
    url_list = []
    # Interows yields idx, values
    for _, row in diff_exp_df.iterrows():
        gene_name, gene_id, genome = row.loc[["gene_name", "gene_id", "genome"]]
        # Name to report with
        if gene_name:
            rep_name = gene_name
        else:
            rep_name = gene_id
        # If multi-genome, append that
        if num_genomes > 1:
            rep_name = rep_name + "_" + genome

        # Attempt to get genome-appropriate url for gene
        url = spipe.get_gene_url(genome, gene_id)
        if url:
            url_list.append(f'<a href="{url}" target="_blank">{rep_name}</a>')
        else:
            url_list.append(f"{rep_name}")

    diff_exp_df.loc[:, "gene_name"] = pd.Series(url_list)

    """ Finalize things for (temp) per-clust gene enrichment html report
    """
    # Number formatting for score, log2_fc, pct1, pct2, and pval_adj
    diff_exp_df["score"] = diff_exp_df["score"].round(1)
    diff_exp_df["log2_FC"] = diff_exp_df["log2_FC"].round(1)
    diff_exp_df["pct1"] = diff_exp_df["pct1"].apply(lambda x: round(x * 100, 1))
    diff_exp_df["pct2"] = diff_exp_df["pct2"].apply(lambda x: round(x * 100, 1))
    diff_exp_df["pval_adj"] = diff_exp_df["pval_adj"].apply(lambda x: f"{x:.1E}")

    # Subset of cols that will be in html
    show_cols = "cluster,gene_name,score,log2_FC,pval_adj".split(",")
    diff_exp_df = diff_exp_df[show_cols]

    # Render and write html table to file
    html_out = diff_exp_df.to_html(
        classes="table",
        table_id="diff-table",
        render_links=True,
        escape=False,
        col_space=1,
        index=False,
    )
    html_out = html_out.replace("dataframe table", "table")

    # Hover over dialog for html table
    tip_rep = '<th style="min-width: 1px;" data-sortable="true" \
         data-toggle="tooltip" data-container="body" data-placement="bottom" title='

    pct1_th_old = '<th style="min-width: 1px;">pct1'
    pct1_tip = '"Percentage of cells expressing gene in cluster."'
    pct1_th_new = tip_rep + pct1_tip + ">Pct1"

    pct2_th_old = '<th style="min-width: 1px;">pct2'
    pct2_tip = '"Percentage of cells expressing gene in opposing clusters."'
    pct2_th_new = tip_rep + pct2_tip + ">Pct2"

    leiden_th_old = '<th style="min-width: 1px;">cluster'
    leiden_tip = '"Clusters generated from the Leiden algorithm."'
    leiden_th_new = tip_rep + leiden_tip + ">Cluster"

    gene_th_old = '<th style="min-width: 1px;">gene_name'
    gene_tip = '"Gene symbol taken from annotation."'
    gene_th_new = tip_rep + gene_tip + ">Gene name"

    score_th_old = '<th style="min-width: 1px;">score'
    score_tip = '"This metric selects genes with the highest log fold change \
        which are also expressed in the greatest number of cells within a cluster. \
        Formula: Log2 FC x (percentage of cells expressing gene in cluster of interest / \
            percentage of cells expressing gene all other clusters)"'
    score_th_new = tip_rep + score_tip + ">Score"

    # Non-scanpy pvalue calc
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ranksums.html
    # scipy.stats.ranksums(x, y)[source]
    #   Compute the Wilcoxon rank-sum statistic for two samples.

    pval_adj_th_old = '<th style="min-width: 1px;">pval_adj'
    pval_adj_tip = '"p-value calculated from Wilcoxon rank-sum test"'
    pval_adj_th_new = tip_rep + pval_adj_tip + ">Pval adj"

    log_FC_th_old = '<th style="min-width: 1px;">log2_FC'
    log_FC_tip = '"Log2 fold change between gene in cluster of \
        interest vs. all other clusters."'
    log_FC_th_new = tip_rep + log_FC_tip + ">Log2 FC"

    html_out = html_out.replace(pct1_th_old, pct1_th_new)
    html_out = html_out.replace(pct2_th_old, pct2_th_new)
    html_out = html_out.replace(gene_th_old, gene_th_new)
    html_out = html_out.replace(log_FC_th_old, log_FC_th_new)
    html_out = html_out.replace(leiden_th_old, leiden_th_new)
    html_out = html_out.replace(score_th_old, score_th_new)
    html_out = html_out.replace(pval_adj_th_old, pval_adj_th_new)
    html_out = html_out.replace(log_FC_th_old, log_FC_th_new)

    return html_out


def plot_barnyard(spipe, samp):
    """Plot the two-species barnyard scatter plot

    Returns nothing (barnyard not in html report)
    """
    if spipe.is_focal_crispr():
        print("Barnyard not yet for focal / crispr")
        return

    # Cell tscp count per species
    cell_data = spipe.read_csv("SFR_SPEC_TSCP_CT", samp)
    genomes = list(cell_data.columns)
    assert len(genomes) > 1, f"Barnyard needs > 1 genomes; {len(genomes)}"
    counts1 = cell_data.iloc[:, 0]
    counts2 = cell_data.iloc[:, 1]

    # Separate cells by species or mixed if not pure enough
    spurity = spipe.get_par_val("dge_species_purity_thresh", as_float=True)
    pfactor = spurity / (1 - spurity)
    cell_type1 = counts1 > (counts2 * pfactor)
    cell_type2 = counts2 > (counts1 * pfactor)
    mixed_cells = ~(cell_type1 | cell_type2)

    # Plotting per se
    fig = plt.figure(figsize=(5, 5), dpi=200)
    ax = fig.add_subplot(111)
    colors = [(0.894, 0.101, 0.109), (0.215, 0.494, 0.721), "gray"]
    s = 5
    fsize = 12
    alpha = 0.25

    plt.scatter(
        counts1[mixed_cells],
        counts2[mixed_cells],
        color=colors[2],
        s=s,
        label=None,
        alpha=alpha,
    )
    plt.scatter(
        counts1[cell_type2],
        counts2[cell_type2],
        color=colors[0],
        s=s,
        label=None,
        alpha=alpha,
    )
    plt.scatter(
        counts1[cell_type1],
        counts2[cell_type1],
        color=colors[1],
        s=s,
        label=None,
        alpha=alpha,
    )
    # These calls just put in the legend
    plt.scatter(
        [],
        [],
        color=colors[0],
        s=10,
        label="%d %s (%0.1f"
        % (sum(cell_type2), genomes[1], 100 * float(sum(cell_type2)) / len(cell_type2))
        + "%)",
    )
    plt.scatter(
        [],
        [],
        color=colors[1],
        s=10,
        label="%d %s (%0.1f"
        % (sum(cell_type1), genomes[0], 100 * float(sum(cell_type1)) / len(cell_type1))
        + "%)",
    )
    plt.scatter(
        [],
        [],
        color=colors[2],
        s=10,
        label="%d Mixed (%0.1f"
        % (sum(mixed_cells), 100 * float(sum(mixed_cells)) / len(mixed_cells))
        + "%)",
    )

    lim = int((counts1 + counts2).max() * 1.1)
    if lim > 30000:
        tickstep = 5000
    else:
        tickstep = 2000
    ax.set_xticks(np.arange(0, lim, tickstep))
    ax.set_yticks(np.arange(0, lim, tickstep))
    ax.set_xticklabels(np.arange(0, lim, tickstep), rotation=90)
    ax.axis([-int(lim / 30.0), lim, -int(lim / 30.0), lim])
    ax.set_xlabel("%s Transcript Counts" % genomes[0], fontsize=fsize)
    ax.set_ylabel("%s Transcript Counts" % genomes[1], fontsize=fsize)
    ax.tick_params(labelsize=fsize)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.legend(fontsize=fsize - 1, handletextpad=0.025, frameon=False)
    ax.get_xaxis().set_major_formatter(
        tickr.FuncFormatter(lambda x, p: format(int(x), ","))
    )
    ax.get_yaxis().set_major_formatter(
        tickr.FuncFormatter(lambda x, p: format(int(x), ","))
    )
    ax.set_title("Barnyard", size=18)

    # Save figure to file
    fname = spipe.filepath("SFR_FIG_BARNYARD", samp)
    fig.savefig(fname, format="png", bbox_inches="tight", dpi=200)

    clear_plot_data()


### -------------------------- UMAP stuff ---------------------------------


def get_umap_pl_cell_subsamp(spipe, samp, adata):
    """Get a subsample of cells (matrix rows) for umap plotting

    If subsampling param is set and number of cells is greater, subsample

    Return adata, possibly with subset of rows (i.e. cells)
    """
    bci_ss = None
    # Subsample target number
    targ_num = spipe.get_par_val("rep_umap_subsamp", as_int=True)

    # Note: adata.raw.X has counts normalized (e.g. to 10e4)
    n_cells = adata.X.shape[0]

    if (targ_num > 0) and (n_cells > targ_num):
        # Min cells per cluster for subsampling (i.e. keep at least this many per cluster)
        min_size = spipe.get_par_val("rep_umap_ss_clus_min", as_int=True)
        spipe.report_run_story2(
            f"Subsampling cells to plot; Target {targ_num} with {min_size} min per cluster"
        )
        # Get subset of indexes from list of cluster id per cell
        clust_alg = spipe.get_par_val("ana_cluster_alg", as_str=True)
        clust_list = list(adata.obs[clust_alg].values)
        cell_ss = utils.label_list_subset(clust_list, targ_num, min_keep=min_size)
        spipe.report_run_story2(
            f"Subsampled {len(cell_ss)} of {n_cells} cells for gex plots"
        )
        # Convert set of array position indexes into barcode indexes for adata
        bci_ss = set()
        for i, bc in enumerate(adata.obs.index):
            if i in cell_ss:
                bci_ss.add(bc)
        adata = adata[adata.obs.index.isin(bci_ss), :]
    else:
        spipe.report_run_story2(f"No cell subsampling for gex plots; All {n_cells}")

    # Save updated UMAP plot subset info
    if not spipe.is_focal_crispr():
        update_umap_file_cell_subsamp(spipe, samp, bci_ss)

    return adata


def get_umap_pl_par_cell_ss(spipe, samp, adata):
    """Possibly mask adata with cell subsamp info from parent

    Return (possibly filtered) adata
    """
    # Get mask from parent
    par_dir = spipe.get_parent_info(key="path")
    umap_df = spipe.read_csv("SFR_CLUST_UMAP", samp, top_dir=par_dir)

    # Barcode index subset
    bci_ss = set(umap_df[umap_df["in_plot_ss"]].index)
    adata = adata[adata.obs.index.isin(bci_ss), :]
    spipe.report_run_story2(f"Subsample umap plot cells from parent; {adata.shape}")
    return adata


def update_umap_file_cell_subsamp(spipe, samp, bci_ss, top_dir=None):
    """Update UMAP file to reflect cell subsample visibility

    Returns nothing; Write info to file
    """
    umap_df = spipe.read_csv("SFR_CLUST_UMAP", samp, top_dir=top_dir)
    # Cell subsample = set of indexes; If no set, assume all are in
    if bci_ss:
        umap_df["in_plot_ss"] = umap_df.index.isin(bci_ss)
    else:
        umap_df["in_plot_ss"] = True

    spipe.write_df(umap_df, "SFR_CLUST_UMAP", samp=samp)


def get_umap_df(spipe, samp, adata, ana_info, targeted=False):
    """Get dataframe with just UMAP plotting info (except gex counts)

    Return dataframe
    """
    focal = spipe.is_focal_crispr()
    if focal:
        pdir = spipe.get_parent_info(key="path")
        umap_df = spipe.read_csv("SFR_CLUST_UMAP", samp, top_dir=pdir)
        umap_df.rename({"umap_X": "x", "umap_Y": "y"}, axis=1, inplace=True)
    else:
        # Collect info to create plots; new dataframe mixing UMAP and cell values
        umap_df = pd.DataFrame(
            adata.obsm["X_umap"], index=adata.obs.index, columns=["x", "y"]
        )
    umap_df["sample"] = adata.obs["sample"]
    umap_df["species"] = adata.obs["species"]
    umap_df["cluster"] = get_cell_cluster_names(spipe, samp, adata, ana_info)

    return umap_df


def plot_umap(spipe, samp, umap_df, adata, ana_info, targeted=False):
    """Plot UMAP clusters

    adata should already be subsampled if needed (i.e. limit plotted cells)

    return json object for plotly
    """
    # Hardcoded parameters
    round_to = 2

    if spipe.get_par_val("rep_save_figs", as_bool=True):
        plot_umap_static(spipe, samp, umap_df)

    focal = spipe.is_focal_crispr()

    # ============ Encoded gex strings for html ==========
    if focal:
        # Only if focal specifically (not crispr)
        if spipe.is_focal_crispr(focal=True):
            # gex data strings from parent html
            lines = rep_pp.get_html_report_gex_block(spipe, samp, parent=True)
            # Includes markers and open/close paren in first/last lines
            story = f"Recycling {len(lines)} gex data lines from parent html"
            spipe.report_run_story2(story)
        else:
            lines = []

        # Zero out non-call matrix elements
        mtx = handle_mtx_zero_thresh(spipe, samp, adata.X)
        cg_mat = csc_matrix(mtx)
        gene_names = adata.var.index.values

        tlines = get_umap_gex_lines(spipe, samp, gene_names, cg_mat)
        story = f"Adding {len(tlines)} focal barcode lines"
        spipe.report_run_story2(story)
        lines += tlines

    else:
        # Convert scipy csr (compressed sparse row) to csc (c s col) matrix
        # Speeds up column slicing
        cg_mat = csc_matrix(adata.layers["norm10k"])
        gene_names = adata.var.index.values
        lines = get_umap_gex_lines(spipe, samp, gene_names, cg_mat)

    # String encoded counts
    tlines = get_umap_count_lines(
        spipe, samp, adata, targeted=targeted, round_to=round_to
    )
    lines += tlines

    # If target enrichment, add this
    if targeted:
        if focal:
            print("target enrichment UMAP not yet for focal / crispr")
        else:
            tlines = get_umap_targeted_lines(spipe, samp, umap_df, round_to=round_to)
            lines += tlines

    # X and Y coords
    # To string; Only rle encoding; no lzw and no value-to-char mapping
    x_vals = umap_df["x"].values.round(round_to)
    y_vals = umap_df["y"].values.round(round_to)
    x_str = num_list_pack_str_rle(x_vals)
    y_str = num_list_pack_str_rle(y_vals)

    # Final string forms
    x_str = f"'{x_str}'"
    y_str = f"'{y_str}'"
    cmap_str = "{\n" + "\n".join(lines) + "\n}"
    # Quoted comma,sep,list strings
    cell_str = ",".join(list(umap_df.index))
    samp_str = ",".join(list(umap_df["sample"].values))
    clust_str = ",".join(list(umap_df["cluster"].values))
    cell_str = f"'{cell_str}'"
    samp_str = f"'{samp_str}'"
    clust_str = f"'{clust_str}'"

    umap_json = {
        "umap_cells": cell_str,
        "umap_x": x_str,
        "umap_y": y_str,
        "umap_cdata": cmap_str,
        "umap_samples": samp_str,
        "umap_clusters": clust_str,
    }
    return umap_json


def plot_umap_static(spipe, samp, umap_df):
    """Plot UMAP clusters to static figures

    return nothing
    """
    # Cols and file key parts; e.g. 'READS' >--> 'SFR_FIG_UMAP_READS'
    # Actual col names may vary, so map by partial match
    k_parts = "CLUS,SAMP".split(",")
    # Only if not crispr
    if not spipe.is_focal_crispr(crispr=True):
        k_parts += "GENE,TSCP,READ".split(",")

    cols = []
    fkeys = []
    for k in k_parts:
        for c in umap_df.columns:
            if k.lower() in c.lower():
                cols.append(c)
                fkeys.append(k)
                break
    assert len(fkeys) == len(
        k_parts
    ), f"Key to col mismatch; Keys {k_parts} Cols {umap_df.columns}"

    spipe.report_run_story2(f"Generating {len(cols)} umap cluster plot figures")
    # For each plot
    for col, k in zip(cols, fkeys):
        fkey = f"SFR_FIG_UMAP_{k}"
        plot_one_umap_static(spipe, samp, fkey, umap_df, col)


def umap_plot_pt_size(spipe, n):
    """Get data-dependant plotting point size for umap plots"""
    if not isinstance(n, int):
        n = len(n)
    if n < 5000:
        outline_size = 0.2
        pt_size = 12000 / n
    else:
        outline_size = 0.1
        pt_size = 120000 / n

    return (pt_size, outline_size)


def plot_one_umap_static(spipe, samp, fkey, umap_df, val_col):
    """Plot one UMAP scatter plot

    fkey = filename key
    umap_df = data, with 'x' 'y' and 'val_col'

    Returns nothing; Saves fig to file
    """
    # Hardcoded parameters; If non-catagorical, log10
    log10 = True
    alpha = 1

    # Size params depend on number of cells (points)
    pt_size, outline_size = umap_plot_pt_size(spipe, len(umap_df))
    # Title = col name; Clean up 'mread' >--> 'read'
    title = val_col.replace("mread", "read").capitalize()
    spipe.report_run_story2(f"Generating and saving umap cluster plot {title}")

    cat = False
    # Counts are not catategorical
    if "count" not in val_col:
        cat = True
        fig = plt.figure(figsize=(5, 5), dpi=200)
        # Unique names
        u_names = sorted(umap_df[val_col].unique(), key=natsort_keygen())
        # Expand legend if too many
        lgd_cols = math.ceil(len(u_names) / 15)
        colors = sns.color_palette(n_colors=len(u_names))
        for i, val_name in enumerate(u_names):
            val_df = umap_df[umap_df[val_col] == val_name]
            plt.scatter(
                val_df["x"],
                val_df["y"],
                s=pt_size,
                alpha=alpha,
                color=colors[i],
                label=val_name,
                linewidths=outline_size,
                edgecolors="white",
            )
    else:
        # Make wider for colorbar
        fig = plt.figure(figsize=(6.3, 5), dpi=200)
        if log10:
            colors = np.log10(umap_df[val_col])
            title += " (log10)"
        else:
            colors = umap_df[val_col]
        s = plt.scatter(
            umap_df["x"],
            umap_df["y"],
            s=pt_size,
            c=colors,
            alpha=alpha,
            cmap=plt.cm.Reds,
            linewidths=outline_size,
            edgecolors="white",
        )
        _ = fig.colorbar(s)

    # gca = get current axis
    ax = plt.gca()
    ax.set(xticklabels=[], yticklabels=[])
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title(title)

    # Save figure to file
    fname = spipe.filepath(fkey, samp)
    if cat:
        lgd = ax.legend(
            bbox_to_anchor=[1.00, 0.5],
            markerscale=2,
            loc="center left",
            ncol=lgd_cols,
            frameon=False,
        )
        plt.savefig(
            fname,
            bbox_extra_artists=[lgd],
            bbox_inches="tight",
            format="png",
            pad_inches=0.25,
        )
    else:
        plt.savefig(fname, bbox_inches="tight", format="png", pad_inches=0.25)

    clear_plot_data()


def handle_mtx_zero_thresh(spipe, samp, mtx):
    """Handle possible zero-masking (sparse) matrix values less than thresh

    Return (sparse) matrix
    """
    if spipe.get_par_val("rep_umap_fb_zero_nc", as_bool=True):
        df = spipe.read_csv("SFR_TSCP_CUTOFF", samp=samp)
        thresh = int(float(df.loc["cell_tscp_cutoff", "value"]))
        if thresh > 1:
            spipe.report_run_story2(f"Masking DGE elements < thresh {thresh}")
            # This changes the matrix; < thresh elements get eliminated
            mtx = utils.mtx_zero_less_than(mtx, thresh)
        else:
            print(f">> handle_mtx_zero_thresh thresh too small {thresh}")
    return mtx


def get_umap_gex_lines(spipe, samp, gene_names, cg_mat):
    """Get (compressed) gene expression for cells

    gene_names = list of gene names
    cg_mat = data matrix

    Return list of lines (one per gene, with data for all plotted cells)
    """
    # Processing takes a while; Updating stuff
    n_genes = len(gene_names)
    update_freq = n_genes / 10
    if n_genes < 2000:
        update_freq = update_freq * 2
    update_freq = int(update_freq)

    n_cells = cg_mat.shape[0]

    spipe.report_run_story2(
        f"Packaging gene expression data for plots ({n_genes} genes, {n_cells} cells)"
    )
    nthreads = spipe.get_par_val("nthreads", as_int=True)
    nthreads = min(4, nthreads)

    # Each column = each gene
    if nthreads > 1:
        spipe.report_run_story2(f"Using {nthreads} threads for multiprocessing")
        # pool = mp.Pool(nthreads, maxtasksperchild=1000)
        pool = mp.Pool(nthreads)
        lines = pool.starmap_async(
            one_umap_gene_str,
            [(cg_mat, gene_names, g, update_freq) for g in range(len(gene_names))],
            chunksize=None,
        ).get()
        # clean up
        pool.close()
        pool.join()

    else:
        spipe.report_run_story2("Not multiprocessing")
        lines = []
        for i in range(len(gene_names)):
            num_str = one_umap_gene_str(cg_mat, gene_names, i, update_freq)
            lines.append(num_str)

    return lines


def get_umap_count_lines(spipe, samp, adata, targeted=False, round_to=2):
    """Get plot count data for cells

    Return list of lines (one per count, with data for all plotted cells)
    """
    lines = []

    # Get counts of genes, tscp, mapped reads and tscp/genes
    xkey = "targeted_" if targeted else ""

    # read_count
    counts = np.log10(adata.obs[f"{xkey}mread_count"].values + 1).round(round_to)
    num_str = num_list_pack_str_rle(counts)
    lines.append(f"'read_count':['Count', 'Log10', 'raw','{num_str}'],")

    if spipe.is_focal_crispr():
        # For focal, 'tscp_count' is gRNA_count
        counts = np.log2(adata.obs[f"{xkey}tscp_count"].values + 1).round(round_to)
        num_str = num_list_pack_str_rle(counts)
        lines.append(f"'gRNA_count':['Count', 'Log2', 'raw','{num_str}'],")

        # For focal, 'gene_count' is number of different gRNA (ndg)
        counts = np.log2(adata.obs[f"{xkey}gene_count"].values + 1).round(round_to)
        num_str = num_list_pack_str_rle(counts)
        lines.append(f"'ndg_count':['Count', 'Log2', 'raw','{num_str}'],")

    else:
        counts = np.log10(adata.obs[f"{xkey}tscp_count"].values + 1).round(round_to)
        num_str = num_list_pack_str_rle(counts)
        lines.append(f"'tscp_count':['Count', 'Log10', 'raw','{num_str}'],")

        # gene_count
        counts = np.log10(adata.obs[f"{xkey}gene_count"].values + 1).round(round_to)
        num_str = num_list_pack_str_rle(counts)
        lines.append(f"'gene_count':['Count', 'Log10', 'raw','{num_str}'],")

    return lines


def get_umap_targeted_lines(spipe, samp, umap_df, round_to=2, verb=True, prefix=""):
    """Get umap target enrichment data lines for plotting

    Return list of encoded strings (i.e. UMAP fields)
    """
    fracs_df = spipe.read_csv("SFR_ENRICHMENT", samp=samp, index_col=0, verb=verb)
    # Get (ordrerd) possible subset of fraction data (umap may be downwampled)
    fracs_df = fracs_df.loc[umap_df.index]

    lines = []
    # For each column
    for col in fracs_df.columns:
        vals = fracs_df[col].values.round(round_to)
        num_str = num_list_pack_str_rle(vals)
        lines.append(f"'{col}':['Fraction', '', 'raw','{num_str}'],")

    return lines


def one_umap_gene_str(cg_mat, gene_names, col, update_freq, round_to=2):
    """Process one gene's worth of data for umap plot

    cg_mat = sparse matrix with normalized DGE values
    gene_names = list of gene names (i.e. matrix columns)
    col = which column to process (int col index)
    update_freq = feedback print control; col Zero makes no prints

    Return str
    """
    if update_freq and ((col + 1) % update_freq) == 0:
        p_str = utils.report_percent_str(col, len(gene_names))
        print(f"# packaged {col:5d} {p_str}")
        sys.stdout.flush()

    # Cell gene counts; Slice of sparse matrix >--> array
    r_slice = np.array(np.squeeze(cg_mat.getcol(col).toarray()))
    # Scale log2
    c_slice = np.log2(r_slice + 1)
    # Round and cast to strings in array (faster than raw python)
    s_slice = np.ceil(c_slice).astype(int).astype(str)
    # Encoded data string
    num_str = int_list_pack_str_lzw(s_slice)

    return f"'{gene_names[col]}':['Expression', 'Log2', 'lzw','{num_str}'],"


##### --------------- number string encoding stuff ---------------------------

# ---- Specalized gex log2'd int case ----
# For string encoding via int-to-char mapping (used for log2'd gex values)
INTSTR2ALPH_DICT = {str(v): k for v, k in enumerate(list("abcdefghijklmnopqrstuvwxyz"))}
ALPH2INTSTR_DICT = {v: k for k, v in INTSTR2ALPH_DICT.items()}


def int_list_pack_str_lzw(nlis):
    """Create compressed string from list of int as str

    nlis = List of str for ints

    This is specialized for gex log2'd values; Separate function for speed

    Values are first run-lengh-encoded after mapping int to char via dict
    Then the rle-encoded str is further encoded via lzw compression

    Values are capped at 26 (length of ALPH2INTSTR_DICT); Above this all the same

    For example: [0,0,4,0,3,3,0,0,0,0,0,1] >--> 'a2dac3a5b' >--> lzw-encoded-str

    return string
    """
    new_lis = []
    prev = nlis[0]
    run = 1
    for v in nlis[1:]:
        if v == prev:
            run += 1
        else:
            # New value; Get dict encoding (or default)
            w = INTSTR2ALPH_DICT.get(prev, "?")
            # Add run-count if multiple row
            if run > 1:
                w += str(run)
            new_lis.append(w)
            prev = v
            run = 1
    # Last one
    w = INTSTR2ALPH_DICT.get(prev, "?")
    if run > 1:
        w += str(run)
    new_lis.append(w)

    # Final encoding
    nst = lzw_encode("".join(new_lis))

    return nst


def num_list_pack_str_rle(nlis, as_str=False, round_to=2):
    """Create compressed string from list of general numbers (as num or str)

    nlis = List of numbers or str

    For example: [0,0,4,0,30,30,0,0,0,0,0,1]  >-->  '0x2,4,0,30x2,0x5,1'

    return string
    """
    # Need list of str; Possibly round and cast
    if not as_str:
        if round_to >= 0:
            nlis = [str(round(x, round_to)) for x in nlis]
        else:
            nlis = [str(x) for x in nlis]

    # vsep and isep are separators for value-counts and items, respectively
    nst = run_len_encode(nlis, trim_zero=True, vsep="x", isep=",", v2c_dict=None)
    return nst


def run_len_encode(
    slis, trim_zero=True, vsep="x", isep=",", v2c_dict=None, v2c_def="?"
):
    """Generate run length encoded string from list of str

    slis = List of str values
    trim_zero = Flag to trim trailing zeros from (str) numbers
    vsep = Value count separator; i.e. 'x' in '2.3x4' (four instances of 2.3)
    isep = Item separator to put together output
    v2c_dict = Dict to map values to char
    v2c_def = Def char when val not found in dic

    Return string
    """
    # If not separators, make sure it's empty str
    vsep = "" if not vsep else vsep
    isep = "" if not isep else isep

    new_lis = []
    prev = slis[0]
    run = 1
    for v in slis[1:]:
        if v == prev:
            # Same value so just increment count
            run += 1
        else:
            # New value. Get run-len-encoded of prev val (str) found run times
            w = rle_encode_one_val(
                prev,
                run,
                trim_zero=trim_zero,
                vsep=vsep,
                v2c_dict=v2c_dict,
                v2c_def=v2c_def,
            )
            new_lis.append(w)
            prev = v
            run = 1
    # Last one
    w = rle_encode_one_val(
        prev, run, trim_zero=trim_zero, vsep=vsep, v2c_dict=v2c_dict, v2c_def=v2c_def
    )
    new_lis.append(w)

    return isep.join(new_lis)


def rle_encode_one_val(w, run, trim_zero=False, vsep="x", v2c_dict=None, v2c_def="?"):
    """Get run length encoding string for one value and run count"""
    if trim_zero:
        w = num_trim_zero_str(w)
    if v2c_dict:
        try:
            w = v2c_dict[w]
        except Exception:
            w = v2c_def
    if run > 1:
        w += vsep + str(run)
    return w


# Regex to strip trailing zeros after decimal; For example
#   2.30 >--> 2.3
#   0.0 >--> 0
#   100 >--> 100
#   from https://stackoverflow.com/questions/44111169/remove-trailing-zeros-after-the-decimal-point-in-python
TRIM_TRAILZERO = re.compile(r"(?:(\.)|(\.\d*?[1-9]\d*?))0+(?=\b|[^0-9])")


def num_trim_zero_str(val):
    """Trim trailing zeros from number as STRING

    Return string
    """
    return TRIM_TRAILZERO.sub(r"\2", val)


def lzw_encode(s):
    """This function compresses a string using the lzw algorithm

    Return string
    """
    dict_encode = {}
    data = list(s)
    out = []
    phrase = data[0]
    code = 256

    for i in range(1, len(data)):
        currChar = data[i]
        try:
            dict_encode["_" + phrase + currChar]
        except Exception:
            if len(phrase) > 1:
                out.append(dict_encode["_" + phrase])
            else:
                out.append(ord(phrase[0]))
            dict_encode["_" + phrase + currChar] = code
            code += 1
            phrase = currChar
        else:
            phrase += currChar
    if len(phrase) > 1:
        out.append(dict_encode["_" + phrase])
    else:
        out.append(ord(phrase[0]))

    for i in range(0, len(out)):
        out[i] = chr(out[i])
    return "".join(out)


def clean_up(spipe, final=False):
    """Clean up report stuff for given sample"""
    # Clean up any open figure
    plt.close("all")

    if spipe.get_par_val("clean_per_step", as_bool=True) or final:
        # Temp files to remove?
        f_lis = []
        utils.check_and_rm_file(f_lis)
