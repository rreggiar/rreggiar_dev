#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import sys
import traceback
import numpy as np
import pandas as pd
import scipy
import scanpy
from collections import defaultdict
from statsmodels.stats.multitest import multipletests


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import pb_const as const
    import pb_irlbpy as irlbpy
    import utils
else:
    import splitpipe.pb_const as const
    import splitpipe.pb_irlbpy as irlbpy
    from splitpipe import utils


# ---------------------------------------------------------------------------
def have_analysis_ins(spipe, verb=True, samp=None):
    """Check if analysis inputs exist

    Return True/False, list of files
    """
    ok = False
    check_files = []
    dge_fkeys = const.out_filepaths(spipe, samp, use_case=True, as_key=True)
    # Files for given or all samples
    if samp:
        check_files.append(spipe.filepath(dge_fkeys["mtx"], samp))
        check_files.append(spipe.filepath("SFR_ASUM_CSV", samp))
    else:
        for samp in spipe.get_samples():
            check_files.append(spipe.filepath(dge_fkeys["mtx"], samp))
            check_files.append(spipe.filepath("SFR_ASUM_CSV", samp))

    # If list is good, we have inputs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def have_analysis_outs(spipe, verb=False, samp=None):
    """Check if analysis output files exist

    Return True/False, list of files
    """
    ok = False
    check_files = []
    # Analysis processing files
    if samp:
        check_files.append(spipe.filepath("SFR_ANA_PROC_DEF", samp))
    else:
        for samp in spipe.get_samples():
            check_files.append(spipe.filepath("SFR_ANA_PROC_DEF", samp))
    # If list is good, we have outs
    if utils.check_infile(check_files, verb=verb):
        ok = True
    return ok, check_files


def run_analysis(spipe):
    """Run analysis step

    Return status
    """
    spipe.report_proc_step("Analysis")

    # Input and output status and list
    # If keep_going, input checks are per sample
    if spipe.keep_going():
        i_ok = True
    else:
        i_ok, i_list = have_analysis_ins(spipe)
    o_ok, o_list = have_analysis_outs(spipe)
    # Run only if don't have output or fresh
    if (not o_ok) or spipe.force_fresh_files():
        if not i_ok:
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for analysis: {bad_list}"
            spipe.set_problem(story)
        else:
            generate_all_analyses(spipe)
    else:
        spipe.report_run_story("Skipping Analysis and using existing outputs")
        spipe.report_run_story2(f"Outputs: {o_list}")
    ok = spipe.no_problems()
    spipe.report_proc_step("Analysis", status=ok)
    return ok


def generate_all_analyses(spipe):
    """Wrapper function to loop through every sample and do analysis

    Return nothing; If any problems, these are reported / remembered
    """
    # Each sample
    samp_list = spipe.get_samples()
    for i, samp in enumerate(samp_list):
        ok = True
        spipe.report_run_story("")
        story = (
            f"Analysis processing sample ({i+1} of {len(samp_list)}) {samp.get_name()}"
        )
        spipe.report_run_story(story)

        # If already not good, skip
        if not samp.get_status(is_good=True):
            stat_step = samp.get_status(step=True)
            story = f"Status of sample {samp.get_name()} bad ({stat_step}); Skipping"
            spipe.report_run_story(story)
            continue

        spipe.report_run_mem(story=f"Analysis samp {i+1}")

        # Have required inputs?
        i_ok, i_list = have_analysis_ins(spipe, samp=samp)
        if not i_ok:
            ok = False
            bad_list = utils.bad_infile_list(i_list)
            story = f"Don't have inputs for analysis: {bad_list}"
            spipe.set_problem(story)
        else:
            # Catch possible problems (e.g. scanpy ...)
            try:
                ok = generate_single_analysis(spipe, samp)
            except Exception as e:
                ok = False
                # Dump traceback, stdout and log
                traceback.print_exc(file=sys.stdout)
                sys.stdout.flush()
                traceback.print_exc(file=spipe.get_log_fh())
                story = f"Failed generate_single_analysis: {e}"
                spipe.set_problem(story)

        # Problem, possibly bail
        if not ok:
            story = f"Analysis problem with sample {samp.get_name()}"
            spipe.set_problem(story)
            samp.set_status(False, "analysis")
            if spipe.keep_going():
                spipe.report_run_story("Keep-going True... So, keep on going!")
                spipe.clear_problems()
                print()
            else:
                break


def generate_single_analysis(spipe, samp):
    """Run analysis for given sample

    Return status
    """
    dge_paths = const.out_filepaths(spipe, samp)

    # Load cell-gene data into anndata object
    # Make sure anndata is constructed fresh from component parts (not from file)
    adata = utils.load_parsebio_dge(dge_paths["dir"], anndata_fresh=True)

    if adata is None:
        story = f"No DGE for {samp.get_name()} ({dge_paths['dir']})"
        spipe.report_warning(story)
        return False

    spipe.report_run_story2(f"Loaded anndata {adata.shape} {dge_paths['dir']}")

    # Cluster status flag var; Assume all good, change if not
    clust_status = "good"

    n_cells, n_genes = adata.shape
    # No clustering for focal barcoding or crispr
    if spipe.is_focal_crispr():
        clust_status = "focal_bc"
        spipe.report_run_story("Focal barcoding; Not clustering")
    else:
        # Only cluster if enough cells (barcode) and genes
        min_bc = spipe.get_par_val("ana_cluster_min_bc", as_int=True)
        if (n_cells >= min_bc) and (n_genes >= min_bc):
            spipe.report_run_story2(
                f"Calling Scanpy to cluster (adata shape {adata.shape})"
            )
            adata, clust_ok = scanpy_cluster_anndata(spipe, adata)
            if not clust_ok:
                clust_status = "failed_clustering"
                story = f"Failed to calculate clusters: {adata.shape}"
                spipe.report_run_story(story)
        else:
            clust_status = "too_few_cells_or_genes"
            story = f"Too few cells / genes to cluster: {adata.shape} < {min_bc}"
            spipe.report_run_story(story)

    # Save anndata to file
    adata.write(dge_paths["anndata"])
    spipe.report_run_story2(f"Wrote anndata: {dge_paths['anndata']}")

    # Get gene enrichment table if cluster data is ok
    if clust_status == "good":
        # Sample cluster assignment
        write_samp_clus_files(spipe, samp, adata)
        # Cluster gene enrichment table
        n_clust = make_clus_gene_tab(spipe, samp, adata)
    else:
        n_clust = 0

    # Save processing details
    write_samp_ana_proc(spipe, samp, clust_status, n_clust, adata.shape)

    return True


def write_samp_ana_proc(spipe, samp, status, n_clust, ad_shape):
    """Write out analysis processing info (for report logic)"""
    fkey = "SFR_ANA_PROC_DEF"
    o_dict = {}
    o_dict["header"] = spipe.get_header_dict(desc="analysis processing", fkey=fkey)
    o_dict["analysis"] = {
        "status": status,
        "cluster_num": n_clust,
        "cluster_alg": spipe.get_par_val("ana_cluster_alg", as_str=True),
        "cluster_min_bc": spipe.get_par_val("ana_cluster_min_bc", as_int=True),
        "cluster_1base": spipe.get_par_val("ana_cluster_1base", as_bool=True),
        "anndata_rows": ad_shape[0],
        "anndata_cols": ad_shape[1],
    }
    fpath = spipe.filepath("SFR_ANA_PROC_DEF", samp)
    utils.write_json(o_dict, fpath)


def write_samp_clus_diff_exp(spipe, samp, diff_exp):
    """Write out cluster diff expression table"""
    # If 1-based cluster numbering, bump if smallest is zero
    if spipe.get_par_val("ana_cluster_1base", as_bool=True):
        if diff_exp["cluster"].min() == 0:
            diff_exp["cluster"] += 1
    # Limit precision for number cols
    cols = "score,log2_FC,pct1,pct2".split(",")
    diff_exp = diff_exp.round(dict(zip(cols, [3] * len(cols))))
    # For scientific notation, need format (not round)
    diff_exp["pval_adj"] = diff_exp["pval_adj"].apply(lambda x: f"{x:0.3e}")
    spipe.write_df(diff_exp, "SFR_CLUST_DIFF_EXP", samp=samp, index=False)


def write_samp_clus_files(spipe, samp, adata):
    """Write cluster info to file(s)"""
    # Cluster assignment
    clust_alg = spipe.get_par_val("ana_cluster_alg", as_str=True)
    clusters = adata.obs[clust_alg].astype(int)
    # If 1-based cluster numbering, bump if smallest is zero
    if spipe.get_par_val("ana_cluster_1base", as_bool=True):
        if clusters.min() == 0:
            clusters += 1
    # Series with data name as cluster alg (e.g. 'leiden'); Set to 'cluster'
    spipe.write_df(
        clusters.rename("cluster"), "SFR_CLUST_ASSIGN", samp=samp, index=True
    )

    # Cluster UMAP (conditional)
    if spipe.get_par_val("ana_save_umap", as_bool=True):
        cols = "umap_X,umap_Y".split(",")
        umap_df = pd.DataFrame(
            adata.obsm["X_umap"], index=adata.obs.index, columns=cols
        )
        umap_df["umap_X"] = umap_df["umap_X"].round(4)
        umap_df["umap_Y"] = umap_df["umap_Y"].round(4)
        spipe.write_df(umap_df, "SFR_CLUST_UMAP", samp=samp, index=True)


# ===================== Scanpy processing code (filter, cluster, umap) ===============


def dump_adata(spipe, adata, story):
    # pass
    # what = adata
    story = f"SCA {adata.shape} {story}"
    spipe.report_run_mem(story=story, verb=2)


def scanpy_cluster_anndata(spipe, adata):
    """Run scanpy processing for clustering

    Return anndata and processing status flag
    """
    # Set vars from parameters
    min_cells = spipe.get_par_val("ana_gene_cell_min", as_int=True)
    rseed = spipe.get_par_val("rseed", as_int=True)
    pca_min = spipe.get_par_val("ana_pca_dim_min", as_int=True)
    pca_max = spipe.get_par_val("ana_pca_dim_max", as_int=True)
    zero_center = spipe.get_par_val("ana_scale_zero_center", as_bool=True)
    use_irlbpy = spipe.get_par_val("ana_pca_irlbpy", as_bool=True)
    scale_max_value = spipe.get_par_val("ana_scale_max_value", as_float=True)
    norm_target_sum = spipe.get_par_val("ana_norm_target_sum", as_float=True)

    # Targeted enrichment special case
    targeted = spipe.is_targeted()
    if targeted:
        pca_max = spipe.get_par_val("ana_pca_dim_targeted_max", as_int=True)
        story = f"Target enrichment; PCA max dim set to {pca_max}"
        spipe.report_run_story2(story)

    # Run through processing steps
    dump_adata(spipe, adata, f"function scanpy_cluster_anndata start, rseed={rseed}")

    scanpy.pp.filter_genes(adata, min_cells=min_cells)
    dump_adata(spipe, adata, f"After filter_genes {min_cells}")

    # Depreciated: https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.normalize_per_cell.html#scanpy-pp-normalize-per-cell
    # scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    # dump_adata(spipe, adata, "After normalize_per_cell")
    # Replace with: https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.normalize_total.html#scanpy.pp.normalize_total
    scanpy.pp.normalize_total(adata, target_sum=norm_target_sum)
    dump_adata(spipe, adata, f"After normalize_total target_sum={norm_target_sum}")

    # Note that we are storing normalized (but not log transformed) counts as raw
    # This will break normal scanpy differential expression (it expects log transformed)
    # Use raw (normalize_total) for cluster gene enrichment table
    # adata.raw = adata
    # dump_adata(spipe, adata, "After setting raw (to normalized_total)")
    adata.layers["norm10k"] = adata.X

    scanpy.pp.log1p(adata)
    dump_adata(spipe, adata, "After log1p")

    # USe defaults (e.g. min_disp=0.5, min_mean=0.0125, max_mean=3)
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html
    scanpy.pp.highly_variable_genes(adata)
    dump_adata(spipe, adata, "After highly_variable_genes")

    # If irlb instead of PCA, the centering is done later (not here)
    if use_irlbpy:
        zero_center = False

    scanpy.pp.scale(adata, max_value=scale_max_value, zero_center=zero_center)
    dump_adata(
        spipe, adata, f"After scale (max={scale_max_value}, zero_center={zero_center})"
    )

    # Number of highly variable
    min_hvar = spipe.get_par_val("ana_pca_min_hvar", as_int=True)
    num_hva = adata.var["highly_variable"].sum()
    if num_hva >= min_hvar:
        use_hv = True
        spipe.report_run_story2(f"PCA using {num_hva} highly variable genes")
    else:
        use_hv = False
        spipe.report_run_story2(
            f"PCA using all genes; Not just highly variable {num_hva}"
        )

    if targeted:
        use_hv = False
        spipe.report_run_story2(
            f"PCA using all genes; Targeted so not highly variable {num_hva}"
        )

    # PCA component dimension bounds
    n_comps = min(pca_max, adata.shape[0] - 1, adata.shape[1] - 1)
    if n_comps < pca_min:
        story = f"PCA dim issue; {n_comps} vs pca_min {pca_min} max {pca_max}, adata {adata.shape}"
        spipe.report_run_story2(story)
        return None, False

    # irlbpy or normal PCA
    if use_irlbpy:
        spipe.report_run_story2(f"PCA irlbpy, n_comps {n_comps}")
        # irlb with zero centering in wrapper (not in scanpy.scale())
        anndata_pca_via_irlbpy(spipe, adata, n_comps, use_hv, zero_center=True)
        dump_adata(spipe, adata, "After irlbpy; Next neighbors")
    else:
        spipe.report_run_story2(f"PCA default, non-irlbpy, n_comps {n_comps}")
        # Default 'solver' is 'arpack' which should pay attention to random seed
        scanpy.tl.pca(
            adata, n_comps=n_comps, use_highly_variable=use_hv, random_state=rseed
        )
        dump_adata(spipe, adata, "After pca; Next neighbors")

    # https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html
    # Default n_neighbors=15
    n_neighbors = spipe.get_par_val("ana_umap_n_neighbors", as_int=True)
    if n_neighbors > n_comps:
        n_neighbors = n_comps

    scanpy.pp.neighbors(
        adata, n_neighbors=n_neighbors, n_pcs=n_comps, random_state=rseed
    )
    dump_adata(spipe, adata, f"After neighbors, n_neighbors={n_neighbors}; Next umap")

    # https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html#scanpy.tl.umap
    # n_components=2 = 2D plot
    # min_dist = minimum dist between points; Default=0.5 (umap-learn def=0.1)
    min_dist = spipe.get_par_val("ana_umap_min_dist", as_float=True)
    min_dist = 0.5
    # spread = "effective scale of embedded points"; Default=1
    spread = 1
    scanpy.tl.umap(
        adata, n_components=2, min_dist=min_dist, spread=spread, random_state=rseed
    )
    dump_adata(spipe, adata, f"After umap, min_dist={min_dist}")

    # Clustering
    clust_alg = spipe.get_par_val("ana_cluster_alg", as_str=True).lower()
    spipe.report_run_story2(f"Clustering ({clust_alg})")
    if clust_alg == "leiden":
        scanpy.tl.leiden(adata, random_state=rseed)
    elif clust_alg == "louvain":
        scanpy.tl.louvain(adata, random_state=rseed)
    else:
        story = f"Analysis unknown (scanpy) cluster algorithm {clust_alg}"
        spipe.set_problem(story)
        return None, False

    dump_adata(spipe, adata, f"After clustering {clust_alg}, done")
    return adata, True


def anndata_pca_via_irlbpy(spipe, adata, n_comps, use_hv, zero_center=True):
    """Call irlbpy code as substitute for anndata PCA call

    adata = anndata, processed up to here (e.g. log, "highly variable")
    n_comps = n PCA components
    use_hv = Flag to use "highly_variable" genes only

    """
    # Get starting matrix; highly variable subset of genes or all
    if use_hv:
        A = adata.X[:, adata.var["highly_variable"]]
    else:
        A = adata.X

    story = f"PCA via irlbpy, input shape {A.shape}, n_comps={n_comps}, use_hv={use_hv}, zero_center={zero_center}"
    spipe.report_run_story2(story)

    """
    Function to call:
        (from https://github.com/airysen/irlbpy/blob/master/irlbpy/irlb.py)

        S = lanczos(A, nval, tol=0.0001, maxit=50, center=None, scale=None, L=None):
    """
    # Increase stringency vs default params (above)
    # maxit = 50
    maxit = 200
    # tol=0.0001
    tol = 0.00001

    if zero_center:
        means = A.mean(0)
        S = irlbpy.lanczos(A, n_comps, tol=tol, maxit=maxit, center=means)
    else:
        S = irlbpy.lanczos(A, n_comps, tol=tol, maxit=maxit)

    # Unpack results into anndata (as if native pca output)
    """
    From irlbpy/irlb.py

        Returns:
            X[0] A j * nu matrix of estimated left singular vectors.
            X[1] A vector of length nu of estimated singular values.
            X[2] A k * nu matrix of estimated right singular vectors.
            X[3] The number of Lanczos iterations run.
            X[4] The number of matrix-vector products run.

        In object:
        return LanczosResult(**{'U': U,
                                's': S[1][0:nu],
                                'V': V,
                                'steps': it,
                                'nmult': mprod
                                })

    From scanpy/_pca.py

        adata.obsm['X_pca'] = X_pca
        adata.uns['pca'] = {}
        adata.uns['pca']['params'] = {
            'zero_center': zero_center,
            'use_highly_variable': use_highly_variable,
        }
        if use_highly_variable:
            adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm['PCs'][adata.var['highly_variable']] = pca_.components_.T
        else:
            adata.varm['PCs'] = pca_.components_.T
        adata.uns['pca']['variance'] = pca_.explained_variance_
        adata.uns['pca']['variance_ratio'] = pca_.explained_variance_ratio_

    """
    # Need to scale; If set X_pca to only S.U, clusters collapse / disperse some
    adata.obsm["X_pca"] = S.U * S.s
    adata.uns["pca"] = {}
    adata.uns["pca"]["params"] = {
        "zero_center": False,
        "use_highly_variable": use_hv,
    }
    if use_hv:
        adata.varm["PCs"] = np.zeros(shape=(adata.n_vars, n_comps))
        adata.varm["PCs"][adata.var["highly_variable"]] = S.V
    else:
        adata.varm["PCs"] = S.V
    # Calc variance explicitly
    adata.uns["pca"]["variance"] = adata.obsm["X_pca"].var(0)
    # NOTE; This requires scaling (of starting matrix) to be correct
    adata.uns["pca"]["variance_ratio"] = adata.uns["pca"]["variance"] / (
        adata.X.shape[1] - 1
    )

    return True


def make_clus_gene_tab(spipe, samp, adata):
    """Generate table with contents of table for per-cluster gene info

    Non-scanpy version of gene enrichment calcs

    adata should have 'norm10k' layer and clustering info

    Output saved to file; Return number of clusters
    """
    clust_alg = spipe.get_par_val("ana_cluster_alg", as_str=True)
    # number of top-scoring genes to keep per cluster
    n_top = spipe.get_par_val("ana_ctab_top_n", as_int=True)
    # Fraction (in / non) percent cells denominator minimum
    fc_denom_min = spipe.get_par_val("ana_ctab_sco_min_den", as_float=True)

    # Calculate differential expression and fraction expressed
    #   for all genes (remaining in adata) across all clusters
    clusters = adata.obs[clust_alg].apply(int)
    # Raw (sparse) matrix (non-scaled counts) and list of genes
    # mtx = adata.raw.X
    # Use gene IDs as names maybe not unique (Series with val = gene_id)
    # gene_ids = adata.raw.var['gene_id']

    mtx = adata.layers["norm10k"]
    gene_ids = adata.var["gene_id"]

    # Create df to hold to-be-reported data indexed by id (for later)
    gene_data = pd.DataFrame(index=list(gene_ids))
    # Copy .values to get strings (else NaN)
    gene_data["gene_id"] = adata.var["gene_id"].values
    gene_data["gene_name"] = adata.var["gene_name"].values
    gene_data["genome"] = adata.var["genome"].values

    # For data collection; Columns = cluster
    cluster_fraction = pd.DataFrame()
    not_cluster_fraction = pd.DataFrame()
    cluster_diff = pd.DataFrame()
    gene_rank_score = pd.DataFrame()

    sorted_clust_list = sorted(clusters.unique())
    story = f"Getting gene enrichment table, {len(sorted_clust_list)} clusters, top {n_top} genes each"
    spipe.report_run_story(story)
    for i in sorted_clust_list:
        spipe.report_run_story3(f"Collecting, cluster {i}")

        # Masks to separate cells
        is_cluster = clusters == i
        not_is_cluster = clusters != i

        # Get cells in the current cluster
        clust_slice = mtx[is_cluster, :]
        # The 'A1' unwraps numpy matrix to simple array
        #   https://numpy.org/doc/stable/reference/generated/numpy.matrix.A1.html
        clust_exp = pd.Series(index=gene_ids, data=clust_slice.sum(0).A1)
        cluster_fraction[i] = pd.Series(
            index=gene_ids, data=(clust_slice > 0).mean(0).A1
        )
        # Convert to transcripts per million (and add 1 pseudocount)
        clust_exp = clust_exp / clust_exp.sum() * 1e6 + 1

        # Get all cells NOT in the current cluster
        not_clust_slice = mtx[not_is_cluster, :]
        not_clust_exp = pd.Series(index=gene_ids, data=not_clust_slice.sum(0).A1)
        not_cluster_fraction[i] = pd.Series(
            index=gene_ids, data=(not_clust_slice > 0).mean(0).A1
        )
        # Convert to transcripts per million (and add 1 pseudocount)
        not_clust_exp = not_clust_exp / not_clust_exp.sum() * 1e6 + 1

        # Get log2 differential expression (cluster vs all other clusters)
        cluster_diff[i] = np.log2(clust_exp / not_clust_exp)

        # Rank score = fold_change * (percent_in / percent_out)
        # Limit how small denominator can be; Values change but rank shouldn't (much?)
        gene_rank_score[i] = (
            cluster_diff[i]
            * cluster_fraction[i]
            / not_cluster_fraction[i].clip(lower=fc_denom_min)
        )


    # This picks the genes with the highest fold-change-score per cluster.
    # This typically works pretty well because the pseudocount
    # forces high average expression in the cluster of interest
    # or else you won't get a high fold change compared to all
    # other clusters.

    # Collect lists of values for top genes of each cluster
    diff_exp = defaultdict(list)
    for i in cluster_diff.columns:
        spipe.report_run_story3(f"Calculating, cluster {i}")

        # Sort on gene rank score, descending for this cluster
        current_rank_scores = gene_rank_score[i].sort_values(ascending=False)
        # Cluster subsets for other data
        current_cluster_diff = cluster_diff[i]
        current_cluster_frac = cluster_fraction[i]
        not_cur_cluster_frac = not_cluster_fraction[i]

        # Masks to separate cells
        is_cluster = clusters == i
        not_is_cluster = clusters != i

        for j in range(n_top):
            # Gene from sorted value index, score at sorted position
            g = current_rank_scores.index[j]
            score = current_rank_scores.iloc[j]

            diff_exp["cluster"].append(i)
            diff_exp["gene_id"].append(gene_data.loc[g, "gene_id"])
            diff_exp["gene_name"].append(gene_data.loc[g, "gene_name"])
            diff_exp["genome"].append(gene_data.loc[g, "genome"])
            diff_exp["score"].append(score)
            diff_exp["log2_FC"].append(current_cluster_diff.loc[g])
            diff_exp["pct1"].append(current_cluster_frac.loc[g])
            diff_exp["pct2"].append(not_cur_cluster_frac.loc[g])

            # Get distributions for gene g for all cells in / out cluster
            # The mtx indexing returns a matrix; Use '.A1' to unwrap into simple array
            exp_cluster = mtx[is_cluster, gene_ids == g].A1
            exp_not_cluster = mtx[not_is_cluster, gene_ids == g].A1

            # Two arrays; No length constraint; whatever in/out of current cluster
            #   https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ranksums.html
            pstat = scipy.stats.ranksums(exp_cluster, exp_not_cluster)[1]
            diff_exp["pval_adj"].append(pstat)

    # Dataframe from collected parts
    diff_exp = pd.DataFrame(diff_exp)

    # Calculate adjusted p-values using bonferroni correction:
    diff_exp["pval_adj"] = multipletests(
        diff_exp.pval_adj, alpha=0.05, method="bonferroni"
    )[1]

    # Save to file here; Rendering in reports.py
    write_samp_clus_diff_exp(spipe, samp, diff_exp)

    # Return number of clusters
    # return len(cluster_diff.columns), cluster_diff, diff_exp

    return len(cluster_diff.columns)
