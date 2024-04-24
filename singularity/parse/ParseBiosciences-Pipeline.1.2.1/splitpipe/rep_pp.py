#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Report (HTML) post processing file manipulation
#


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


### ------------------- HTML file post-proc -------------------


def post_proc_html_report(spipe, samp, ofname=""):
    """Possibly post process new html report (writing a new one)

    Return status
    """
    new_lines = []

    if spipe.is_tcr():
        new_lines = post_proc_tcr_html(spipe, samp)
    elif spipe.is_focal_crispr():
        new_lines = post_proc_focal_html(spipe, samp)
    else:
        return True

    if not new_lines:
        spipe.set_problem(f"Post-process html merge for sample {samp.get_name()}")
        return False

    # If keeping temp files, save default html (may not exist)
    if spipe.keep_temp_files():
        ofname = spipe.filepath("SF_ASUM_HTML", samp)
        dst = spipe.filepath("TMP_SF_ASUM_HTML", samp)
        utils.copy_file(ofname, dst, verb=False)

    # Save
    if not ofname:
        ofname = spipe.filepath("SF_ASUM_HTML", samp)
    with open(ofname, "w") as OUTFILE:
        print("\n".join(new_lines), file=OUTFILE)
        spipe.report_run_story2(f"Post processed html: {ofname}")

    return True


def post_proc_focal_html(spipe, samp):
    """Post process new and parent html files together

    Return list of lines
    """
    spipe.report_run_story2("Post processing html (focal)")
    # Current (raw focal / crispr) and parent html as list of lines
    cur_lines = get_html_report_lines(spipe, samp, parent=False)
    par_lines = get_html_report_lines(spipe, samp, parent=True)

    # Parent: Want almost all, including header, pages 1 and 2; Not footer
    # Do not want UMAP cmap2, so get up and downstream parts around this
    p_up_lines = line_list_extract_block(par_lines, e_key="block_end_umap_cmap1")
    # Don't actually want first line (with start key)
    p_dn_lines = line_list_extract_block(
        par_lines, s_key="block_end_umap_cmap2", e_key="block_end_page_Page2"
    )[1:]
    # Parent (header) lines have javascript state var assignments
    # Set focal / crispr variable true for output
    if spipe.is_focal_crispr(focal=True):
        p_up_lines = simple_sub_list_lines(
            p_up_lines, "const is_focal = false", "const is_focal = true"
        )
    else:
        p_up_lines = simple_sub_list_lines(
            p_up_lines, "const is_crispr = false", "const is_crispr = true"
        )

    # Want all of focal page1; Will be page3 in output
    c_pg3_lines = line_list_extract_block(
        cur_lines, s_key="block_start_page_Page1", stoe=True
    )
    # UMAP colormap; Starts as cmap1 but will be cmap2 in output
    c_cmap2_lines = line_list_extract_block(
        cur_lines, s_key="block_start_umap_cmap1", stoe=True
    )

    # Page 1 becomes page 3
    c_pg3_lines = simple_sub_list_lines(c_pg3_lines, "Page1", "Page3")
    c_pg3_lines = simple_sub_list_lines(
        c_pg3_lines, "pg1_sum_title_line", "pg3_sum_title_line"
    )
    c_pg3_lines = simple_sub_list_lines(c_pg3_lines, "pg1-plotly", "pg3-plotly")
    # Color map 1 becomes color map 2 (Var change but also update landmark comments)
    c_cmap2_lines = simple_sub_list_lines(
        c_cmap2_lines, "const umap_cmap_dict1", "const umap_cmap_dict2"
    )
    c_cmap2_lines = simple_sub_list_lines(c_cmap2_lines, "_umap_cmap1", "_umap_cmap2")

    # For crispr (not focal) remove saturation and bc rank plots
    if spipe.is_focal_crispr(crispr=True):
        # Inverse removes blocks
        # Javascript (plotly) for bc-rank, tscp sat and gene sat plots are grouped
        c_pg3_lines = line_list_extract_block(
            c_pg3_lines, s_key="block_start_js_bc_rank_plot", stoe=True, inverse=True
        )
        c_pg3_lines = line_list_extract_block(
            c_pg3_lines, s_key="block_start_js_tscp_sat_plot", stoe=True, inverse=True
        )
        c_pg3_lines = line_list_extract_block(
            c_pg3_lines, s_key="block_start_js_gene_sat_plot", stoe=True, inverse=True
        )
        # Div sections of report
        c_pg3_lines = line_list_extract_block(
            c_pg3_lines, s_key="block_start_div_id_bc_rank_box", stoe=True, inverse=True
        )
        c_pg3_lines = line_list_extract_block(
            c_pg3_lines,
            s_key="block_start_div_id_tscp_sat_box",
            stoe=True,
            inverse=True,
        )
        c_pg3_lines = line_list_extract_block(
            c_pg3_lines,
            s_key="block_start_div_id_gene_sat_box",
            stoe=True,
            inverse=True,
        )

    # Need end of body and footer
    foot_lines = ["</body>"] + line_list_extract_block(par_lines, s_key="<footer>")

    # Build output lines from parts
    new_lines = p_up_lines + c_cmap2_lines + p_dn_lines + c_pg3_lines + foot_lines

    # Update nav blocks to include Page3
    new_lines = update_nav_list_lines(new_lines, "Summary", "Clustering", "CRISPR")

    return new_lines


def post_proc_tcr_html(spipe, samp):
    """Post process merging of new html and parent for TCR

    Returns list of lines for new html
    """
    spipe.report_run_story2("Post processing html (TCR)")
    # Current html as list of lines
    cur_lines = get_html_report_lines(spipe, samp, parent=False)
    # Have parent
    if spipe.have_parent():
        par_lines = get_html_report_lines(spipe, samp, parent=True)

        # Want all of parent html (e.g. header) up until end of page2
        par_lines = line_list_extract_block(par_lines, e_key="block_end_page_Page2")
        # Current TCR want only Page3 to end (no header)
        cur_lines = line_list_extract_block(cur_lines, s_key="block_start_page_Page3")

        # Set TCR flags
        par_lines = simple_sub_list_lines(
            par_lines, "const is_tcr_parent = false", "const is_tcr_parent = true"
        )
        par_lines = simple_sub_list_lines(
            par_lines, "const is_tcr_only = true", "const is_tcr_only = false"
        )

        new_lines = par_lines + cur_lines

        # Update nav blocks to include TCR Page3
        new_lines = update_nav_list_lines(new_lines, "Summary", "Clustering", "TCR")

    # No parent
    else:
        # Update with added TCR only; Other links gone
        new_lines = update_nav_list_lines(cur_lines, "", "", "TCR Summary")

    return new_lines


def get_html_report_lines(spipe, samp, parent=False, verb=True):
    """Load HTML report file into list of lines

    Return list of lines
    """
    top_dir = None
    if parent:
        top_dir = spipe.get_parent_info(key="path")
    fname = spipe.filepath("SF_ASUM_HTML", samp, top_dir=top_dir)
    lines = utils.lines_from_fname(fname, comment=None)
    if verb:
        spipe.report_run_story2(f"Loaded {len(lines)} lines from {fname}")
    return lines


### ------------------- HTML file manipulations -------------------


def line_list_block_coords(lines, s_key=None, e_key=None, inclusive=True, stoe=False):
    """Search list of lines for block between start, end keys

    inclusive = flag to include lines with key blocks
    stoe = flag to set end key via start key (e.g. 'start_xyz' >--> 'end_xyz')

    Return tuple(start, end)
    """
    if stoe:
        assert isinstance(s_key, str), f"stoe true but not s_key {type(s_key)}"
        e_key = s_key.replace("_start_", "_end_")

    b_st = b_en = -1
    for i, line in enumerate(lines):
        if s_key and (b_st < 0) and s_key in line:
            b_st = i
            if not inclusive:
                b_st += 1
        if e_key and (b_en < 0) and e_key in line:
            b_en = i
            if inclusive:
                b_en += 1

    return (b_st, b_en)


def line_list_extract_block(
    lines, s_key=None, e_key=None, inclusive=True, inverse=False, stoe=False
):
    """Search list of lines for block between start, end keys

    inclusive = flag to keep lines with key blocks
    inverse = flag go invert what's saved; If inverse=True, block removed
    stoe = flag to set end key via start key (e.g. 'start_xyz' >--> 'end_xyz')

    Return list of lines
    """
    if stoe:
        assert isinstance(s_key, str), f"stoe true but not s_key {type(s_key)}"
        e_key = s_key.replace("_start_", "_end_")

    b_st, b_en = line_list_block_coords(
        lines, s_key=s_key, e_key=e_key, inclusive=inclusive, stoe=stoe
    )

    # Set all bounds for slicing; b_st and b_en are < 0 if not set/found
    # If expected block coords not found, fail
    if s_key and e_key:
        if (b_st < 0) or (b_en < 0):
            return []
    elif s_key:
        if b_st < 0:
            return []
        b_en = len(lines)
    elif e_key:
        if b_en < 0:
            return []
        b_st = 0

    if inverse:
        new_lines = lines[:b_st] + lines[b_en:]
    else:
        new_lines = lines[b_st:b_en]

    return new_lines


def simple_sub_list_lines(line_list, s_pat, r_pat):
    """Apply simple substitution of patterns in list of lines

    Return list of lines
    """
    new_lines = []
    for line in line_list:
        new_lines.append(line.replace(s_pat, r_pat))
    return new_lines


def update_nav_list_lines(line_list, page1, page2, page3):
    """Modify nav bar list in lines from html

    page1 = text for page 1 link
    page2 = text for page 2 link
    page3 = text for page 3 link

    Return updated line list
    """
    # Search key words
    start_key = "block_start_list_navbar"
    end_key = "block_end_list_navbar"
    new_lines = []

    in_block = False
    for line in line_list:
        # Block start marker
        if start_key in line:
            in_block = True
            new_lines.append(line)
            continue

        # Block end
        if end_key in line:
            in_block = False

            # Now add what's real
            if page1:
                new_lines.append(one_nav_link_line(1, page1))
            if page2:
                new_lines.append(one_nav_link_line(2, page2))
            if page3:
                new_lines.append(one_nav_link_line(3, page3))

            # Finish with end marker
            new_lines.append(line)
            continue

        # Ignore actual lines in block
        if in_block:
            continue

        new_lines.append(line)

    return new_lines


def one_nav_link_line(page_num, new_txt):
    """Get single nav link line

    Lines need to look like this:

        <li><a href="#" onclick="return show_page('Page1');">Summary</a></li>

    Return string
    """
    new_line = '<li><a href="#" onclick="return show_page('
    new_line += f"'Page{page_num}');"
    new_line += '">'
    new_line += f"{new_txt}</a></li>"

    return new_line
