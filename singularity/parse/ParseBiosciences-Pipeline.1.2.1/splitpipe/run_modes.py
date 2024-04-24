#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

from collections import OrderedDict

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


# ---------------------- Mode (run step) stuff ----------------------
# List of run modes and calls; e.g. mkref >--> mkref.run_mkref(spipe)
# Note: keys here define the 'cannonical' spelling for modes
RUN_MODE_DEFS = """
    # Add all as needed for parsing, etc; No actual run step
    all         None                    All
    
    # mkref is special, single mode and use case, 1 step
    mkref	    mkref.run_mkref         MkRef

    # Combine mode is special; Subset of downstream steps follow
    comb	    combine.run_combine     Combine

    # Normal 'all' mode starts with steps here
    inqc	    inqc.run_inqc           InQC
    pre	        align.run_preprocess    PreProcess
    align	    align.run_align         Align
    post	    align.run_postprocess   PostProcess
    mol	        molinfo.run_molinfo     Molinfo
    # Choice here; Either TCR or DGE, not both
    tcr	        tcr.run_tcr             TCR
    dge	        dge.run_dge             DGE

    ana	        analysis.run_analysis   Analysis
    rep	        report.run_reports      Report
    clean	    clean.run_clean         Clean
    """


def _make_mode_run_dict():
    """Create dict with run [mode] = [call, name]"""
    run_dict = OrderedDict()
    for line in RUN_MODE_DEFS.split("\n"):
        parts = line.split("#")[0].split()
        if len(parts) > 2:
            run_dict[parts[0]] = [parts[1], parts[2]]
    return run_dict


# Create once as global
MODE_RUN_DICT = _make_mode_run_dict()
MODE_ORDER_LIST = list(MODE_RUN_DICT.keys())


def parse_mode(mode):
    """Attempt to parse given mode into cannonical string"""
    mode_str = ""
    # Given string can be any case and longer
    for m in MODE_ORDER_LIST:
        if mode.lower().startswith(m):
            mode_str = m
            break
    return mode_str


def run_mode_index(mode, toxic=True):
    """Get run mode index (i.e. order or running)

    Return -1 if fail, else 0-base index
    """
    r_ind = -1
    p_mode = parse_mode(mode)
    if p_mode:
        r_ind = MODE_ORDER_LIST.index(p_mode)
    elif toxic:
        raise ValueError(f"Unrecognized mode: {mode}")
    return r_ind


def compare_mode_order(mode1, mode2):
    """Compare run steps with respect to run order; before vs after

    Return -1, 0, 1
    """
    step1 = run_mode_index(mode1)
    step2 = run_mode_index(mode2)
    if step1 < step2:
        return -1
    elif step1 > step2:
        return 1
    return 0


def raw_run_step_list(spipe):
    """Get list of run steps (aka "modes") to be called

    Uses given mode, use_case, etc in spipe to decide full-run steps

    Return list
    """
    run_steps = []

    # Mkref is special
    if spipe.is_mkref():
        run_steps = ["mkref"]

    # Combine mode
    elif spipe.is_combine():
        # Start at combine
        run_steps = ["combine"]
        # Next step depends on DGE / TCR
        if spipe.is_tcr():
            next_step = "report"
        else:
            next_step = "analysis"
        # Append rest of run steps
        next_ind = run_mode_index(next_step)
        run_steps += MODE_ORDER_LIST[next_ind:]

    # Normal (e.g. 'all') mode
    else:
        # input QC and preprocessing
        run_steps = ["inQC", "preproc"]
        if spipe.is_tcr():
            run_steps += ["TCR"]
        else:
            run_steps += ["align", "post", "molinfo", "DGE", "analysis"]
        # Report to end
        next_ind = run_mode_index("report")
        run_steps += MODE_ORDER_LIST[next_ind:]

    # Make sure every step is the cannonical name
    run_steps = [parse_mode(mode) for mode in run_steps]
    return run_steps


def run_step_to_call_name(run_step, get_call=False, get_name=False):
    """Get function call fro given run step or list of steps

    Return str or list[str]
    """
    # If neither call or name is specified, return both
    if (not get_call) and (not get_name):
        get_call = get_name = True

    if isinstance(run_step, list):
        call_list = []
        for mode in run_step:
            call = run_step_to_call_name(mode, get_call=get_call, get_name=get_name)
            if not call:
                return []
            call_list.append(call)
        return call_list

    else:
        answer = ""
        mode = parse_mode(run_step)
        if mode in MODE_RUN_DICT:
            # Dict has tuples with (call, name)
            if get_call and get_name:
                answer = MODE_RUN_DICT[mode]
            elif get_call:
                answer = MODE_RUN_DICT[mode][0]
            elif get_name:
                answer = MODE_RUN_DICT[mode][1]
        return answer


# ---------------------- Use case stuff -------------------------------------


def decide_use_case(spipe):
    """Decide use case from parameter settings

    Return use_case string
    """
    # Decide based on other parameters
    is_targeted = spipe.is_targeted()
    have_parent = spipe.have_parent()
    is_crispr = spipe.is_focal_crispr(crispr=True)
    is_focal = spipe.is_focal_crispr(focal=True)
    is_tcr = spipe.get_par_val("tcr_analysis", as_bool=True)

    # Different combos of input settings determine use case
    if is_tcr:
        if have_parent:
            use_case = "tcr_parent"
        else:
            use_case = "tcr_only"
    elif is_crispr:
        use_case = "crispr"
    elif is_focal:
        use_case = "focal_bc"
    elif is_targeted:
        if have_parent:
            use_case = "target_parent"
        else:
            use_case = "target_only"
    else:
        use_case = "normal"

    return use_case
