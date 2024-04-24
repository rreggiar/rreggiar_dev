#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import os
import re


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


# ----------------------------------------------------------------------------
def get_mem_script_path():
    """Get path to memory script

    Return str
    """
    # Get top level path
    try:
        path = os.path.dirname(__file__)
    except Exception:
        path = "."

    # Specific location
    call_path = f"{path}/scripts/mem_loop.py"
    return call_path


def init_mem_profile(spipe, num=0, pid=None, seconds=2):
    """Initialize memory profile tracing

    num = (max) number of sampling points for memory script

    Return status
    """
    if not spipe.get_par_val("log_memory_use", as_bool=True):
        spipe.report_run_story("Not running memory profile tracing")
        return True

    spipe.report_run_story("Starting memory profile tracing")
    fname = spipe.filepath("PF_MEM_PROFILE", None)
    # Possibly bumping file versions?
    b = utils.bump_logs(fname)
    story = f"(bumped {b} prior version(s))" if b else ""
    spipe.report_run_story(f"Memory log {fname} {story}")

    # Header
    ver = spipe.get_version()
    utils.file_header(ofile=fname, story="Run time memory profile", ver=ver)
    with open(fname, "a") as OUTFILE:
        story = f"{spipe.get_sys_mem()} RAM and {spipe.get_sys_cpu()} CPUs"
        print(story, file=OUTFILE)

    # If not given PID, get it
    if not pid:
        pid = os.getpid()
    # Limit number of samples so not forever; 1 day = 86400 seconds
    if not num:
        num = 86400 * 7
    # Command line to run
    call_path = get_mem_script_path()
    cmd_str = f"{call_path} \
        --outfile {fname} \
        --number {num} \
        --seconds {seconds} \
        --pid {pid} \
        --quiet \
        --append \
        --prefix 'Tstep' \
        "
    # Clean to single space (not printing as above)
    cmd_str = re.sub(" +", " ", cmd_str)

    # Launch detached; Only need pass/fail and story (no output)
    spipe.report_run_story2(f"Call: {cmd_str}")
    ok, _, story = utils.call_exec(cmd_str, detach=True, verb=False)
    if not ok:
        spipe.set_problem("Failed to start memory profile script")
        spipe.set_problem(story)

    return ok


def report_to_mem_profile(spipe, story):
    """Initialize memory profile tracing

    Return status
    """
    if not spipe.get_par_val("log_memory_use", as_bool=True):
        return True

    ok = True
    fname = spipe.filepath("PF_MEM_PROFILE", None)
    try:
        OUTFILE = open(fname, "a")
        print(f"Tnote\t{story}", file=OUTFILE)
    except Exception:
        ok = False

    return ok
