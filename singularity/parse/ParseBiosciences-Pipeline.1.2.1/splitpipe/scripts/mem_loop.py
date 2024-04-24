#!/usr/bin/env python
# 2021-10-13 RTK
# 2023-09-06 RTK; v0.2 add --prefix and --outfile
# 2023-12-06 RTK; v0.3 add --quiet and --append
#
# Call to report memory use in a loop
#

import datetime
import argparse
import psutil
import sys
import os
import time


__version__ = "V0.3; RTK 2023-12-06"

DEF_STEP = 2  # Default sampling step; Seconds

DEF_UPDATE = 15  # Default update step number if output to file


# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Simple memory monitor tool")
    parser.add_argument(
        "-n", "--number", type=int, default=0, help="Number of sample steps (then stop)"
    )
    parser.add_argument(
        "-s",
        "--seconds",
        type=float,
        default=DEF_STEP,
        help="Sleep length (seconds; fractional OK)",
    )
    parser.add_argument(
        "--no_count", action="store_true", help="Don't print count for each sample step"
    )
    parser.add_argument(
        "--no_sys", action="store_true", help="Don't report system memory"
    )
    parser.add_argument(
        "--no_proc", action="store_true", help="Don't report process memory"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="No reported feedback; Only time data lines",
    )
    parser.add_argument("-p", "--pid", type=int, help="Process id")
    parser.add_argument(
        "-q", "--pname", help="Process name; Given must match process name"
    )
    parser.add_argument(
        "-Q",
        "--partname",
        help="Partial process name; Give a substring of process name(s)",
    )
    parser.add_argument(
        "--no_kids",
        action="store_true",
        help="Do not include any children processes; Default = include",
    )
    parser.add_argument("--prefix", default="", help="Prefix for output lines")
    parser.add_argument(
        "--update",
        type=int,
        default=DEF_UPDATE,
        help="Update frequency if output to file; Number of steps; Zero for none",
    )
    parser.add_argument("-o", "--outfile", help="Output filename")
    parser.add_argument("-a", "--append", action="store_true", help="Append outfile")
    parser.add_argument(
        "-V", "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.parse_args()

    number = args.number
    seconds = args.seconds
    children = False if args.no_kids else True

    pref = ""
    pid = -1
    pname = ""
    if args.pid or args.pname or args.partname:
        pid, pname = find_process(args)
        if pid < 0:
            print("Process to monitor not found")
            return
        pref = f"{pname} ({pid})"

    out_proc = False if args.no_proc else True
    out_sys = False if args.no_sys else True

    # Feedback
    if number:
        story = f"# Total {number} samples, {seconds} second steps"
    else:
        story = f"# Infinite loop every {seconds} seconds"
    if not args.quiet:
        print(story)

    # If output file, init write here
    # Note that if file is *held* open to write, other processes cannot also
    #   write to the file. Each step is thus *appended* to any output file so
    #   other processes may write to output file (i.e. inject marker text)
    if args.outfile:
        mode = "a" if args.append else "w"
        OUTFILE = open(args.outfile, mode)
        print(story, file=OUTFILE)
        OUTFILE.close()
        if not args.quiet:
            print(f"# Output to {args.outfile} ({mode})")
    else:
        if not args.quiet:
            print("# Output to sdtout")

    n = 0
    while True:
        n += 1
        if number and (n > number):
            break

        # Update? Only if not writing to file and non-zero step
        if args.outfile and args.update and ((n % args.update) == 0):
            if not args.quiet:
                print(f"# Step {n}")

        # If specific process, make sure it's still alive
        if pid > -1:
            if not find_proc_by_pid(pid, verb=False):
                if not args.quiet:
                    print(f"Process {pid} doesn't exist; done")
                break

        mstory = report_mem_str(
            pref=pref, out_proc=out_proc, out_sys=out_sys, pid=pid, children=children
        )
        tstory = get_time_str()

        o_parts = []
        if args.prefix:
            o_parts.append(args.prefix)
        if not args.no_count:
            o_parts.append(str(n))
        o_parts.append(tstory)
        o_parts.append(mstory)
        story = "\t".join(o_parts)

        # Appending file?
        if args.outfile:
            OUTFILE = open(args.outfile, "a")
            print(story, file=OUTFILE)
            OUTFILE.close()
        else:
            print(story)
            sys.stdout.flush()

        time.sleep(seconds)

    # Clean up
    if args.outfile and (not args.quiet):
        print("New file:", args.outfile)


def find_process(comargs, single=True):
    """Try to find process by id or name

    If found, return tuple (pid, name); If not found (-1, '')
    """
    pid = -1
    pname = ""
    proc = None
    if comargs.pid:
        proc = find_proc_by_pid(comargs.pid, verb=True)
    elif comargs.pname or comargs.partname:
        if comargs.partname:
            m_name = comargs.partname
        else:
            m_name = comargs.pname
        # Go over all processes comparing names
        matches = []
        for t_proc in psutil.process_iter():
            if comargs.partname:
                if m_name in t_proc.name():
                    matches.append(t_proc)
            else:
                if m_name == t_proc.name():
                    matches.append(t_proc)
        if matches:
            if (len(matches) > 1) and single:
                print(f"Multiple ({len(matches)}) processes match name '{m_name}':")
                for t_proc in matches:
                    print(f"    {t_proc.pid}\t{t_proc.name()}")
            else:
                proc = matches[0]
        else:
            print(f"Could not find process matching name '{m_name}'")
    if proc:
        pid = proc.pid
        pname = proc.name()
    return (pid, pname)


def find_proc_by_pid(pid, verb=True):
    """Attempt to get process by id

    Return process object
    """
    proc = None
    try:
        proc = psutil.Process(pid)
    except Exception as e:
        if verb:
            print(f"Could not find process with pid {pid}")
            print(e)
    return proc


def get_proc_mem_info(pid=-1, children=True, verb=False):
    """Get memory info, possibly for given process and children

    Return tuple (rss, vms)
    """
    # Get process
    if pid >= 0:
        try:
            proc = psutil.Process(pid)
        except Exception as e:
            if verb:
                print(f"Failed to find process pid={pid}; {e}")
            return None
    else:
        proc = psutil.Process(os.getpid())

    # Memory parts; https://psutil.readthedocs.io/en/latest/index.html?highlight=memory_info#psutil.Process.memory_info
    # rss = “Resident Set Size”, this is the non-swapped physical memory a process has used.
    rss = proc.memory_info().rss
    # vms = “Virtual Memory Size”, this is the total amount of virtual memory used by the process.
    vms = proc.memory_info().vms

    # Include any children processes?
    if children:
        for kid in proc.children(recursive=True):
            rss += kid.memory_info().rss
            vms += kid.memory_info().vms

    return (rss, vms)


def report_mem_str(
    pref="Memory status:", out_proc=True, out_sys=True, pid=-1, children=True
):
    """Get 'nice' string with (process) memory use and (sys) resources

    return str
    """
    mem_info = get_proc_mem_info(pid=pid, children=children)
    if not mem_info:
        return f"Failed to get memory info for pid={pid}"
    rss, vms = mem_info

    proc_tot = rss + vms
    rss_str = report_bytes_str(rss, full=False)
    vms_str = report_bytes_str(vms, full=False)
    proc_str = report_bytes_str(proc_tot, full=False)
    per_str = f"{psutil.virtual_memory().percent}%"
    avail_str = report_bytes_str(psutil.virtual_memory().available, full=False)

    rep_str = pref
    if out_proc and out_sys:
        rep_str = f"{pref} {proc_str} proc, {rss_str} rss, {vms_str} vms, {per_str} usd, {avail_str} avail"
    elif out_proc:
        rep_str = f"{pref} Process {proc_str}, {rss_str} rss, {vms_str} vms"
    elif out_sys:
        rep_str = f"{pref} System {per_str} used, {avail_str} available"
    return rep_str


def report_bytes_str(n, full=True):
    """Get 'nice' string to report bytes

    return str
    """
    # Based on http://code.activestate.com/recipes/578019
    n = int(n)
    symbols = ("K", "M", "G", "T", "P", "E", "Z", "Y")
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            val = float(n) / prefix[s]
            rep_str = f"{val:.1f}{s}"
            if full:
                rep_str += f" ({n})"
            # return '%.1f%s' % (value, s)
            return rep_str

    return f"{n}B"


def get_time_str(stime=None, date=True):
    """Get (pretty?) string for time

    stime = start time to get delta from; If none, use current time
    """
    if stime:
        t = datetime.datetime.now() - stime
    else:
        t = datetime.datetime.now()
    # Truncate precision via cutting string
    t_str = str(t)[:-4]
    if not date:
        t_str = t_str.split()[1]
    return t_str


# Actually call program
if __name__ == "__main__":
    main()
