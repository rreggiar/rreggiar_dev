#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Misc util funcitons

import os
import getpass
import sys
import io
import re
import psutil
import shutil
import subprocess
import datetime
import time
import json
import csv
import gzip
import numpy as np
import pandas as pd
import scipy
import scipy.io as sio
import scanpy
import anndata
import random
from collections import OrderedDict

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
else:
    from splitpipe import bcutils



# ---------------------------------------------------------------------------

# Reverse complement
DNA_RC_DICT = dict(zip(list("NACGT"), list("NTGCA")))


def reverse_complement(seq):
    """Get sequence reverse compliment"""
    return "".join([DNA_RC_DICT[s] for s in seq][::-1])


# Regex
# Matches float/int with optional sign and chars before / after
#   This part (?<![+-]) is not a capture; Instead prevents \D matching +-
#   From https://stackoverflow.com/questions/58841472/regex-for-digits-and-plus-minus-sign
NUM_REGEX = re.compile(r"\D*(?<![+-])([+-]?([0-9]*)(\.([0-9]+))?)\D*$")


def str_to_num(
    string, as_float=False, commas=True, regex=False, fail_val=None, nan_ok=True
):
    """Attempt to get number from string

    string = string to parse for number
    as_float = flag to return float; Else return int
    commas = flag to allow commas in number
    regex = flag to use (above) regular expression; Find number in string
    fail_val = value to return on failure
    nan_ok = flag to allow NaN values

    Return tuple (success, value)
    """
    ok = False
    # Optional clean up / extraction
    if isinstance(string, str):
        # Strip any commas?
        if commas:
            string = string.replace(",", "")
        # regex; Allows surrounding chars (e.g. 'abc-122.4x' >--> -122.4)
        if regex:
            mat = NUM_REGEX.match(string)
            if mat:
                string = mat[1]
    # Cast
    try:
        # Cast to float first regardless; then maybe int
        n = float(string)
        if not as_float:
            n = int(string)
        ok = True
    except Exception:
        if fail_val is None:
            n = string
        else:
            n = fail_val

    # Check NaN
    if ok and np.isnan(n) and (not nan_ok):
        ok = False
        if fail_val is not None:
            n = fail_val

    return ok, n


def whole_float_to_int(number):
    """Convert whole number float to int, else leave alone"""
    if number.is_integer():
        return int(number)
    else:
        return number


def str_to_bool(string, one_letter=False):
    """Parse string as boolean

    Return tuple (success, value)
    """
    ok = False
    answer = None
    if string:
        yes_list = ["TRUE", "YES", "GOOD"]
        no_list = ["FALSE", "NO", "BAD"]
        if one_letter:
            yes_list += ["T", "Y"]
            no_list += ["F", "N"]
        # In case non-string
        try:
            if str(string).upper() in yes_list:
                answer = True
                ok = True
            elif str(string).upper() in no_list:
                answer = False
                ok = True
        except Exception:
            pass
    return ok, answer


def list_min_max(vlist, as_float=False, commas=True):
    """Get min and max value from list (iterable; set works)

    Return number of items checked, min, max
    """
    n = 0
    min_val = max_val = None
    for v in vlist:
        ok, val = str_to_num(v, as_float=as_float, commas=commas)
        if ok:
            n += 1
            if min_val is None:
                min_val = max_val = val
            else:
                min_val = min(val, min_val)
                max_val = max(val, max_val)
    return n, min_val, max_val


def list_find_match(v_list, m_list, nocase=True, startswith=False, strip=True):
    """Look for string or list of strings in list

    v_list = list of strings to check
    m_list = single or list of strings to check for match
    nocase = flag to ignore case
    startswith = flag to check start

    If string is found, return that; No match empty string

    Return str, index
    """
    # Ensure match is a list of strings
    if isinstance(m_list, str):
        m_list = [m_list]
    m2_list = []
    for s in m_list:
        s = str(s)
        if strip:
            s = s.strip()
        if nocase:
            s = s.upper()
        m2_list.append(s)

    # Screen value list
    for i, v in enumerate(v_list):
        v_str = str(v)
        if strip:
            v_str = v_str.strip()
        if nocase:
            v_str = v_str.upper()
        hit = False
        for s in m2_list:
            hit = False
            if startswith:
                if v_str.startswith(s):
                    hit = True
            else:
                if v_str == s:
                    hit = True
            if hit:
                return v, i

    return "", -1


def parse_range_list(rlis, warn=True):
    """Parse range list; maybe [], [end], [start,end]

    first and last default zero; If given, set as int

    Return status,first,last
    """
    ok = True
    first = last = 0
    if rlis:
        if len(rlis) > 1:
            if (len(rlis) > 2) and warn:
                print(f"Extra range args ignored; {rlis}")
            ok, first = str_to_num(rlis[0])
            if ok:
                ok, last = str_to_num(rlis[1])
        else:
            ok, last = str_to_num(rlis[0])
    return ok, first, last


def pad_line(line, plen):
    """Pad line to given width"""
    # Format to length
    fmt_str = f"{{line:{plen}s}}"
    new_line = fmt_str.format(line=line)
    return new_line[:plen]


def print_now(story, newline=True, file=None):
    """Print and immediately flush buffer"""
    OUTFILE, _ = get_out_fh(ofname=file)
    if newline:
        print(story, file=OUTFILE)
    else:
        print(story, end="", file=OUTFILE)
    OUTFILE.flush()


def print_story(story, banner=True, header="NOTE:", newlines=True, file=None):
    """Print story to get attention"""
    ostring = ""
    bstring = "=" * 80 + "\n"
    if newlines:
        ostring += "\n"
    if banner:
        ostring += bstring
    if header:
        ostring += header + "\n"
    ostring += story + "\n"
    if banner:
        ostring += bstring
    if newlines:
        ostring += "\n"
    print_now(ostring, newline=True, file=file)


def report_range_str(r_ends, arg2=None, lower=False):
    """Get string reporting range restrictions"""
    # Collect possible first and last limits
    first = last = ""
    # Given int(s)
    if isinstance(r_ends, int):
        if isinstance(arg2, int):
            first = r_ends
            last = arg2
        else:
            last = r_ends
    # Else given list or tuple
    elif r_ends:
        if len(r_ends) > 1:
            first, last = r_ends[:2]
        else:
            last = r_ends[0]
    # Now compose story
    if first and last:
        story = f"From {first} to {last}"
    elif last:
        story = f"Range to max {last}"
    else:
        story = "No range restrictions"
    # Lowercase?
    if lower:
        story = story.lower()
    return story


def report_num_str(num, big_num=1, round_to=3, nan="NAN?", zero=True):
    """Format number for reporting

    big_num = threshold for numbers to be "big" and rounded to int
    round_to = rounding precision for non-big numbers
    nan     = string for non-number arguments; If None, crash on non-number
    zero    = flag to report zero as just that, without any decimals

    Return string
    """
    try:
        n = float(num)
        if n >= big_num:
            s = "{:,}".format(round(n))
        else:
            if (n == 0) and zero:
                s = "0"
            else:
                fmt = "{:." + str(round_to) + "f}"
                s = fmt.format(n)
    except Exception:
        if nan is None:
            story = "report_number_string given non-number to format: {0}".format(num)
            raise ValueError(story)
        s = nan
    return s


def report_percent_str(num, den=1, round_to=2, zero=np.NAN, perchar=True):
    """Get pretty precent string

    num = numerator; If this is only argument, can convert fract to percent
    den = denominator
    round_to = decimal places
    perchar = flag to include percent char %

    Return string
    """
    try:
        pc = "%" if perchar else ""
        num = float(num) * 10**round_to
        den = float(den) / 100 * 10**round_to
        rep_str = f"{round(num/den, round_to)}{pc}"
    except Exception:
        rep_str = f"{zero}"
    return rep_str


def get_time(stime=None, date=True, as_str=True):
    """Get (pretty?) string for time

    stime = start time to get delta from; If none, use current time
    """
    if stime:
        t = datetime.datetime.now(datetime.timezone.utc) - stime
    else:
        t = datetime.datetime.now(datetime.timezone.utc)
    # String?
    if as_str:
        # truncate precision via cutting string
        t_str = str(t)[:-4]
        if not date:
            t_str = t_str.split()[1]
        t = t_str
    return t


def report_time(stime=None, pref="# Time", msg="", sep="\t", ofname=None):
    """Report current time

    pref    = Prefix to print
    msg     = Message
    """
    parts = [pref]
    if msg:
        parts.append(msg)
    parts.append(get_time(stime=stime, as_str=True))
    print_now(sep.join(parts), file=ofname)


# One letter chars for months
MONTH_1CHAR = list("JFRPMUYASOND")
LETTERS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")


def ymdm_timestamp(month_letter=False, min_letter=True):
    """Get year, month, day, minute encoded timestamp

    month_letter = Flag to encode month with 1-letter code; Days single digit
    min_letter = Flag to encode minutes with letters

    Format as YYMmDdMm:
        Year = two digits YY
        Month = two digit or single letter
        Day = two digit or single digit if paired with month letter code
        Minutes = digits or letter code

    Return string
    """
    # Universal time
    t = datetime.datetime.now(datetime.timezone.utc)
    # Only last two digits of year
    year = str(t.year)[-2:]
    # Month as one char code; Day as simple int
    if month_letter:
        month = MONTH_1CHAR[t.month - 1]
        day = str(t.day)
    else:
        month = f"{t.month:02d}"
        day = f"{t.day:02d}"
    # Encode minutes into two chars (52-letter alphabet)
    minutes = t.hour * 60 + t.minute
    if min_letter:
        min1 = int(minutes / 52)
        min2 = minutes % 52
        minutes = LETTERS[min1] + LETTERS[min2]
    else:
        minutes = str(minutes)

    return year + month + day + minutes


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


def report_disk_str(path, pref="Disk status:"):
    """Get 'nice' string with disk space numbers

    return str
    """
    # Expand any env vars, etc to full explicit path; Strip leading slash
    fpath = fname_full_path(path)
    if not fpath:
        return f"Bad path (report_disk_str): {path}"

    if fpath == "/":
        top = "/"
    else:
        top = "/" + fpath[1:].split("/")[0] + "/"
    # Actual disk usage
    total, used, free = shutil.disk_usage(fpath)
    per_u = report_percent_str(used / total)
    per_f = report_percent_str(free / total)
    t_bytes = report_bytes_str(total, full=False)
    u_bytes = report_bytes_str(used, full=False)
    f_bytes = report_bytes_str(free, full=False)
    rep_str = f"{pref} '{top}' {t_bytes} total, {u_bytes} ({per_u}) used, {f_bytes} ({per_f}) free"
    return rep_str


def report_file_info_str(fname, name=True, basename=False, size=True, do_time=True):
    """Get info about file (if exists)

    return str
    """
    rep_str = ""
    if os.path.exists(fname):
        if name:
            if basename:
                rep_str = os.path.basename(fname)
            else:
                rep_str = f"{fname}"
        if size:
            bts = os.path.getsize(fname)
            bt_str = report_bytes_str(bts, full=False)
            prefix = ", " if rep_str else ""
            rep_str += f"{prefix}{bt_str}, ({bts})"
        if do_time:
            date = time.strftime("%m/%d/%Y", time.gmtime(os.path.getctime(fname)))
            prefix = ", " if rep_str else ""
            rep_str += f"{prefix}{date}"
    return rep_str


def report_dict_str(
    dic, sep="\t", pref="", unwrap="", uw_simp=False, usep="\t", k_fmt=None, v_fmt=None
):
    """Get dictionary contents into formatted string

    sep     = Separator for key and val
    pref    = Prefix (for each line)
    unwrap  = Flag / str to write iterable values on per line; If empty, one line
    uw_simp = Flag to unwrap list/dict of simple types
    usep    = Separator for unwrapping (not between key and val)
    k_fmt   = Format string for key; "{k}" (key; Should add separator if sep='')
    v_fmt   = Format string for val; "{v}"

    Return string
    """
    str_rows = []
    # If no format strings, make them
    if not k_fmt:
        k_fmt = "{k}" + sep
    if not v_fmt:
        v_fmt = "{v}"
    # Add prefix to key format string
    if pref:
        k_fmt = pref + k_fmt

    # Unwrap 'simple' types
    # UW_SIMP_TYPES = (int, float, str, bool)
    UW_SIMP_TYPES = (int, float, bool)

    # Process dict
    for k, v in dic.items():
        # Determine if unwrapping
        do_unwrap = False
        if unwrap:
            # Get list of values
            uw_list = None
            if isinstance(v, list):
                uw_list = v
            elif isinstance(v, dict):
                uw_list = v.values()
            # Check list of values by type
            if uw_list:
                n_non_simp = 0
                for vv in uw_list:
                    if not isinstance(vv, UW_SIMP_TYPES):
                        n_non_simp += 1
                # If any non-simple types, or flag, then unwrap
                if n_non_simp or uw_simp:
                    do_unwrap = True

        # Actually process values
        if do_unwrap and isinstance(v, list):
            for i, v_item in enumerate(v):
                i_of = f"{unwrap}_{i+1}_of_{len(v)}" + usep
                o_str = k_fmt.format(k=k) + i_of + v_fmt.format(v=v_item)
                str_rows.append(o_str)
        elif do_unwrap and isinstance(v, dict):
            for i, v_key in enumerate(v.keys()):
                i_of = f"{unwrap}_{i+1}_of_{len(v)}" + usep
                o_str = (
                    k_fmt.format(k=k) + i_of + f"'{v_key}': " + v_fmt.format(v=v[v_key])
                )
                str_rows.append(o_str)
        else:
            o_str = k_fmt.format(k=k) + v_fmt.format(v=v)
            str_rows.append(o_str)

    return "\n".join(str_rows)


def report_list_str(lis, sep="\t", pref="", unwrap="item", v_fmt=None):
    """Get list or dict contents into formatted string

    sep     = Separator for key and val
    pref    = Prefix (for each line)
    unwrap  = String for each unwrapped item
    v_fmt   = Format string for val; "{v}"

    Return string
    """
    str_rows = []
    # If no format strings, make them
    if not v_fmt:
        v_fmt = "{v}"

    if not isinstance(lis, list):
        lis = list(lis)

    for i, v_item in enumerate(lis):
        i_of = f"{unwrap}_{i+1}_of_{len(lis)}" + sep
        o_str = pref + i_of + v_fmt.format(v=v_item)
        str_rows.append(o_str)

    return "\n".join(str_rows)


def report_df_info_str(df, name=None, dims=True, mem=True):
    """Reports dataframe properties

    Returns string
    """
    row, col = df.shape
    bts = df.memory_usage().sum()
    memory = report_bytes_str(bts, full=False)
    if name:
        rep_str = f"DF '{name}': "
    else:
        rep_str = "Dataframe: "
    if dims and mem:
        rep_str += f"Dims {row} x {col}, Mem {memory} ({bts})"
    elif dims:
        rep_str += f"Dims {row} x {col}"
    elif mem:
        rep_str += f"Mem {memory} ({bts})"
    return rep_str


def report_df_str(df, pref="# ", index_name="", n_rows=0, row_num=True):
    """Get string from dataframe for reporting

    df      = dataframe (or series)
    pref    = prefix string (start of each row)
    n_rows  = number of rows; If zero, dump all
    row_num = cap on max rows to print (output string)

    return string
    """
    # In case series, turn to dataframe
    df = pd_series_to_df(df)
    if index_name:
        df.index.name = index_name
    # String version of dataframe, broken by line
    df_lines = df.to_string().split("\n")
    new_lines = []
    for i, line in enumerate(df_lines):
        if n_rows and (i > n_rows):
            break
        if row_num:
            line = f"Row_{i+1}_of_{len(df_lines)}  " + line
        if pref:
            line = pref + line
        new_lines.append(line)
    return "\n".join(new_lines)


def check_df_col_dup_counts(df, col, max_count=1, pref=""):
    """Check if value counts of given col exceed limit

    For violations, creates string with <value> <count> per line

    Return count, str
    """
    over_str = ""
    col_counts = df[col].value_counts()
    dup_df = pd.DataFrame(col_counts[col_counts > max_count])
    n_over = len(dup_df)
    # If any over the limit, format as table (df) output
    if n_over:
        # Max name length to justify output; Extra for more space
        max_nlen = max([len(col)] + list(dup_df.index.str.len())) + 2
        # Cook up rows of table then into str
        rows = [f"{pref}{col.ljust(max_nlen)} Count"]
        for name, count in zip(dup_df.index, dup_df[col]):
            rows.append(f"{pref}{name.ljust(max_nlen)} {count}")
        over_str = "\n".join(rows)
    return n_over, over_str


def check_df_col_str_content(df, col, min_len=0, max_len=0, only_char=None, pref=""):
    """Check if dataframe column str content violates constraints

    df, col = dataframe and column
    min_len = minimum length constraint
    max_len = maximum length constraint
    only_char = string with only allowed characters
    pref = Prefix for lines of output string

    For violations, creates string with <index> <violation> per line

    Return count, str
    """
    bad_str = ""
    n_bad = 0

    min_set = set()
    max_set = set()
    bad_set = set()

    # Lengths
    if min_len:
        filt_df = df[df[col].str.len() < min_len]
        if len(filt_df):
            min_set.update(filt_df.index)
    if max_len:
        filt_df = df[df[col].str.len() > max_len]
        if len(filt_df):
            max_set.update(filt_df.index)
    # Restricted chars
    if only_char:
        pattern = f"^[{only_char}]+$"
        filt_df = df[~df[col].str.contains(pattern, regex=True)]
        if len(filt_df):
            bad_set.update(filt_df.index)

    v_set = min_set | max_set | bad_set
    n_bad = len(v_set)
    # Cook up string of violations
    if n_bad:
        # Max name, val lengths to justify output; Extra for more space
        max_nlen = max([len("Index")] + [len(str(i)) for i in v_set]) + 2
        bad_df = df[df[col].isin(v_set)]
        # Value plus four because want space and two quotes
        max_vlen = max([len("Value")] + list(bad_df[col].str.len())) + 4
        # Header then rest of violation strings
        v_list = [f"{pref}{'Index'.ljust(max_nlen)} {col} violation(s)"]
        for id in v_set:
            val = df.loc[id, col]
            probs = []
            probs.append("Too short" if id in min_set else "")
            probs.append("Too long" if id in max_set else "")
            probs.append("Bad char" if id in bad_set else "")
            v_str = ", ".join([p for p in probs if p])
            v_list.append(
                f"{pref}{str(id).ljust(max_nlen)}'{val.ljust(max_vlen)}' {v_str}"
            )
        bad_str = "\n".join(v_list)

    return n_bad, bad_str


def pd_series_to_df(df, cols=None):
    """Convert pandas series to dataframe, setting column, index names

    Return dataframe
    """
    in_series = False
    if isinstance(df, pd.Series):
        in_series = True
        df = pd.DataFrame(df)
        # Column names; As list or comma,sep,string
        if cols:
            if isinstance(cols, str):
                cols = cols.split(",")
            # If starting from series, set index and single col
            if in_series:
                df.index.name = cols[0]
                df.columns = cols[1:]
            else:
                df.columns = cols
    return df


### ----------------------------------------- file stuff  -----------------------------------------


def write_json(odict, fname, indent=4):
    """Write dict to json file"""
    with open(fname, "w") as OUTFILE:
        json.dump(odict, OUTFILE, indent=indent)
    return fname


def read_json(fname):
    """Parse json file"""
    with open(fname) as INFILE:
        ndict = json.load(INFILE)
    return ndict


def check_read_json(fname, verb=False):
    """Parse json file if found

    Return tuple (status, dict)
    """
    ndict = {}
    ok = check_infile(fname, verb=verb)
    if ok:
        ndict = read_json(fname)
    return ok, ndict


def get_header_dict(ver="", vernum="", desc="", fkey=""):
    """Get dict with header info (for json files)"""
    ndict = {}
    if desc:
        ndict["description"] = desc
    if fkey:
        ndict["filekey"] = fkey
    if ver:
        ndict["version"] = ver
    if vernum:
        ndict["ver_number"] = vernum
    ndict["date"] = get_time(date=True, as_str=True).split()[0]
    return ndict


def ver_at_least_ver(ver1, ver2):
    """Compare versions and return True if ver1 is at least ver2

    versions like 'v0.9.6' or 'v1.0.3p'

    Return boolean
    """
    # Get numbers parts only
    v1_num = v2_num = ""
    if ver1:
        v1_num = re.search(r"\s*([\d.]+)", ver1).group(1)
    if ver2:
        v2_num = re.search(r"\s*([\d.]+)", ver2).group(1)
    # Compare number parts
    if v1_num and v2_num:
        answer = True
        v1_nums = [int(i) for i in v1_num.split(".")]
        v2_nums = [int(i) for i in v2_num.split(".")]
        for i in range(len(v1_nums)):
            # Query greater than target = True
            if v1_nums[i] > v2_nums[i]:
                break
            # Query less than target = False
            if v1_nums[i] < v2_nums[i]:
                answer = False
                break
    else:
        answer = False
    return answer


def write_df(
    df,
    fname,
    sep=",",
    header=True,
    index=False,
    cols=None,
    mode="w",
    gzip=False,
    round_to=4,
    quoting=False,
    na_rep="NaN",
    verb=True,
):
    """Write pandas dataframe (or series) content to named file
    Wrapper around pandas IO

    df      = dataframe (or series)
    fname   = output filename
    sep     = field separator; passed to pandas
    header  = flag for header lines; passed to pandas
    index   = flag for index col; passed to pandas
    cols    = column names to write; list of str or comma,sep,string
    mode    = file open mode; 'a' for append
    gzip    = flag to write with gzip compression
    round_to = limit on precision; < 0 = ignore
    quoting = flag to quote strings
    na_rep  = how to represent NaN values (default = '')
    verb    = flag to report (return) string with results

    Return string: if verb, story of what's written, else fname
    """
    # Cast to dataframe (in case series); This handles col names
    df = pd_series_to_df(df, cols=cols)

    compression = "gzip" if gzip else "infer"
    quoting = csv.QUOTE_MINIMAL if quoting else csv.QUOTE_NONE
    # Limit precision
    ffmt_str = None
    if round_to >= 0:
        ffmt_str = f"%.{round_to}f"
    df.to_csv(
        fname,
        sep=sep,
        header=header,
        index=index,
        mode=mode,
        compression=compression,
        quoting=quoting,
        float_format=ffmt_str,
        na_rep=na_rep,
    )

    if verb:
        rows, cols = df.shape
        if index:
            cols += 1
        what = "Appended" if mode == "a" else "Wrote"
        fname = f"{what} {rows} x {cols} table to {fname}"
    return fname


def read_csv(
    fname,
    sep=",",
    index_col=None,
    names=None,
    usecols=None,
    head_if=None,
    header=0,
    dtype=None,
    verb=True,
):
    """Wrapper around pandas read function

    Return dataframe or None
    """
    if names:
        if isinstance(names, str):
            names = names.split(",")
    if usecols:
        if isinstance(usecols, str):
            usecols = usecols.split(",")

    # Header conditional on given name(s)
    if head_if:
        if isinstance(head_if, str):
            head_if = head_if.split(",")
        header = None

    df = None
    try:
        df = pd.read_csv(
            fname,
            sep=sep,
            comment="#",
            header=header,
            index_col=index_col,
            names=names,
            usecols=usecols,
            dtype=dtype,
        )
        # Conditional header, check for given name(s) in first row
        if head_if:
            head_found = False
            for col in head_if:
                if col in df.iloc[0].values:
                    head_found = True
                    break
            if head_found:
                # Drop first row and set columns; drop original header
                df = df.rename(columns=df.iloc[0]).drop(df.index[0])
    except Exception as e:
        if verb:
            print(f"Dataframe load failed; {e}")

    return df


def get_df_cutoff_counts(df, col, cutoff=1):
    """Get dataframe row count with col val of at least cutoff

    cutoff may be a list

    Return list
    """
    if not isinstance(cutoff, list):
        cutoff = [cutoff]
    counts = []
    for cut in cutoff:
        n = len(df[df[col] >= cut])
        counts.append(n)
    return counts


def guess_df_col_by_values(df, vals, max_missing=-1, min_frac=-1):
    """Guess dataframe column by values that should be contained

    Return column
    """
    vset = set(vals)
    best_col = ""
    best_missing = len(vset)
    # Scan columns looking for the fewest non-found (missing) vals
    for c in df.columns:
        cset = set(df[c].values)
        n_missing = len(vset - cset)
        if n_missing < best_missing:
            best_col = c
            best_missing = n_missing
            if best_missing == 0:
                break
    # Possibly qualify match
    if max_missing >= 0:
        if best_missing > max_missing:
            best_col = ""
    if min_frac >= 0:
        if best_missing / len(vset) < min_frac:
            best_col = ""
    return best_col


def read_mtx(fname, csc=False, verb=True):
    """Read sparse matrix from file

    csc = Flag for Compressed Sparse Col matrix; Else csr (Comp Sp Row)

    Return scipy matrix
    """
    mtx = None
    try:
        mtx = sio.mmread(fname)
        if csc:
            mtx = scipy.sparse.csc_matrix(mtx)
        else:
            mtx = scipy.sparse.csr_matrix(mtx)
        if verb:
            print(f"# Loaded matrix from {fname}")
    except Exception as e:
        if verb:
            print(f"Matrix load failed; {e}")
    return mtx


def write_mtx(mtx, fname, verb=True):
    """Write matrix in scipy "Matrix Market" sparse format

    Return string: if verb, story of what's written, else fname
    """
    rows, cols = mtx.shape
    comment = f"Rows=cells ({rows}), Cols=genes ({cols})"
    sio.mmwrite(fname, mtx, comment=comment)
    if verb:
        fname = f"Wrote {rows} by {cols} table to {fname}"
    return fname


def mtx_zero_less_than(mtx, thresh, copy=False):
    """Zero out scipy sparse matrix values less than threshold

    Return matrix
    """
    if copy:
        mtx = mtx.copy()
    try:
        # Because sparse, some tricks to make this efficient
        # https://seanlaw.github.io/2019/02/27/set-values-in-sparse-matrix/
        nonzero_mask = np.array(mtx[mtx.nonzero()] < thresh)[0]
        rows = mtx.nonzero()[0][nonzero_mask]
        cols = mtx.nonzero()[1][nonzero_mask]
        mtx[rows, cols] = 0
        # In place
        mtx.eliminate_zeros()
    except Exception as e:
        print(f"mtx_zero_less_than exception; {e}")
    return mtx


def mtx_stat_str(mtx, min_list=None):
    """Collect some stats about scipy sparse matrix

    Return string
    """
    if not min_list:
        min_list = [1, 3, 10]
    assert isinstance(min_list, list), f"Expected list got {type(min_list)}"

    row, col = mtx.shape
    total = row * col
    story = f"Matrix {row} x {col} = {total} total"
    for thresh in min_list:
        num = mtx[mtx >= thresh].shape[1]
        per_s = report_percent_str(num, total)
        story += f"; {num} ({per_s}) >= {thresh}"
    return story


def make_subdir(sub_dir, verb=True, toxic=True):
    """Check / create subdir

    sub_dir = filename or list of filenames
    verb = Flag to report problem
    toxic = Flag to crash on problem

    Returns status of subdir; exists or not?
    """
    # If given list, recur for each item
    if isinstance(sub_dir, list):
        f_lis = []
        for f in sub_dir:
            name = make_subdir(f, verb=verb, toxic=toxic)
            if name:
                f_lis.append(name)
            else:
                f_lis = None
                break
        return f_lis

    # Should be path at this point
    assert isinstance(sub_dir, str), f"Expected str got {type(sub_dir)}"
    # Only try to create if doesn't exist
    if not os.path.exists(sub_dir):
        try:
            os.makedirs(sub_dir)
        except Exception:
            story = f"Problem creating subdir |{sub_dir}|"
            if toxic:
                raise IOError(story)
            else:
                if verb:
                    print(story)
    return os.path.isdir(sub_dir)


def fnames_for_path(
    path="",
    depth=-1,
    get_files=True,
    get_subdirs=False,
    regex=None,
    name=None,
    start=None,
    ext=None,
    contain=None,
    qual_path=False,
    verb=False,
):
    """Collect qualifying files under path

    path        File path to look in and under
    depth       Max depth to keep
    get_files   Flag to consider files (vs subdirs)
    get_subdirs Flag to consider subdirs (vs files)
    regex       Regex for names to match
    name        Name for exact name match
    start       Prefix for name match
    ext         Extension for name match
    contain     Substring name match
    qual_path   Flag to inlcude path in match qualification tests

    return list of filenames
    """
    # No path default to current
    if not path:
        path = "."
    # In case restricting depth
    p_count = path.count("/")
    # If extension given, make sure it starts with dot
    if ext:
        if ext[0] != ".":
            ext = "." + ext
    # Find all files and (dis)qualify what to collect
    f_list = []
    for wpath, subdirs, files in os.walk(path):
        if get_files and get_subdirs:
            targs = subdirs + files
        elif get_files:
            targs = files
        elif get_subdirs:
            targs = subdirs
        # Each target name
        for fname in targs:
            # Qualify full path or just name?
            fullpath = os.path.join(wpath, fname)
            qualname = fullpath if qual_path else fname

            keep = True
            # Regex on possibly full qualify name
            if regex:
                if not re.match(regex, qualname):
                    keep = False
            # Name and start only for file itself
            if keep and name:
                if fname != name:
                    keep = False
            if keep and start:
                if not fname.startswith(start):
                    keep = False
            # Extension at end
            if keep and ext:
                if ext != os.path.splitext(qualname)[1]:
                    keep = False
            # Contain can match full qualify name
            if keep and contain:
                if contain not in qualname:
                    keep = False
            # Depth of file from top
            if keep and depth >= 0:
                f_count = fullpath.count("/") - p_count
                if f_count > depth:
                    keep = False
            if keep:
                f_list.append(fullpath)
    if verb:
        print("# Found", len(f_list), "matching files in", path)
    return f_list


def fname_list_to_dict(
    flist, regex=None, regex_key=1, ext=None, contain=None, qual_path=False
):
    """Get dict from list of files

    Return dict
    """
    new_dict = {}
    for f in flist:
        # Get parts to qualify Qualify
        parts = f.split("/")
        fname = parts[-1]
        bname, fext = filename, file_ext = os.path.splitext(fname)

        # Default dict key = basename
        key = bname
        # Qualify full path or just name?
        qualname = f if qual_path else fname
        keep = True
        # Regex on possibly full qualify name
        if regex:
            mat = re.match(regex, qualname)
            if mat:
                key = mat[regex_key]
            else:
                keep = False
        # Extension
        if keep and ext:
            if ext != fext:
                keep = False
        # Contain can match full qualify name
        if keep and contain:
            if contain not in qualname:
                keep = False

        if keep:
            new_dict[key] = f
    return new_dict


def bad_infile_list(flist, path=None):
    """Collect only 'bad' (non-exist/non-open) filenames from list"""
    bads = []
    for fname in flist:
        if not check_infile(fname, path=path, verb=False):
            bads.append(fname)
    return bads


def check_infile(fname, path=None, verb=True, toxic=False, dir_slash=True):
    """Check named input file(s) can be opened for read

    fname = filename or list of filenames
    path = Path to prepend
    verb = Flag to report problem
    toxic = Flag to crash on problem
    dir_slash = flag to make sure fname ends '/' if is a dir

    Returns filename if readable; None if not
    """
    # Check if given empty str or None first
    if not fname:
        story = "Given empty filename to check"
        if toxic:
            raise IOError(story)
        if verb:
            print(story)
        return None

    # If given filename list, recur for each item
    if isinstance(fname, list):
        f_lis = []
        for f in fname:
            name = check_infile(
                f, path=path, verb=verb, toxic=toxic, dir_slash=dir_slash
            )
            if name:
                f_lis.append(name)
            else:
                f_lis = None
                break
        return f_lis

    # At this point, single file name string
    assert isinstance(fname, str), f"Expected str got {type(fname)}"
    r_fname = None

    # Path to prepend?
    if path:
        f_path = fname_full_path(path, dir_slash=True)
        if f_path:
            fname = f_path + fname
        else:
            story = f"Given path is bad: '{path}'"
            if toxic:
                raise IOError(story)
            if verb:
                print(story)
            return None

    # Expand any env vars; User for ~/ and vars for $home
    fname = os.path.expanduser(fname)
    fname = os.path.expandvars(fname)
    # Read access means input is good
    if os.access(fname, os.R_OK):
        r_fname = fname
        # If dir, expand to explicit path
        if os.path.isdir(r_fname):
            r_fname = fname_full_path(r_fname, dir_slash=dir_slash)
    else:
        story = f"Cannot read file: '{fname}'"
        if toxic:
            raise IOError(story)
        if verb:
            print(story)
    return r_fname


def fname_full_path(fname, dir_slash=True, strip=True, exists=True):
    """Expand filename to full path

    fname = filename or list of filenames
    dir_slash = flag to make sure fname ends '/' if is a dir
    strip = flag to strip surrounding whitespace
    exists = flag to only return existing path

    return filename or None
    """
    # Check if given empty str or None first
    if not fname:
        return None

    # If given filename list, recur for each item
    if isinstance(fname, list):
        f_lis = []
        for f in fname:
            name = fname_full_path(f, dir_slash=dir_slash, exists=exists)
            if name:
                f_lis.append(name)
            else:
                f_lis = None
                break
        return f_lis

    # Single file name string
    assert isinstance(fname, str), f"Expected str got {type(fname)}"
    e_fname = None
    if strip:
        fname = fname.strip()
    # Expand any env vars; User for ~/ and vars for $home
    fname = os.path.expanduser(fname)
    fname = os.path.expandvars(fname)
    # This will strip out repeated, trailing slashes;
    #   '/path//to//x' >--> '/path/to/x'
    #   '/path//to//x/' >--> '/path/to/x'
    fname = os.path.normpath(fname)

    # Must be real file?
    if exists:
        if os.path.exists(fname):
            e_fname = fname
            # If dir, maybe need trailing slash
            if os.path.isdir(e_fname) and dir_slash:
                # Special case (e.g. normpath leaves '/' untouched)
                if e_fname[-1] != "/":
                    e_fname += "/"
    else:
        e_fname = fname
        if dir_slash:
            e_fname = e_fname + "/"
    return e_fname


def same_path(fname, sname):
    """Compare paths to see if same"""
    fpath = fname_full_path(fname, exists=False)
    spath = fname_full_path(sname, exists=False)
    return bool(fpath == spath)


def fname_dir_base(fname, dir_slash=True):
    """Return tuple with dir and basename"""
    dpart = os.path.dirname(fname)
    if dir_slash:
        dpart += "/"
    bname = os.path.basename(fname)
    return (dpart, bname)


def clean_vname_chars(
    name, to_dash=False, to_under=True, dash_ok=False, int_first=False, verb=True
):
    """Clean up given name to variable-like form

    to_dash = flag to substitute 'bad' chars to dash '-'
    to_under = flag to substitute 'bad' chars to underscore '_'
    dash_ok = flag to allow dash in name; else 'bad'
    int_first = flag to allow name to start with int; else 'bad'

    Return string or list
    """
    # If list, process each
    if isinstance(name, list):
        new_names = []
        for n in name:
            n_name = clean_vname_chars(
                n, to_dash=to_dash, to_under=to_under, int_first=int_first, verb=verb
            )
            if not n_name:
                return
            new_names.append(n_name)
        return new_names
    # In case naked number, check and cast to string
    ok, num = str_to_num(name)
    if ok:
        name = str(name)
    if not isinstance(name, str):
        story = f"name not string or list; type: {0}".format(type(name))
        raise Exception(story)
    # Replace non-target chars with target
    if to_dash:
        n_name = re.sub(r"[^a-zA-Z0-9-]", "-", name)
    elif to_under:
        if dash_ok:
            n_name = re.sub(r"[^a-zA-Z0-9-]", "_", name)
        else:
            n_name = re.sub(r"[^a-zA-Z0-9]", "_", name)
    else:
        n_name = re.sub(r"[^a-zA-Z0-9]", "", name)
    # Check first char; Can only be letter or underscore; Int if flag
    if int_first:
        regex = r"[^a-zA-Z_0-9]"
    else:
        regex = r"[^a-zA-Z_]"
    if re.match(regex, n_name[0]):
        if verb:
            print(
                f"Failed to clean '{name}' (clean_vname_chars, under={to_under}, dash={to_dash}, int_first={int_first})"
            )
        return False
    return n_name


def clean_fname_chars(fname, new="_", no_dupe_new=False):
    """Clean given name to only given regex characters

    fname = filename only; Should not be full path ('/' get cleaned)
    new = replacement char
    no_dup_new = flag to disallow repeated new replacements

    Return str
    """
    # Minimally return replacement char
    fnew = new
    fname = str(fname)
    if len(fname) > 0:
        # Replace any non-alpha numeric plus dash with new char
        fnew = re.sub("[^a-zA-Z0-9-]", new, fname)
        # Remove duplicate new chars?
        if no_dupe_new:
            # From: https://stackoverflow.com/questions/49695477/removing-specific-duplicated-characters-from-a-string-in-python
            pattern = "(?P<char>[" + re.escape(new) + "])(?P=char)+"
            fnew = re.sub(pattern, r"\1", fnew)
        # Remove new char if this is first
        if len(fnew) > 0:
            if fnew[0] == new:
                fnew = fnew[1:]
    return fnew


def check_file_gz(fname, path=None, verb=True, toxic=True):
    """Check if a gzip file exists; May be gz or not

    path = Path to prepend
    verb = Flag to report problem
    toxic = Flag to crash on problem

    Returns filename if exists; None if not
    """
    # Try filename as given; No verb, toxic until later; No need to open
    r_fname = check_infile(fname, path=path, verb=False, toxic=False)
    if not r_fname:
        # If given .gz file, strip this; Else add .gz and try again
        if fname.endswith(".gz"):
            n_fname = fname[:-3]
        else:
            n_fname = fname + ".gz"
        r_fname = check_infile(n_fname, path=path, verb=False, toxic=False)

    if not r_fname:
        story = f"Cannot read file: '{fname}' path '{path}'"
        if toxic:
            raise IOError(story)
        if verb:
            print(story)

    return r_fname


def check_exec(excom, timeout=None):
    """Check if named executable command exists; i.e. 'which excom'

    Return tuple (status, path, story)
    """
    com_args = f"which {excom}"
    ok, callout, story = call_exec(com_args, timeout=timeout, verb=False)
    callout = callout.strip()
    return ok, callout, story


def get_exec_path_and_ver(exe, ver_call, ver_line=1, timeout=None):
    """Attempt to get full path and version info for executable

    exe = excutable command to call
    ver_call = arguments to get version info
    ver_line = line of output to keep for version info (1-based)

    Issues to report are returned in story variable

    return tuple (status, path, version, story)
    """
    ver = ""
    ok, path, story = check_exec(exe, timeout=timeout)
    if ok:
        com_args = f"{exe} {ver_call}"
        ok, output, story = call_exec(com_args, timeout=timeout, verb=False)
        if ok:
            lines = output.split("\n")
            try:
                ver = lines[ver_line - 1]
            except Exception as e:
                ok = False
                story = f"Parse version |{com_args}'| ver_line={ver_line} failed {e}"
        else:
            story = f"Command failed: |{com_args}|"
    else:
        story = f"Command not found: |{exe}|"

    return (ok, path, ver, story)


def call_exec(
    com_args,
    shell=True,
    detach=False,
    timeout=None,
    verb=True,
    vok=100,
    verr_all=True,
    vpref="# Subproc: ",
    rep_time="Duration",
    ofname="",
):
    """Call subprocess executable as if via command line.

    com_args    = Command line given as string or list; If list join as strings
    shell       = Flag to call with shell (so env vars good, etc); Usually needed
    detach      = Flag to call exec in background; Don't wait for completion
    timeout     = timeout arg to limit subprocess
    verb        = Flag to report story; i.e. verbose mode
    vok         = (verbose) limit on output length if status is ok; 0 = no limit
    verr_all    = (verbose) flag to report all output on error
    vpref       = (verbose) prefix string for prints; '' or None for no prefix
    rep_time    = (verbose) string to report time; If None, not reported
    ofname      = Name for output file (vs stdout)

    Both stdout and stderr returned together as output (unless detach)
    Output is decoded to string with utf-8 and error='replace';
        This will replace (destroy) any unicode chars in output

    Return tuple (status, output, story)
    """
    status = True
    callout = out_story = ""
    story_lines = []

    # Args to string in case list
    if isinstance(com_args, list):
        com_args = " ".join([str(x) for x in com_args])
    if not isinstance(com_args, str):
        story = f"com_args not string or list; type: {0}".format(type(com_args))
        raise Exception(story)
    # Clean up multiple spaces between args (just better to report)
    com_args = " ".join(com_args.split())

    # Output story; Keep and maybe print
    OUTFILE, new_open = get_out_fh(ofname=ofname)
    story = f"{vpref}Call:   {com_args}"
    if verb:
        print_now(story, file=OUTFILE)
    story_lines.append(story)

    s_time = datetime.datetime.now(datetime.timezone.utc)
    # Call and catch output or problem
    try:
        # Detached cannot use run; Need direct call
        if detach:
            # Based on https://stackoverflow.com/questions/3516007/run-process-and-dont-wait
            # No outputs and no parent-child file handle sharing
            # Use shell true and com_args as str (OS dependent?... matter?)
            callproc = subprocess.Popen(
                com_args,
                shell=True,
                stdin=None,
                stdout=None,
                stderr=None,
                close_fds=True,
            )
        else:
            # stdout, stderr, text args to catch output together
            # See https://docs.python.org/3/library/subprocess.html#using-the-subprocess-module
            # However, text can fail to decode, so better to catch bytes here
            callproc = subprocess.run(
                com_args,
                shell=shell,
                timeout=timeout,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                # text=True,
                # encoding="utf-8",
            )
            # Decode stdout and stderr to str; Try to get before checking failure
            callout = callproc.stdout.decode("utf-8", errors="replace")
            # This will raise an exception if subprocess call failed
            callproc.check_returncode()

        out_story = "Success"

    except Exception as e:
        status = False
        out_story = f"Failed: {e}"
        # If verbose all error, remove output length limit
        if verr_all:
            vok = 0

    # Report story; Collect to return, maybe print here
    story1 = f"{vpref}Status: {out_story}"
    story_lines.append(story1)
    if vok > 0:
        story2 = f"{vpref}Output: (first {vok} char)\n{callout[:vok]}"
    else:
        story2 = f"{vpref}Output: (full)\n{callout}"
    story_lines.append(story2)
    t_str = get_time(stime=s_time, as_str=True)
    story3 = f"{vpref}{rep_time} {t_str}"
    story_lines.append(story3)
    if verb:
        print_now(story1, file=OUTFILE)
        print_now(story2, file=OUTFILE)
        if rep_time:
            print_now(story3, file=OUTFILE)

    if new_open:
        OUTFILE.close()

    return status, callout, "\n".join(story_lines)


def wrap_exe_com_args(command, line_continue=True, as_str=True):
    """Take command line arguments, clean up and 'wrap' one per line

    Return str or list
    """
    # Replace multiple spaces with single
    command = re.sub(" +", " ", command)
    # Split on command line option starts
    parts = command.split(" -")
    sing_lines = []
    for i, line in enumerate(parts):
        if i > 0:
            line = "-" + line
        if line_continue:
            trail = " \\" if i < (len(parts) - 1) else ""
            line += trail
        sing_lines.append(line)

    # Wrapped string or raw list
    if as_str:
        sing_lines = "\n".join(sing_lines)
    return sing_lines


def get_dir_file_lines(path, fname, verb=False, toxic=False):
    """Collect file content matching fname in/under given dir path

    Return dict with [path] = list of file lines
    """
    m_files = {}
    files = fnames_for_path(path=path, fname=fname)
    # Get conteent for matching files; Only non-comment, non-blank lines
    for f in files:
        lines = lines_from_fname(
            f, verb=verb, comment="#", strip=True, all_lines=False, toxic=toxic
        )
        if lines:
            m_files[f] = lines
    return m_files


def lines_from_fname(
    fname, verb=False, comment="#", strip=True, all_lines=True, max_keep=0, toxic=False
):
    """Attempt to parse given file and return list of lines

    fname = file
    verb = flag to report loading
    comment = char to strip lines past; None for full lines
    strip = flag to strip the end of input lines (i.e. rstrip)
    all_lines = flag to keep all lines, even if empty
    toxic = flag to crash if can't parse file

    Return list of file line (strings)
    """
    lines = []
    try:
        n = n_keep = 0
        with open(fname) as INFILE:
            for line in INFILE:
                n += 1
                if comment:
                    line = line.split(comment)[0]
                if strip:
                    line = line.rstrip()
                if not line:
                    if not all_lines:
                        continue
                n_keep += 1
                if max_keep and (n_keep > max_keep):
                    break
                lines.append(line)
        if verb:
            if all_lines:
                print(f"Read {len(lines)} lines from {fname}")
            else:
                print(f"Collected {len(lines)} lines from {fname} ({n} input lines)")
    except Exception as e:
        story = f"Failed to read from {fname}; {e}"
        if toxic:
            raise ValueError(story)
        if verb:
            print(story)
    return lines


def get_out_fh(ofname="", mode="a"):
    """Get output file handle; Input = filename, file handle, or nothing

    If filename is real and not already a file handle, open with mode
    If filename is a file handle, use that
    If fname = is nothing, use stdout

    Return tuple (File handle, Flag if newly opened)
    """
    # Default = stdout
    FH = sys.stdout
    new_open = False
    if ofname:
        # File handle already; use this
        if isinstance(ofname, io.IOBase):
            FH = ofname
        # Try to open file
        else:
            FH = open(ofname, mode)
            new_open = True
    return FH, new_open


def close_out_fh(fh, verb=True):
    """Close passed file handle, unless stdout

    verb = Flag to report story
    """
    if fh and (fh != sys.stdout):
        if verb:
            if fh.mode == "a":
                print("Updated file:", fh.name)
            else:
                print("New file:", fh.name)
        fh.close()


def file_header(
    ofile=None, mode="w", story=None, tblank=True, pref="# ", ver="", date=True
):
    """Write file header to file or string

    ofile = file or string; If string, attempt to open as file (mode)

    Return string
    """
    if not pref:
        pref = ""
    lines = []
    if story:
        lines.append(f"{pref}{story}")
    if ver:
        lines.append(f"{pref}{ver}")
    if date:
        # Get just date part of '<date> <time>'
        date = get_time(date=True, as_str=True).split()[0]
        lines.append(f"{pref}{date}")
    if tblank:
        lines.append(f"{pref}")
    string = "\n".join(lines)
    # If there's an output file, write to that
    if ofile:
        if isinstance(ofile, str):
            with open(ofile, mode) as OUTFILE:
                print(string, file=OUTFILE)
        else:
            print(string, file=ofile)
    return string


def dump_dict(
    dic,
    ofname="",
    pref="",
    sep="\t",
    kwidth=20,
    sort=True,
    head=True,
    otype=False,
    max_out=0,
):
    """Dump dictionary contents, to file or stdout"""
    OUTFILE, new_open = get_out_fh(ofname=ofname, mode="w")

    keys = dic.keys()
    if sort:
        keys = sorted(keys)
    if head:
        fbar = "-" * kwidth
        print(f"{fbar} {len(keys)} dict items {fbar}", file=OUTFILE)
    if kwidth:
        fmt_str = "{pref}{k:" + str(kwidth) + "}{sep}{v}"
    else:
        fmt_str = "{pref}{k}{sep}{v}"
    for i, k in enumerate(keys):
        v = dic[k]
        if otype:
            v = type(v)
        print(fmt_str.format(pref=pref, k=k, sep=sep, v=v), file=OUTFILE)
        if max_out and i >= max_out:
            break

    if new_open:
        close_out_fh(OUTFILE)


def dict_lower(dic, keys=True, vals=True):
    """Convert dict keys/values to lower case

    Return new dict
    """
    lc_dic = {}
    for k, v in dic.items():
        if keys and isinstance(k, str):
            k = k.lower()
        if vals and isinstance(v, str):
            v = v.lower()
        lc_dic[k] = v
    return lc_dic


def dict_format_vals(dic, ffunc=None):
    """Apply format to values in dictionary

    Return given dic
    """
    if not ffunc:
        ffunc = report_num_str()
    for k in dic.keys():
        dic[k] = ffunc()
    return dic


def dict_val_count(t_dict, val):
    """Return count of value in dict"""
    n = sum(v == val for v in t_dict.values())
    return n


def get_nested_dict_val(t_dict, key, sep="."):
    """Attempt to recursively get value from dict (of dicts...)"""
    parts = key.split(sep)
    t_dict = t_dict.get(parts[0])
    if t_dict:
        if len(parts) > 1:
            t_dict = get_nested_dict_val(t_dict, sep.join(parts[1:]))
    return t_dict


def sep_dic_by_dup_vals(t_dict):
    """Separate dict into duplicated and unique value parts

    Return tuple of dicts (unique, duplicate)
    """
    u_dict = {}
    d_dict = {}
    for k, v in t_dict.items():
        if dict_val_count(t_dict, v) > 1:
            d_dict[k] = v
        else:
            u_dict[k] = v
    return (u_dict, d_dict)


def dict_to_kv_dict_list(sdic, key="name", val="val", pretty=True):
    """Create list of single named {key: value} dicts from single starting dict

    This is one format for data passed to javascript via jinja2
    Also cleans up numbers for pretty printing

    Return list of dicts
    """
    dict_list = []
    for k, v in sdic.items():
        if pretty:
            # Parent dir is a string
            if not k.startswith("Parent"):
                v = report_num_str(v)

        kv_dict = {key: k, val: v}
        dict_list.append(kv_dict)
    return dict_list


def check_and_gzip_file(fname, pigz=True, comp_lev=4, verb=True):
    """If input file exists, try to compress

    pigz = flag to use 'pigz' command if available
    comp_lev = compression level flag; If < 1, use default
    """
    ok = True
    # If given filename list, recur for each item
    if isinstance(fname, list):
        for f in fname:
            if not check_and_gzip_file(f, pigz=pigz, comp_lev=comp_lev):
                ok = False
                break
    # Single file name string
    else:
        # Only do anything if file exists
        if check_infile(fname, verb=False, toxic=False):
            # pigz can directly replace gzip for compression
            #   (pigz = parallel gzip)
            have_pig, _, _ = check_exec("pigz")
            if pigz and have_pig:
                command = f"pigz --force {fname}"
            else:
                command = f"gzip --force {fname}"
            # Compression flag? (pigz allows level 0, but only >0 here)
            if comp_lev:
                command += f" -{comp_lev}"

            # Execute
            ok, _, ex_story = call_exec(command, verb=verb)
            if not ok:
                print(f"Subprocess failed: |{command}|")
                print(ex_story)
    return ok


def check_and_rm_file(fname, isdir=False, verb=True):
    """Check if file exists and try to delete

    fname = filename or list of filenames
    isdir = flag to (recursively) remove dir

    Return number of file(s) deleted; -1 if issue
    """
    n = 0
    # If given filename list, recur for each item
    if isinstance(fname, list):
        for f in fname:
            if check_and_rm_file(f, isdir=isdir) < 0:
                n = -1
                break
            n += 1
    # Single file name string
    else:
        assert isinstance(fname, str), f"Expected str got {type(fname)}"
        if check_infile(fname, verb=False, toxic=False):
            try:
                if isdir:
                    shutil.rmtree(fname)
                else:
                    os.remove(fname)
                n = 1
            except Exception:
                n = -1
    return n


def copy_file(src, dst, move=False, verb=True):
    """Copy (or move) file from src to dst

    src = source filename or list of filenames
    dst = destination name or dir

    Returns filename if successful; None if not
    """
    # If given filename list, recur for each item
    if isinstance(src, list):
        assert os.path.isdir(dst), f"To copy file list, dst must be dir {dst}"
        f_lis = []
        for f in src:
            name = copy_file(f, dst, move=move, verb=verb)
            if name:
                f_lis.append(name)
            else:
                f_lis = None
                break
        return f_lis

    # Single file name string
    assert isinstance(src, str), f"Expected source str got {type(src)}"
    com = "Move" if move else "Copy"
    o_str = None
    try:
        if move:
            o_str = shutil.move(src, dst)
        else:
            o_str = shutil.copy(src, dst)
        if verb:
            o_str = f"{com} {src} to {dst}"
    except Exception as e:
        if verb:
            print(f"Failed {com}: |{e}|")
    return o_str


def bump_logs(logname, verb=False, full_verb=False):
    """Bump version numbers of log file(s), like Linux /var/log/ files

    logname = first (version free, uncompress) filename

    Return number files considered
    """
    # Base part and dir path (with added slash for later)
    bname = os.path.basename(logname)
    dpath = os.path.dirname(logname) + "/"
    # All files that start with basename
    # List for '123.log' matches '123.log', 123.log.1' '123.log.2.gz' ...
    flist = fnames_for_path(path=dpath, start=bname, verb=verb)
    for fname in sorted(flist, reverse=True):
        # List of suffix parts beyond basename
        suf_parts = os.path.basename(fname).split(".")
        # Has .log, number and gz; e.g. 'something.log.4.gz'
        if len(suf_parts) > 3:
            ok, num = str_to_num(suf_parts[2], fail_val=0)
            new_name = dpath + f"{bname}.{num+1}.gz"
            copy_file(fname, new_name, move=True, verb=full_verb)
        # Only .log and number; e.g. 'something.log.1'
        elif len(suf_parts) > 2:
            ok, num = str_to_num(suf_parts[2], fail_val=0)
            new_name = dpath + f"{bname}.{num+1}"
            copy_file(fname, new_name, move=True, verb=full_verb)
            check_and_gzip_file(new_name, pigz=False, verb=False)
        # Only .log
        else:
            new_name = dpath + f"{bname}.1"
            copy_file(fname, new_name, move=True, verb=full_verb)
    return len(flist)


### ----------------------------------- fasta / fastq stuff  -------------------------------------


def get_fq2_from_fq1(fq1):
    """Try to get fq2 (R2) name from fq1 (R1) filename"""
    fq2 = ""
    parts = fq1.split("R1")
    # Make sure 'R1' is in name and only sub last instance
    #   e.g. '/pathR1/examR1.fq' >--> '/pathR1/examR2.fq'
    if len(parts) > 1:
        fq2 = "R1".join(parts[:-1]) + "R2" + parts[-1]
    return fq2


def fasta_kmer_sets(fpath, kmer_len=30, kmer_step=1):
    """Collect k-mers from fasta file

    kmer_len = word size

    Return tuple of sets, forward and reverse compliment
    """
    records = fasta_name_seq_dict(fpath)
    if records:
        kmer_set = set()
        rc_kmer_set = set()
        for seq in records.values():
            for i in range(0, len(seq) - kmer_len, kmer_step):
                kmer_set.add(seq[i : i + kmer_len])
            # Rev comp strand
            seq = reverse_complement(seq)
            for i in range(0, len(seq) - kmer_len, kmer_step):
                rc_kmer_set.add(seq[i : i + kmer_len])
        # tuple of sets
        return kmer_set, rc_kmer_set
    # Failed
    return None


def fasta_name_seq_dict(fpath, upper=True, one_name=True, dup_suf=""):
    """Get sequences from fasta file; fasta can be gzip

    upper    = flag to make all seq uppercase
    one_name = flag to keep first fasta header 'word' as name
    dup_suff = suffix for duplicates; Else duplicates not kept

    Return dictionary of records with key = name, val = seq
    """
    records = OrderedDict()
    tup_list = fasta_name_seq_list(fpath, upper=upper, one_name=one_name)
    # If keeping duplicates, make dict with names to use
    if dup_suf:
        # First get count of any dups
        name_dups = {}
        # Loop over [name, seq] tuples
        for name, _ in tup_list:
            if name in name_dups:
                name_dups[name] += 1
            else:
                name_dups[name] = 1
        # Make lists of names to use for each starting name
        for name, count in name_dups.items():
            if count > 1:
                n_list = []
                for i in range(count):
                    n_list.append(f"{name}{dup_suf}{i+1}")
                name_dups[name] = sorted(n_list, reverse=True)
            else:
                name_dups[name] = [name]
    # List to dict
    for name, seq in tup_list:
        # Raw name collision and fixing this?
        if dup_suf:
            name = name_dups[name].pop()
        records[name] = seq
    return records


def fasta_name_seq_list(fpath, upper=True, one_name=False):
    """Get sequences from fasta file; fasta can be gzip

    upper    = flag to make all seq uppercase
    one_name = flag to keep first fasta header 'word' as name

    Return list of (name, seq) tuples
    """
    if fpath.endswith("gz"):
        in_gz = True
        INFILE = gzip.open(fpath, "rb")
    else:
        in_gz = False
        INFILE = open(fpath)

    records = []
    name = seq = ""
    while True:
        if in_gz:
            line = INFILE.readline().decode("utf-8")
        else:
            line = INFILE.readline()
        if not line:
            break

        if line[0] == ">":
            if seq:
                records.append((name, seq))
                seq = ""
            if one_name:
                name = line[1:].strip().split()[0]
            else:
                name = line[1:].strip().replace(" ", "_")
        else:
            if upper:
                seq += line.strip().upper()
            else:
                seq += line.strip()

    if seq:
        records.append((name, seq))
    INFILE.close()
    return records


def load_fastq_recs(fpath, number, skip=0, head=True, seq=True, quality=True):
    """Read seqs from fastq / fastq.gz

    number = number to collect; If 0, collect all
    skip = number of initial records to skip
    head = flag to keep header lines too; Else empty list
    seq = flag to keep seq lines too; Else empty list
    quality = flag to keep quality lines too; Else empty list

    Return tuple of lists ([headers], [seqs], [qual])
    """
    if fpath.endswith(".gz"):
        infile = gzip.open(fpath, "rb")
        is_gz = True
    else:
        infile = open(fpath, "r")
        is_gz = False

    n = nkept = 0
    head_list = []
    seq_list = []
    qual_list = []
    while True:
        headline = infile.readline()
        seqline = infile.readline()
        # Don't use strand line
        _ = infile.readline()
        qualine = infile.readline()
        if not headline:
            break

        n += 1
        if n < skip:
            continue

        if is_gz:
            headline = headline.decode()
            seqline = seqline.decode()
            qualine = qualine.decode()

        if head:
            # Strip leading @ from header
            head_list.append(headline.strip()[1:])
        if seq:
            seq_list.append(seqline.strip())
        if quality:
            qual_list.append(qualine.strip())

        nkept += 1
        if number and (nkept >= number):
            break

    infile.close()
    return (head_list, seq_list, qual_list)


def get_sys_login_info():
    """Get system info

    Return tuple of node, machine, sysname, uname
    """
    sysname, node, _, _, machine = os.uname()
    uname = getpass.getuser()

    return node, machine, sysname, uname


def get_sys_command_line(sep=" ", wrap=False):
    """Return current system command line

    sep = separator for args
    wrap = flag to wrap on options starting with '-' and add '\' '\n' to string

    Return string
    """
    comline = sep.join(sys.argv)
    if wrap:
        # Wrap on ' -' so that line ends ' \' and new one starts '-'
        from_str = sep + "-"
        to_str = " \\\n-"
        comline = comline.replace(from_str, to_str)
    return comline


### -----------------------------    -----------------------------------


def parse_def_dict(fname, sep=",", keys=None, keep_list=False, verb=False):
    """Parse file simple 'definition' file into dict

    fname = file
    sep = separator between key value on lines
    keys = list (or comma,sep,string) of keywords to keep / check for
    keep_list = flag to keep duplicates in a list (vs just replace)
    verb = flag to report issue (i.e. keys don't match expectation)

    return dict of parsed lines
    """
    # Get only non-comment, non-empty lines from file
    lines = lines_from_fname(fname, verb=False, comment="#", all_lines=False)
    if not lines:
        return None
    # Collect dict from lines with at least two separated parts
    new_dict = {}
    for line in lines:
        parts = line.split(sep)
        if len(parts) > 1:
            new_key = parts[0]
            new_val = sep.join(parts[1:])
            # Keeping list of duplicated keys and have already?
            if keep_list and (new_key in new_dict):
                val = new_dict[new_key]
                # Append if alreayd have a list, else create new list
                if isinstance(val, list):
                    new_dict[new_key].append(new_val)
                else:
                    new_dict[new_key] = [val, new_val]
            # Just collect
            else:
                new_dict[new_key] = new_val

    # If given list of keys, make sure all / only these returned
    if keys:
        k_dict = {}
        # Make list; May be given comma delimited string
        if isinstance(keys, str):
            keys = keys.split(",")
        for k, v in new_dict.items():
            if k in keys:
                k_dict[k] = v
        # Didn't find all expected keys
        if len(k_dict) < len(keys):
            if verb:
                print(f"Failed to get expected {len(keys)} keys from file {fname}")
                print("Expecting:", keys)
                print("File has: ", list(new_dict.keys()))
            new_dict = None
        else:
            # Replace new with key dic
            new_dict = k_dict
    return new_dict


def parse_sample_def(fpath):
    """Parse sample definition from file

    Return tuple (name, wells, type)
    """
    name = wells = stype = ""
    with open(fpath) as infile:
        for line in infile:
            parts = line.strip().split(",")
            if parts[0] == "Name":
                name = parts[1]
            if parts[0] == "Wells":
                wells = parts[1]
            if parts[0] == "Type":
                stype = parts[1]
    return (name, wells, stype)


### ----------------------------------------- DGE -----------------------------------------


def load_anndata(fname):
    """Load scanpy anndata structure

    Return anndata or none
    """
    adata = anndata.read_h5ad(fname)
    if not isinstance(adata, scanpy.AnnData):
        adata = None
    return adata


def guess_parsebio_dge_fnames(dpath):
    """Attempt to collect DGE file names for given path

    dpath = path to DGE output dir

    Return tuple status, dict
    """
    mtx_set = set("count_matrix.mtx".split(","))
    var_set = set("all_genes.csv,target_genes.csv,all_guides.csv".split(","))
    obs_set = set("cell_metadata.csv".split(","))
    anndata_set = set("anndata.h5ad".split(","))

    new_dict = {"path": dpath, "mtx": "", "obs": "", "var": "", "anndata": ""}
    fpaths = fnames_for_path(dpath, depth=1)
    for one_path in fpaths:
        bname = os.path.basename(one_path)
        if bname in mtx_set:
            new_dict["mtx"] = one_path
        elif bname in var_set:
            new_dict["var"] = one_path
        elif bname in obs_set:
            new_dict["obs"] = one_path
        elif bname in anndata_set:
            new_dict["anndata"] = one_path
    # All good if all three parts for anndata
    ok = bool(new_dict["mtx"] and new_dict["var"] and new_dict["obs"])
    return ok, new_dict


def load_parsebio_dge(
    dpath,
    as_dict=False,
    anndata_fresh=False,
    gene_name=True,
    geno_pref=True,
    make_unique=True,
    verb=True,
    pref="# ",
):
    """
    Load DGE data from given path

    dpath = directory path with DGE files
    as_dict = flag to return dict with component structs (matrix, dataframes)
    anndata_fresh = flag to rebuild anndata from parts
    gene_name = flag to use gene name (symbol) as key, else use ID (accession)
    geno_pref = flag to prefix gene name with genome; If more than one genome, do this anyway
    make_unque = flag to make gene (var) names unique

    Return anndata or dict or None
    """
    # Dict with filenames
    ok, f_dict = guess_parsebio_dge_fnames(dpath)
    if not ok:
        if verb:
            print(f"{pref}Failed to find DGE files in {dpath}")
        return None

    # If already have anndata, don't want dict, and don't need fresh
    if f_dict["anndata"] and (not as_dict) and (not anndata_fresh):
        if verb:
            print(f"{pref}Loading {f_dict['anndata']}")
        return load_anndata(f_dict["anndata"])

    part_dict = {"mtx": None, "obs": None, "var": None}
    if verb:
        print(f"{pref}Loading anndata component files from {dpath}")

    # matrix
    part_dict["mtx"] = scipy.sparse.csr_matrix(sio.mmread(f_dict["mtx"]))

    # Cell metadata file currently like this (v0.9.7; Not targeted):
    #   bc_wells,sample,species,gene_count,tscp_count,read_count,bc1_well,bc2_well,bc3_well,bc1_wind,bc2_wind,bc3_wind
    #   70_71_89,Diseased_donor644,hg38,432,2329,5513,F10,F11,H5,70,71,89
    #   06_59_82,Healthy_donor630,hg38,347,1741,4574,A6,E11,G10,6,59,82
    #    ...
    df = read_csv(f_dict["obs"], index_col="bc_wells")
    # Discard all well-related columns; Make sure counts are int with zeros
    keep_cols = []
    count_cols = []
    for col in df.columns:
        if col.endswith("_count"):
            count_cols.append(col)
        if col.endswith("_well") or col.endswith("_wind"):
            continue
        keep_cols.append(col)
    df = df[keep_cols]
    # Set missing values
    df[count_cols] = df[count_cols].fillna(0).astype(int)
    df["sample"] = df["sample"].fillna("all-well")
    # IMPORTANT Set sample and species cols explictily to 'category'
    # If dtype cast to str, this breaks h5 file read / write
    cat_cols = ["sample", "species"]
    df[cat_cols] = df[cat_cols].astype("category")
    part_dict["obs"] = df

    # Gene file currently like this:
    #    ,gene_id,gene_name,genome
    #    0,ENSG00000000003,TSPAN6,hg38
    #    1,ENSG00000000419,DprintePM1,hg38
    #    ...
    # Note: 'gene_name' is actually gene 'symbol' (i.e. uppercase abbreviation)
    part_dict["var"] = read_csv(f_dict["var"])

    if as_dict:
        answer = part_dict
    else:
        answer = anndata_from_parts(
            part_dict["mtx"], part_dict["var"], part_dict["obs"]
        )
    return answer


def init_scanpy_params(verbosity=None, verb=True, pref="# "):
    """Initialize scanpy logging info"""
    if verbosity:
        scanpy.settings.verbosity = verbosity
        if verb:
            print(f"{pref}Set scanpy verbosity to {verbosity}")


def anndata_from_parts(
    count_mtx,
    genes_df,
    cells_df,
    gene_name=True,
    geno_pref=True,
    make_unique=True,
    verb=True,
    pref="# ",
):
    """Create scanpy anndata from component parts

    count_mtx = (sparse) count matrix
    genes_df = dataframe with genes, guides info (anndata 'var')
    cells_df = dataframe with cell info (anndata 'obs')
    gene_name = flag to use gene name (symbol) as key, else use ID (accession)
    geno_pref = flag to prefix gene name with genome; If more than one genome, do this anyway
    make_unque = flag to make var names unique

    returns scanpy anndata structure
    """
    if verb:
        print(f"{pref}Creating anndata with matrix {count_mtx.shape}")
    # Get AnnData structure
    adata = scanpy.AnnData(X=count_mtx, dtype=int)
    x_row, x_col = adata.shape

    # Cook up index as genome + gene(id or name)
    if gene_name:
        # Replace any missing name with ID
        genes_df["gene_name"] = genes_df["gene_name"].fillna(genes_df["gene_id"])
        genes_df["gene"] = genes_df["gene_name"]
    else:
        genes_df["gene"] = genes_df["gene_id"]

    # Append genome to 'gene' name if there are more than one
    if geno_pref or (len(genes_df["genome"].unique()) > 1):
        genes_df["gene"] = genes_df["gene"] + "_" + genes_df["genome"]

    genes_df.set_index("gene", inplace=True)
    n_genes = len(genes_df)

    # Cells
    n_cells = len(cells_df)

    # If transposed change to rows = cells, cols = genes
    if (x_row == n_cells) and (x_col == n_genes):
        if verb:
            print(
                f"{pref}Matrix {x_row} x {x_col}, cell {n_cells}, gene {n_genes} so transposing"
            )
        adata.transpose
        x_row, x_col = adata.shape

    if (x_row != n_cells) or (x_col != n_genes):
        story = f"Matrix {x_row} x {x_col} / cell {n_cells} / gene {n_genes} count mismatch "
        raise ValueError(story)

    # Cells = observations, Genes = variables
    adata.obs = cells_df
    adata.var = genes_df

    if make_unique:
        adata.var_names_make_unique()

    return adata


def load_tscp_assign_df(fname, fpath, sub_ind=-1, verb=True):
    """Load transcript assignment dataframe from file

    fname = filename
    fpath = file path
    sub_ind = index for (combine mode) sublibrary; If >= 0, add num to barcode

    Adds explicit 'well_ind' 1-based int column derived from well index field

    v0.9.4 file looks like:

        bc_wells,bc_ints,genome,gene,gene_name,count,exonic,rt_type,cell_barcode,polyN
        01_01_02,48_48_49,hg38,ENSG00000115163,CENPA,3,False,R,CTGCTTTG_GATAGACA_GCCACATA,CCAAGTTTCA
        01_01_02,48_48_49,hg38,ENSG00000163346,PBXIP1,2,True,R,CTGCTTTG_GATAGACA_GCCACATA,TTTTCAGGGT

    Return dataframe
    """
    tas_df = None

    csv_path = check_file_gz(fname, path=fpath, verb=verb)
    if not csv_path:
        if verb:
            print(f"Bad transcript assignement file (w/wo .gz) {fpath}/{fname}")
    else:
        try:
            # Pandas can read with/without .gz; compression='infer' is default
            # No index; barcode triples (i.e. 'bc_wells') should be just a column
            tas_df = read_csv(csv_path, index_col=None)

            # Well is from first barcode index
            tas_df["bc1_wind"] = tas_df["bc_wells"].apply(bcutils.cind_to_bc1_part)

            # If given index (combine mode), append 1-based to barcode seqs, indexes
            if sub_ind >= 0:
                tas_df = bcutils.add_suf_to_df_bci(tas_df, sub_ind + 1)
        except Exception:
            if verb:
                print(f"Problem loading / intializing {fpath}/{fname}")

    return tas_df


def label_list_subset(in_list, n_targ, min_keep=1):
    """Get subset of labels in given list, both fix min and proportional

    in_list = list with "labels" ... anything that can be put in a set
    n_targ = target number for final subsample
    min_keep = minimum number for each label; This preceeds proportional sampling

    Because of rounding, the number of returned indexes may differ from n_targ
    Also, if min_keep * number of lables > n_targ, result is bigger than n_targ

    Return set of list indexes to keep
    """
    n_all = len(in_list)
    # Start with all indexes valid
    keep_set = set(range(n_all))

    # Collect indexes for each label into sets
    lab_sets = {}
    for i, label in enumerate(in_list):
        if label in lab_sets:
            lab_sets[label].add(i)
        else:
            lab_sets[label] = {i}

    # Get total fixed count; These are at most min_keep per set
    # (Not simply n_labels * min_keep because some may be too small)
    n_fixed = 0
    for lset in lab_sets.values():
        n_fixed += min(min_keep, len(lset))

    # Total number beyond fixed set that can be subsampled
    n_ss_tot = n_all - n_fixed
    # Number of subsample to keep beyond fixed set to reach target
    n_ss_keep = n_targ - n_fixed

    if n_ss_keep > 0:
        prob_kill = 1 - n_ss_keep / n_ss_tot
        for _, lset in lab_sets.items():
            n_extra = len(lset) - min_keep
            if n_extra > 0:
                n_kill = round(n_extra * prob_kill)
                kill_lis = random.sample(sorted(lset), k=n_kill)
                kill_set = set(kill_lis)
                keep_set = keep_set - kill_set

    return keep_set
