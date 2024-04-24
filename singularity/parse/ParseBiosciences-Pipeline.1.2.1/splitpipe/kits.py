#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Kit and chemistry related functions

import numpy as np
import pandas as pd


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
    import utils
else:
    from splitpipe import bcutils
    from splitpipe import utils


### -------------------------- Kit settings / parsing  ----------------------

KIT_CHEM_DEFS = """
    kit         chem  rows    cols    plates  bc1           bc2     bc3    ktype
    WT_mini     v1    1       12      1       n24_v4        v1      v1     normal
    WT          v1    4       12      1       v2            v1      v1     normal
    WT_mega     v1    8       12      1       n198_v5       v1      v1     normal

    WT_mini     v2    1       12      1       n24_v4        v1      v1     normal
    WT          v2    4       12      1       n99_v5        v1      v1     normal
    WT_mega     v2    8       12      1       n198_v5       v1      v1     normal

#    WT_mini     v3    1       12      1       n24_R1_v3_2   v1      R3_v3  normal
#    WT          v3    4       12      1       n99_R1_v3_2   v1      R3_v3  normal
#    WT_mega     v3    8       12      1       n198_R1_v3_2  v1      R3_v3  normal
    WT_mini     v3    1       12      1       n26_R1_v3_3   v1      R3_v3  normal
    WT          v3    4       12      1       n102_R1_v3_3  v1      R3_v3  normal
    WT_mega     v3    8       12      1       n208_R1_v3_3  v1      R3_v3  normal
"""


def _make_kit_chem_tab_df():
    """Get kit / chemistry table as dataframe"""
    head_row = []
    rows = []
    for line in KIT_CHEM_DEFS.split("\n"):
        parts = line.split("#")[0].split()
        if parts:
            if head_row:
                rows.append(parts)
            else:
                head_row = parts

    kit_df = pd.DataFrame(rows, columns=head_row)
    num_cols = "rows,cols,plates".split(",")
    kit_df[num_cols] = kit_df[num_cols].apply(pd.to_numeric, axis=1)
    kit_df["nwells"] = kit_df[num_cols].product(axis=1)
    return kit_df


# Create once as global
KIT_CHEM_TAB = _make_kit_chem_tab_df()


def get_kit_chem_ddata(kit, chem="v2"):
    """Get definition data for given kit and chemistry

    Return dict
    """
    new_dict = {}
    chem = parse_chemistry(chem)
    kit, _ = parse_kit(kit, chem=chem)
    query = f"kit == '{kit}' and chem == '{chem}'"
    drow = KIT_CHEM_TAB.query(query)
    if len(drow):
        new_dict = drow.iloc[0].to_dict()
    return new_dict


def parse_chemistry(chem, as_num=False):
    """Parse given chemistry to standard form

    Return str
    """
    answer = ""
    # Get number
    ok, n = utils.str_to_num(chem, regex=True)
    if ok:
        if as_num:
            answer = n
        else:
            chem = f"v{n}"
            if chem in KIT_CHEM_TAB["chem"].values:
                answer = chem
    return answer


def parse_kit(kit, chem=None, unique=True):
    """Attempt to parse given kit into cannonical form

    Returns cannonical named kit and story

    If name without chem is ambiguous and unique flag, fail

    Return tuple (str, str)
    """
    if not kit:
        return "", "<empty kit>"
    if is_custom_kit(kit):
        return "custom", ""

    # Make sure kit is string (e.g. in case given int)
    kit_s = str(kit)
    # If given chemistry, make sure it's good
    if chem:
        chem_s = parse_chemistry(chem)
        if not chem_s:
            return None, f"Bad chemistry ({chem}) given for parse_kit({kit})"
        chem = chem_s

    # If no chemistry, unique kit-chem not required
    if not chem:
        unique = False

    story = ""
    k_list = []
    df = KIT_CHEM_TAB

    # Try to parse number first; If not int, empty list
    k_list = kit_chem_from_int(kit_s, chem=chem)
    if len(k_list) < 1:
        # Lowercase and maybe dash not underscore (e.g. WT-mega >--> WT_mega)
        kit_s = kit_s.replace("-", "_").lower()
        df = df[df["kit"].str.lower() == kit_s]
        if chem:
            df = df[df["chem"] == str(chem)]
        k_list = [list(t) for t in df[["kit", "chem"]].values]

    # Different story for different number of matches
    if len(k_list) < 1:
        kit_s = ""
        story = f"Kit '{kit}' unrecognized"
    elif len(k_list) > 1:
        if unique:
            kit_s = ""
            story = f"Kit '{kit}' ambiguous; [kit,chem] = {k_list}"
        else:
            # First one, as all in list match except for chem
            kit_s, chem = k_list[0]
    else:
        kit_s, chem = k_list[0]

    return kit_s, story


def kit_chem_from_int(kint, chem=None, verb=True):
    """Get kit from int based on number of wells (e.g. 48 >--> WT), maybe chemistry

    chem = chemistry specification

    Return list of lists [[kit,chem], [kit,chem], ...]
    """
    if is_custom_kit(kint):
        chem = ""

    k_list = []
    # If given chemistry, make sure it's good
    if chem:
        chem_s = parse_chemistry(chem)
        if not chem_s:
            if verb:
                print(f"Bad chemistry ({chem}) given for kit_chem_from_int({kint})")
            return k_list
        chem = chem_s

    # Make sure starting with int; No regex, only simple int should work
    ok, n = utils.str_to_num(kint, regex=False)
    if ok:
        # Subset of rows by wells (int) and maybe chem (str)
        df = KIT_CHEM_TAB[KIT_CHEM_TAB["nwells"] == n]
        if chem:
            df = df[df["chem"] == chem]
        k_list = [list(t) for t in df[["kit", "chem"]].values]

    return k_list


def same_chemistry(chem1, chem2):
    """Check if two chemistries are the same"""
    # Parse to cannonical (str) form
    pchem_1 = parse_chemistry(chem1)
    pchem_2 = parse_chemistry(chem2)
    # Have to be something and agree
    return bool(pchem_1 and (pchem_1 == pchem_2))


def same_kit(kit1, kit2):
    """Check if two kits are the same"""
    # Parse to cannonical (str) form; Don't need story part
    pkit_1, _ = parse_kit(kit1)
    pkit_2, _ = parse_kit(kit2)
    # Have to have something and agree
    return bool(pkit_1 and (pkit_1 == pkit_2))


### ----------------------------- Kit functions -----------------------------


def is_custom_kit(kit):
    """Is given kit custom?"""
    custom = False
    if str(kit).lower().startswith("custom"):
        custom = True
    else:
        ok, n = utils.str_to_num(kit, regex=False)
        if ok and (n == 1):
            custom = True
    return custom


def kit_name_list(chem="v2", as_list=False, ktype="normal"):
    """Get list of kits

    as_list = flag to return list with [kit,chem] items
    ktype = filter for kit type column

    Return list[kit,] or list of [[kit, chem],]
    """
    # Filter table on any restrictions
    df = KIT_CHEM_TAB
    if chem:
        chem = parse_chemistry(chem)
        df = df[df["chem"] == str(chem)]
    if ktype:
        df = df[df["ktype"] == ktype]
    # Returning list of only kit or [kit, chem]
    if as_list:
        k_list = [list(t) for t in df[["kit", "chem"]].values]
    else:
        k_list = list(df["kit"].unique())
    return k_list


def kit_chem_list(kit):
    """Get list of chemistry for kit

    Return list[chem]
    """
    kit, _ = parse_kit(kit, chem=None)
    df = KIT_CHEM_TAB
    df = df[df["kit"] == kit]
    return list(df["chem"])


def kit_num_rows_cols_plates(kit, chem):
    """Number of rows, cols, plates for kit

    Return tuple(row, col, plate)
    """
    rows = cols = plates = 0
    dd_dict = get_kit_chem_ddata(kit, chem=chem)
    if dd_dict:
        rows = dd_dict["rows"]
        cols = dd_dict["cols"]
        plates = dd_dict["plates"]
    return rows, cols, plates


def kit_num_wells(kit, chem):
    """Get number of wells for kit

    Return int
    """
    if is_custom_kit(kit):
        return 1

    kit, _ = parse_kit(kit, chem=chem)
    # Filter table on any restrictions
    df = KIT_CHEM_TAB
    df = df[df["kit"] == kit]
    if chem:
        chem = parse_chemistry(chem)
        df = df[df["chem"] == str(chem)]
    # Only get number if unambiguous
    k_set = set(df["nwells"])
    nwell = list(k_set)[0] if len(k_set) == 1 else 0

    return nwell


def kit_bc_set_list(kit, chem, as_tup=False):
    """Get list of barcode set names for kit and chemistry

    kit = kit name
    chem = chemistry version
    as_tup = flag to return list of [bcX, name] items vs names

    Return list of barcode names per round
    """
    bc_list = []
    # Parse first
    chem = parse_chemistry(chem)
    kit_s, _ = parse_kit(kit, chem=chem)
    if kit_s and chem:
        df = KIT_CHEM_TAB
        row = df[(df["kit"] == kit_s) & (df["chem"] == str(chem))]
        cols = [f"bc{r}" for r in "123"]
        bc_list = list(df.loc[row.index, cols].values.flatten())
        # Tuple as [bcX, name]? This form matches 'bc_round_set' param
        if as_tup:
            new_list = []
            for col, bc in zip(cols, bc_list):
                new_list.append([col, bc])
            bc_list = new_list
    return bc_list


def kit_round1_bc_dict(chem):
    """Get dict mapping kit to name of round1 barcode set name

    Return dict[kit] = r1bc
    """
    new_dict = {}
    k_list = kit_name_list(as_list=True, chem=chem)
    for kit, chem in k_list:
        bc_list = kit_bc_set_list(kit, chem)
        new_dict[kit] = bc_list[0]
    return new_dict


def kit_ok_combine_sublib_count(spipe, from_par=True):
    """Check if number of sublibs is OK for given kit

    from_par = flag to get info from parameters (vs structs, etc)

    Return status
    """
    ok = False
    story = ""

    kit_s, kit_n, _ = spipe.get_kit()
    n_sublib = spipe.num_sublibs(from_par=from_par)
    # Need at least 2
    if n_sublib > 1:
        max_sublib = 100000
        # Custom kit = no checks
        if kit_n == 1:
            pass
        elif kit_n == 12:
            max_sublib = spipe.get_par_val("comb_max_sublibs_mini", as_int=True)
        elif kit_n == 48:
            max_sublib = spipe.get_par_val("comb_max_sublibs_wt", as_int=True)
        elif kit_n == 96:
            max_sublib = spipe.get_par_val("comb_max_sublibs_mega", as_int=True)
        else:
            raise ValueError(f"Unrecognized kit: {kit_s} {kit_n}")

        # Too many?
        if n_sublib > max_sublib:
            story = f"Too many sublibs = {n_sublib} for kit {kit_s}; Max = {max_sublib}"
            spipe.set_problem(story)
            ok = False
        else:
            ok = True

    return ok


def report_valid_kits(chem=None, ktype="normal", pad=True):
    """Report kit collection;

    If not given chemistry, report that too

    Returns nothing
    """
    # List of kit names, possibly limited to chem and ktype
    k_list = kit_name_list(as_list=False, chem=chem, ktype=ktype)

    if pad:
        print()
    print(f"There are {len(k_list)} installed kits:\n")
    for kit in k_list:
        n_wells = kit_num_wells(kit, chem)
        c_list = kit_chem_list(kit)
        chemistry = ", ".join(c_list)
        print(f"    {kit:12s} {n_wells} wells,   Chemistry: {chemistry}")
    if pad:
        print()


#### ------------------------ Kit barcode scoring ---------------------------


def guess_kit_from_bc(spipe):
    """Guess kit via score against barcode sets

    Return tuple (kit, best_score, list of all scores)
    """
    chem = spipe.get_chemistry()
    kit_bc_names = kit_round1_bc_dict(chem)
    story = f"Scoring fastq data against {len(kit_bc_names)} kits, chemistry {chem}"
    spipe.report_run_story(story)

    if not spipe.have_in_fastq(fq1=True, fq2=True):
        spipe.set_problem("Cannot score kit without fastq input")
        return (0, 0, None)

    # R2 sample barcode seq counts (doesn't have to be log10, but tuned for that)
    ok, fq2_df = spipe.get_fastq_samp_df("r2", split_bc=True)
    if not ok:
        spipe.set_problem("Failed to get fastq sample reads")
        return (0, 0, None)

    counts = np.log10(fq2_df["bc1"].value_counts())

    # Screen kits and barcode sets
    all_scores = {}
    bc_paths = bcutils.barcode_filepaths()
    remain_frac = spipe.get_par_val("kit_score_remain_frac", as_float=True)
    best_kit = best_score = None
    for kit, bc in kit_bc_names.items():
        # Only need perfect match sequences for scoring
        data_fname = bc_paths[bc]["data_file"]
        bc_df = bcutils.read_bc_data_csv(data_fname)
        bc_set = set(bc_df["sequence"])
        score = bcutils.bc_count_kit_score(counts, bc_set, remain_frac=remain_frac)
        spipe.report_run_story(f"  {kit} ({bc}) = {round(score, 3)}")
        if (not best_score) or (score > best_score):
            best_score = score
            best_kit = kit
        all_scores[kit] = [bc, score]
    spipe.report_run_story(f"Best scoring kit = {best_kit}, {round(best_score, 3)}")

    return best_kit, best_score, all_scores


def report_kit_score_str(spipe, pref="# Kit_score ", unset="<Not set>"):
    """Get lines reporting kit scoring specifications

    Return str
    """
    if spipe.get_mode("comb,mkref"):
        return ""

    str_rows = []

    ks_dict = spipe._get_kit_specs()
    if ks_dict and "kit_bc_scoring" in ks_dict:
        # Kit and kit score info
        kit_s, _, _ = spipe.get_kit()
        given_kit = unset
        score_dict = {}
        ks_info = ks_dict["kit_bc_scoring"]

        if ks_info:
            if ks_info["given_kit"]:
                given_kit = ks_info["given_kit"]
            score_dict = ks_info["kit_bc_scores"]
            if not score_dict:
                score_dict = {}

            # print(f"# Kit_score Processing: {ks_info['processing']}", file=OUTFILE)
            str_rows.append(f"{pref}Processing:     {ks_info['processing']}")

        # print(f"# Kit_score Given kit: {given_kit}", file=OUTFILE)
        str_rows.append(f"{pref}Given kit:      {given_kit}")

        # print(
        #   f"# Kit_score Total kits screened: {len(score_dict)}", file=OUTFILE
        # )
        str_rows.append(f"{pref}Total screened: {len(score_dict)}")
        for kit_name, (bc_name, score) in score_dict.items():
            score = str(round(float(score), 3))
            k_story = f"{kit_name} ({bc_name})"
            arrow = "    <---- using" if kit_name == kit_s else ""
            str_rows.append(f"{pref}Score {k_story:<24} = {score}{arrow}")

    return "\n".join(str_rows)
