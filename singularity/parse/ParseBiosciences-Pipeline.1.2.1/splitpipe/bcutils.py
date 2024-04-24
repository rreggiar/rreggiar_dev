#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Barcode funcitons

import re
import os
import json
import numpy as np
import pandas as pd
from collections import defaultdict


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcinfo
    import plates
    import utils
else:
    from splitpipe import bcinfo
    from splitpipe import plates
    from splitpipe import utils


### ------------------------------ Barcode setup --------------------------------


def get_bc_info(bc_rounds, path, amp_seq, cell_thresh=0):
    """Load and init barcode set info structure

    Return BCInfo struct, story
    """
    bc_info = None
    story = ""
    try:
        bc_info = bcinfo.BCInfo(bc_rounds, path, amp_seq, cell_thresh=cell_thresh)
    except Exception as e:
        story = f"Failed to init barcodes: {e}"
    return bc_info, story


def get_bc_round_list(bc_args, as_int=False, names=False):
    """Get list of barcode rounds from (comline) arg list

    bc_args = list of list like [[1, 'v1'], ['bc2', 'v3'], ...
    as_int = flag to return list of numbers, not raw strings

    Return list (e.g. [1, 'bc2', ...])
    """
    rlist = []
    for bc_list in bc_args:
        if names:
            r = bc_list[1]
        else:
            r = bc_list[0]
            # Parse int or just collect?
            if as_int:
                ok, r = utils.str_to_num(r, regex=True)
                if not ok:
                    return []
        rlist.append(r)
    return rlist


### ------------------------------ Barcode IO --------------------------------
# Barcode file matching regex
BC_DATA_NAME_REGEX = re.compile(r".*/?bc_data_(\w+).csv")
BC_DICT_NAME_REGEX = re.compile(r".*/?bc_dict_(\w+).json")
BC_JSTR_NAME_REGEX = re.compile(r".*/?bc_data_(\w+).jstring")


def barcode_path():
    """Get path for installed barcodes"""
    try:
        path = os.path.dirname(__file__) + "/barcodes/"
    except Exception:
        path = "./barcodes/"
    return path


def barcode_filepaths(bc_path=None, jstring=False):
    """Get filenames for barcodes

    jstring = Flag to get (combined data and dict) encoded jstring files

    If jstring, return dict[bc] = file
    else dict[bc] = {'data_file': file, 'dict_file': file}

    Return dict with [name] = dict of file paths
    """
    if not bc_path:
        bc_path = barcode_path()

    new_dict = {}
    # Legit barcode sets are either jstring or have both csv data and json dict files
    if jstring:
        jstr_files = utils.fnames_for_path(bc_path, ext="jstring", depth=0)
        new_dict = utils.fname_list_to_dict(
            jstr_files, regex=BC_JSTR_NAME_REGEX, regex_key=1
        )
    else:
        csv_files = utils.fnames_for_path(bc_path, ext="csv", depth=0)
        json_files = utils.fnames_for_path(bc_path, ext="json", depth=0)
        csv_dict = utils.fname_list_to_dict(
            csv_files, regex=BC_DATA_NAME_REGEX, regex_key=1
        )
        json_dict = utils.fname_list_to_dict(
            json_files, regex=BC_DICT_NAME_REGEX, regex_key=1
        )
        # Return dict has names for both files
        for name in sorted(csv_dict.keys()):
            if name in json_dict:
                new_dict[name] = {
                    "data_file": csv_dict[name],
                    "dict_file": json_dict[name],
                }
    return new_dict


def read_bc_dict(fname):
    """Load barcode edit dict

    Return dict
    """
    with open(fname, "r") as INFILE:
        bc_dict = json.load(INFILE)
        if not bc_dict:
            print(f"Failed json.load for {fname}")
            return None

    # Top level has int keys and holds default dicts
    new_dict = {}
    for k, _ in bc_dict.items():
        new_dict[int(k)] = defaultdict(list, bc_dict[k])
    return new_dict


def read_bc_data_csv(fname):
    """Load barcode dataframe, check / update cols

    Return dataframe
    """
    # Don't want any index; Will set after checking all cols
    df = utils.read_csv(fname, index_col=None)

    cols = set("bci,sequence,uid,well,stype".split(","))
    if not set(cols).issubset(df.columns):
        print(f"Barcode data file problem {fname}")
        print(f"Expected cols {sorted(cols)}")
        print(f"Loaded cols   {sorted(df.columns)}")
        return None

    if not df["bci"].is_unique:
        print(f"Barcode data file problem {fname}")
        print("Barcode indexes ('bci') are not unique")
        return None

    df.set_index("bci", inplace=True)
    df.index.name = "bci"

    # Make sure plate col is here, as string
    if "plate" not in df.columns:
        df["plate"] = 1
    df["plate"] = df["plate"].astype(str)

    # Combined well + plate + stype have to be unique
    max_dup = max(df.groupby("well,plate,stype".split(",")).size())
    if max_dup > 1:
        print(f"Barcode data file problem {fname}")
        print("Well + stype + plate combinations are not unique")
        return None

    # Expand col collection with well-derived ones
    df = add_well_cols_to_bc_df(df)

    return df


### --------------------- Barcode dataframe functions ------------------------


def get_bc_well_list(df, plate=True, unique=True):
    """Get list of (single) well strings (e.g. A4,D10) from bc dataframe

    plate = flag to get the well + plate version (unique well)
    unique = flag to get unique list; Else full number of rows

    Return list (of str)
    """
    # plate-well names are unique across plates; No prefix if only one plate
    col = "well_plate" if plate else "well"
    if unique:
        wlist = list(df[col].unique())
    else:
        wlist = list(df[col])
    return wlist


def get_bc_num_wells(df):
    """Get number of wells for bc dataframe"""
    wlist = get_bc_well_list(df, unique=True)
    return len(wlist)


def get_bc_num_plates(df):
    """Get number of plates for bc dataframe"""
    return len(df["plate"].unique())


def get_bc_bci_list(df):
    """Get list of barcode indexes from bc dataframe"""
    nlist = list(df.index.values)
    return nlist


def get_bc_wind_list(df, zero_pad=True, unique=True):
    """Get list of well indexes from bc dataframe

    zero_pad = flag to get zero-padded index
    unique = flag to get unique list; Else full number of rows

    Return list (of str)
    """
    col = "wind_zpad" if zero_pad else "wind"
    if unique:
        wlist = list(df[col].unique())
    else:
        wlist = list(df[col])
    return wlist


def get_bc_rows_cols(df, as_list=False):
    """Get row and col for plate wells; Counts or lists

    as_list = Flag to return lists of rows, cols rather than just dims

    Return tuple (rows,cols) as int,int or list,list
    """
    rows = list(df["well_row"].unique())
    cols = list(df["well_col"].unique())
    # Make sure cols are sorted numeric
    cols = [str(c) for c in sorted([int(c) for c in cols])]
    # Returning just counts?
    if not as_list:
        rows = len(rows)
        cols = len(cols)
    return rows, cols


def get_bc_df_story(df):
    """Get story for given barcode dataframe

    Return string
    """
    n_bc = len(df)
    n_well = get_bc_num_wells(df)
    n_plate = get_bc_num_plates(df)
    rows, cols = get_bc_rows_cols(df, as_list=True)
    # Format so different ones align somewhat (e.g. for bc_list option)
    story = f"{n_bc:3d} barcodes"
    if n_plate > 1:
        story += f", {n_plate:3d} plates"
    story += (
        f", {n_well:3d} wells; rows {rows[0]}-{rows[-1]}, cols {cols[0]}-{cols[-1]}"
    )
    return story


def add_well_cols_to_bc_df(df):
    """Add well-derived columns to barcode dataframe

    Return (same) dataframe
    """
    # Split (string) row and col parts from well
    df["well_row"] = df["well"].apply(lambda s: s[0])
    df["well_col"] = df["well"].apply(lambda s: s[1:])

    n_plate = get_bc_num_plates(df)
    # Add plate + well, which is unique if multi plates; Prefix if multi plate
    if n_plate > 1:
        df["well_plate"] = "p" + df["plate"].astype(str) + "_" + df["well"]
        df["well_row"] = "p" + df["plate"].astype(str) + "_" + df["well_row"]
    else:
        df["well_plate"] = df["well"]

    # Mapping wells to 1-based indexes; Unique for wells across plates
    wells = plates.wells_for_96plate()
    wints = plates.wells_for_96plate(as_well=False, as_int=True)

    well_to_wint = dict(zip(wells, wints))

    # Collect well index; Some wells may be missing in df, so explicit loop
    wint_list = []
    for plate in df["plate"].unique():
        plate_off = 96 * (int(plate) - 1)
        for well in df[df["plate"] == plate]["well"]:
            wint = well_to_wint[well] + plate_off
            wint_list.append(wint)
    df["wind"] = wint_list

    # Zero padded string version
    if max(df["wind"]) > 99:
        df["wind_zpad"] = df["wind"].apply(lambda n: f"{n:03d}")
    else:
        df["wind_zpad"] = df["wind"].apply(lambda n: f"{n:02d}")

    # Well index is used as a string
    df["wind"] = df["wind"].astype(str)

    return df


### ------------------------- Read amplicon structure -------------------------
# For parsing amplicon representation into BC and polyN start positions
DNA_BASES = set(list("ACGT"))


def bc_amp_part_index(b):
    """Part index for given character

    If normal DNA base, -1
    If 'N' then = 0 for polyN
    If int, then int N for barcode-N
    """
    if b in DNA_BASES:
        index = -1
    else:
        if b == "N":
            index = 0
        else:
            index = int(b)
    return index


def get_amp_part_bounds(seq):
    """Get start and end coords (for slice) for non-base parts in seq

    Assumes amplicon like below, with BC as digit, (optional) polyN as N
    bc_amp_seq = 'NNNNNNNNNN33333333GTGGCCGATGTTTCGCATCGGCGTACGACT22222222ATCCACGTGCTTGAGACTGTGG11111111'

    Result lists are indexed [0]=N; [1]=round1; [2]=round2; [3]=round3

    Return arrays of start and end coords (for string slicing)
    """
    non_dna = set(seq.upper()) - DNA_BASES
    num_parts = len(non_dna)
    # If no N given, bump size
    if "N" not in non_dna:
        num_parts += 1
    part_ends = [-1] * num_parts
    part_starts = [-1] * num_parts
    # Scan each char; Annotate non-base parts
    prev_pi = -1
    i = 0
    while i < len(seq):
        pi = bc_amp_part_index(seq[i])
        # Change of part
        if pi != prev_pi:
            if pi >= 0:
                part_starts[pi] = i
            if prev_pi >= 0:
                part_ends[prev_pi] = i
        prev_pi = pi
        i += 1

    # Last end coord
    if part_ends[prev_pi] < 0:
        part_ends[prev_pi] = i
    return part_starts, part_ends


### ------------------------------ Barcode correction --------------------------------


def get_perfect_bc_counts(fq_df, bc_info, drop_valid_cols=True):
    """Evaluate sampled barcode sequences for perfect barcode counts

    fq_df = dataframe with barcode seqs
    bc_info = BCInfo struct
    drop_valid_cols = flag to delete added valid bc cols

    Return tuple (dict-with-count[perfect-bc-combo], threshold value, stats dict)
    """
    assert isinstance(fq_df, pd.DataFrame), f"Expected DataFrame got {type(fq_df)}"

    # Make sure barcode subseqs are split out
    fq_df = fq2_df_split_bc(fq_df, bc_info.get_amp_seq())

    # Sets a bit faster for membership checking
    bc1_seq_set = bc_info.get_bc_seqs(1, as_set=True)
    bc2_seq_set = bc_info.get_bc_seqs(2, as_set=True)
    bc3_seq_set = bc_info.get_bc_seqs(3, as_set=True)

    # Boolean perfect match of subseq to seq sets
    fq_df["bc1_valid"] = fq_df["bc1"].apply(lambda s: s in bc1_seq_set)
    fq_df["bc2_valid"] = fq_df["bc2"].apply(lambda s: s in bc2_seq_set)
    fq_df["bc3_valid"] = fq_df["bc3"].apply(lambda s: s in bc3_seq_set)

    # Counts of (perfect) barcode combinations
    counts = (
        fq_df.query("bc1_valid & bc2_valid & bc3_valid")
        .groupby(["bc1", "bc2", "bc3"])
        .size()
        .sort_values(ascending=False)
    )

    reads_in_cells_thresh = bc_info.get_cell_thresh()
    calc = counts.iloc[
        abs(counts.cumsum() / counts.sum() - reads_in_cells_thresh).values.argmin()
    ]
    count_threshold = max(5, calc)

    # Stats for perfect sequences
    stat_dict = {
        "bc1_perfrac": round(fq_df["bc1_valid"].sum() / len(fq_df), 3),
        "bc2_perfrac": round(fq_df["bc2_valid"].sum() / len(fq_df), 3),
        "bc3_perfrac": round(fq_df["bc3_valid"].sum() / len(fq_df), 3),
    }

    # Clean up dataframe?
    if drop_valid_cols:
        drop_cols = "bc1_valid,bc2_valid,bc3_valid".split(",")
        fq_df.drop(drop_cols, axis=1, inplace=True)

    return counts.to_dict(), count_threshold, stat_dict


def correct_barcodes(bc1, bc2, bc3, counts, count_thresh, bc_dicts, max_d):
    """Correct barcode sequences using pre-calculated values

    Return barcodes, edit distance max, edit distance sum
    """
    bc1_matches, edit_dist1, found1 = get_min_edit_dists(bc1, bc_dicts[1], max_d)
    bc2_matches, edit_dist2, found2 = get_min_edit_dists(bc2, bc_dicts[2], max_d)
    bc3_matches, edit_dist3, found3 = get_min_edit_dists(bc3, bc_dicts[3], max_d)
    # Valid if all found
    if found1 and found2 and found3:
        ed_max = max(edit_dist1, edit_dist2, edit_dist3)
        ed_sum = edit_dist1 + edit_dist2 + edit_dist3
    else:
        return "", "", "", -1, -1

    bc1_fixed = bc2_fixed = bc3_fixed = ""
    # Perfect match == done
    if 0 == edit_dist1 == edit_dist2 == edit_dist3:
        bc1_fixed = bc1_matches[0]
        bc2_fixed = bc2_matches[0]
        bc3_fixed = bc3_matches[0]
    # Look for single over-threshold combination
    else:
        matches = 0
        for bc1_m in bc1_matches:
            for bc2_m in bc2_matches:
                for bc3_m in bc3_matches:
                    try:
                        cur_counts = counts[(bc1_m, bc2_m, bc3_m)]
                    except Exception:
                        cur_counts = 0
                    if cur_counts > count_thresh:
                        bc1_fixed = bc1_m
                        bc2_fixed = bc2_m
                        bc3_fixed = bc3_m
                        matches += 1
                        # More than one possibility == don't use
                        if matches > 1:
                            bc1_fixed = bc2_fixed = bc3_fixed = ""
                            break

    return bc1_fixed, bc2_fixed, bc3_fixed, ed_max, ed_sum


def get_min_edit_dists(bc, edit_dict, max_d):
    """Returns a list of nearest edit dist seqs

    Input 8nt barcode, edit_dist_dictionary

    Output <list of nearest edit distance seqs>, <edit dist>
    """
    bc_matches = edit_dict[0].get(bc, [])
    edit_dist = 0
    found = bool(bc_matches)
    while (not found) and (edit_dist < max_d):
        edit_dist += 1
        bc_matches = edit_dict[edit_dist].get(bc, [])
        if bc_matches:
            found = True
            break
        if edit_dist >= max_d:
            break
    return bc_matches, edit_dist, found


def fq2_df_split_bc(df, amp_seq, fresh=False):
    """Add split barcode parts to fq2 (sampled) dataframe

    Return dataframe
    """
    bc_starts, bc_ends = get_amp_part_bounds(amp_seq)
    col_set = set(df.columns)
    seqs = df["seq"]
    quals = df["qual"]
    # Each barcode, starting from zero (0 = poly-N)
    n = 0
    for st, en in zip(bc_starts, bc_ends):
        if (en - st) > 0:
            bc_col = f"bc{n}"
            if fresh or (bc_col not in col_set):
                df[bc_col] = seqs.str.slice(bc_starts[n], bc_ends[n])

            qc_col = f"qc{n}"
            if fresh or (qc_col not in col_set):
                df[qc_col] = quals.str.slice(bc_starts[n], bc_ends[n])
        n += 1

    return df


def fq2_df_Q30_stats(fq2_df):
    """Add Q30 stats

    Return (updated given) dict
    """
    new_stats = {}
    if "qc0" in fq2_df.columns:
        new_stats["polyN_Q30"] = round(np.mean(fq2_df["qc0"].apply(seq_qual_score)), 3)
    new_stats["bc1_Q30"] = round(np.mean(fq2_df["qc1"].apply(seq_qual_score)), 3)
    new_stats["bc2_Q30"] = round(np.mean(fq2_df["qc2"].apply(seq_qual_score)), 3)
    new_stats["bc3_Q30"] = round(np.mean(fq2_df["qc3"].apply(seq_qual_score)), 3)
    return new_stats


def seq_qual_score(qual):
    """Get Q30 score for list of quality srings"""
    # Convert seq quality string into mean Q score
    # https://www.biostars.org/p/9463767/
    #
    # Pipeline = mean fraction "good" (Q30) bases
    m = np.mean([ord(c) > 62 for c in qual])
    return m


### ------------------------------ Barcode triples --------------------------------


def sep_df_bci_cols(df, bc_info=None, bc_wind=False, bc_well=False, which=None):
    """Separate dataframe composite barcode string into new columns

    df = cell dataframe with barcode well indexes 'bc_wells' as string triples '10_48_08'
    bc_info = BCInfo struct; Need for wells
    bc_wind = Flag to add well indexes (e.g 10, 48, 08)
    bc_well = Flag to add well names (e.g. A10, D12, A8)
    which = rounds to add; Single or list; If none, all rounds added

    Return dataframe
    """

    # If no list or set given, set all
    if not which:
        which = [1, 2, 3]
    if not isinstance(which, list):
        which = [which]

    # Get source of well info; Index or column
    if df.index.name == "bc_wells":
        bcw_df = pd.Series(df.index, index=df.index)
    else:
        bcw_df = df["bc_wells"]

    wind_cols = []
    # Well indexes; Always need (e.g. for well), maybe don't want to keep
    if bc_wind or bc_well:
        if 1 in which:
            df["bc1_wind"] = bcw_df.apply(cind_to_bc1_part)
            wind_cols.append("bc1_wind")
        if 2 in which:
            df["bc2_wind"] = bcw_df.apply(cind_to_bc2_part)
            wind_cols.append("bc2_wind")
        if 3 in which:
            df["bc3_wind"] = bcw_df.apply(cind_to_bc3_part)
            wind_cols.append("bc3_wind")

    # Well labels Barcode indexes
    if bc_well:
        assert isinstance(
            bc_info, bcinfo.BCInfo
        ), f"Expected BCInfo got {type(bc_info)}"
        if 1 in which:
            ww_map = bc_info.get_windzp_to_well(1)
            df["bc1_well"] = df["bc1_wind"].apply(lambda k: ww_map[k])
        if 2 in which:
            ww_map = bc_info.get_windzp_to_well(2)
            df["bc2_well"] = df["bc2_wind"].apply(lambda k: ww_map[k])
        if 3 in which:
            ww_map = bc_info.get_windzp_to_well(3)
            df["bc3_well"] = df["bc3_wind"].apply(lambda k: ww_map[k])

    # Clean up. Drop well indexes or reformat so no leading zero (but still string)
    if not bc_wind:
        df.drop(columns=wind_cols, inplace=True)
    else:
        for col in wind_cols:
            df[col] = df[col].apply(lambda s: f"{int(s):d}")

    return df


# Composite index (e.g. well, barcode): b1_b2_b3__sub
def cind_to_bc1_part(r):
    return r.split("__")[0].split("_")[0]


def cind_to_bc2_part(r):
    return r.split("__")[0].split("_")[1]


def cind_to_bc3_part(r):
    return r.split("__")[0].split("_")[2]


def cind_to_sublib(r):
    slib = ""
    parts = r.split("__")
    if len(parts) > 1:
        slib = parts[1]
    return slib


# v0.9.7 header line as: w1_w2_w3__pt__i1_i2_i3__bc1_bc2_bc3__polyN__dstamp_sindex
# Combined well indexes, primer type, bc indexes, bc seqs, poly N, data stamp and sequencing index
# Record names like:
#   20_15_07__T__20_63_55__TTCAACAT_CTTCACAA_GGGCGAAG__TTATTCTTTG__220814EO_CAGATC


# Well index parts
def fqrec_to_bc_winds(r):
    return r.split("__")[0]


def fqrec_to_bc1_wind(r):
    return fqrec_to_bc_winds(r).split("_")[0]


def fqrec_to_ptype(r):
    return r.split("__")[1]


# Barcode index parts
def fqrec_to_bc_bcis(r):
    return r.split("__")[2]


def fqrec_to_bc1_bci(r):
    return fqrec_to_bc_bcis(r).split("_")[0]


# Barcode sequence parts
def fqrec_to_bc_seqs(r):
    return r.split("__")[3]


def fqrec_to_bc1_seq(r):
    return r.split("__")[3].split("_")[0]


### ------------------------------ misc --------------------------------


def add_suf_to_df_bci(df, suffix, pref="__s", index=False):
    """Add index suffice to barcode cols in dataframe

    df = dataframe with barcode columns
    suffix = suffix to add to barcodes
    pref = string to append before appending suffix
    index = flag to apply suffix to index

    Return dataframe or nothing
    """
    # list of columns to modify
    cols = "bc_wells,bc_ints,cell_barcode".split(",")

    # Index
    if index:
        # If cols also, check index name is in list
        if cols:
            match = False
            for c in cols:
                if df.index.name == c:
                    match = True
                    break
            if not match:
                print(
                    f"Problem adding suffix to df; {cols} don't match index {df.index.name}"
                )
                return None

        df.index = df.index + pref + str(suffix)
    else:
        # Each listed column
        for c in cols:
            if c not in df.columns:
                print(f"Problem adding suffix to df; {c} not in columns {df.columns}")
                return None
            df[c] = df[c] + pref + str(suffix)
    return df


def write_bc_user_data(spipe, verb=True):
    """Write user-cleaned barcode data to output file

    Returns nothing
    """
    bc_info = spipe.get_bc_info()
    # Header part
    fpath = spipe.filepath("PF_BARCODE_DATA", None)
    kit_s, kit_n, kit_source = spipe.get_kit()
    story = f"barcode data; kit {kit_s} {kit_n}"
    ver = spipe.get_version()
    utils.file_header(ofile=fpath, story=story, ver=ver)

    # Cleaning barcodes to subset of type codes?
    kill_set = None
    if spipe.get_par_val("bc_save_data_clean", as_bool=True):
        kill_set = set("X")
    # Subset of cols to save
    out_cols = ["sequence", "uid", "well", "plate", "stype", "wind", "bc_round"]

    header = True
    name_list = bc_info.get_bc_name(as_list=True)
    for i, name in enumerate(name_list):
        if not name:
            continue

        bc_df = bc_info.get_bc_data(i)
        with open(fpath, "a") as OUTFILE:
            print(f"# Round {i} {name}", file=OUTFILE)

        # Clean up barcode output set?
        if kill_set:
            bc_df = bc_df[~bc_df["stype"].isin(kill_set)]

        story = utils.write_df(
            bc_df[out_cols],
            fpath,
            sep=",",
            header=header,
            index=True,
            mode="a",
            verb=False,
        )
        # Header for first barcodes only
        header = False

    if verb:
        spipe.report_run_story(f"Saved barcode data {fpath}")


def write_val_bc_counts(spipe, valid_counts, valid_2_counts=None, verb=True):
    """Write barcode data with valid bc counts

    valid_counts = list[rnd] of dicts[seq] with valid barcode counts
    valid_2_counts = dict[(seq1, seq2)] = count

    Returns nothing
    """
    bc_info = spipe.get_bc_info()

    # Collect dataframes with added counts col
    df_list = []
    for i in range(3):
        df = bc_info.get_bc_data(i + 1)
        df["valid_count"] = df["sequence"].map(valid_counts[i + 1])
        df_list.append(df)
    bc_df = pd.concat(df_list, axis=0)

    # Collapse into counts per well
    bc_df["round_well"] = bc_df["bc_round"].astype(str) + "_" + bc_df["well"]
    bc_df = bc_df.groupby("round_well")["valid_count"].sum().reset_index()
    # Rename col to simple count
    bc_df = bc_df.rename(columns={"valid_count": "count"})

    # Final as <round_well> <count>
    bc_df.set_index("round_well", inplace=True)
    spipe.write_df(bc_df, "PF_VAL_BC_DATA")

    # 2-bc counts?
    if valid_2_counts is not None:
        seq1_to_well = bc_info.get_seq_to_well(1)
        seq2_to_well = bc_info.get_seq_to_well(2)
        new_dict = {}
        for seq_tup, n in valid_2_counts.items():
            well1 = seq1_to_well[seq_tup[0]]
            well2 = seq2_to_well[seq_tup[1]]
            new_dict[well1 + "_" + well2] = n
        df = pd.DataFrame.from_dict(new_dict, orient="index", columns=["count"])
        df = df.sort_index()
        df.index.name = "round_1_2"
        spipe.write_df(df, "PF_VAL_2BC_DATA")


# -------------------------- BC kit scoring --------------------------------


def bc_count_kit_score(counts, bc_set, remain_frac=0.2, verb=False):
    """Calculate the score of fastq barcode seq counts for given barcode set

    counts = ranked barcode counts; e.g. series from df['bc1'].value_counts()
    bc_set = set of barcode sequences

    Return score
    """
    # Split counts into first and remaining parts
    # First part of counts = top N-bc counts
    f_tot = len(bc_set)
    f_counts = counts.iloc[:f_tot]
    f_num = len(set(f_counts.index) & bc_set)
    f_frac = f_num / f_tot
    f_mean = f_counts.mean()

    # Remaining part of counts via slice and threshold
    r_counts = counts.iloc[f_tot:]
    # Thresh = fraction of values in first part
    r_thresh = f_mean * remain_frac
    r_counts = r_counts[r_counts >= r_thresh]
    r_tot = len(r_counts.index)

    # None case
    if r_tot < 1:
        r_num = 0
        r_mean = r_frac = 0
    else:
        r_num = len(set(r_counts.index) & bc_set)
        r_frac = r_num / r_tot
        r_mean = r_counts.mean()

    if verb:
        print(f_tot, f_mean, r_thresh, sep="\t")
        print(f_num, f_frac, f_mean, sep="\t")
        print(r_num, r_frac, r_mean, sep="\t")

    score = f_frac - r_frac
    return score
