#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import re
import numpy as np
import pandas as pd


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


# ---------------------------------------------------------------------------
class WellPlate:
    """Class for well plates; e.g. 96-well ..."""

    def __init__(self, n_rows=0, n_cols=0, name=""):
        # Default is 96 well
        self._n_rows = n_rows if n_rows > 0 else 8
        self._n_cols = n_cols if n_cols > 0 else 12
        self._n_rows = int(self._n_rows)
        self._n_cols = int(self._n_cols)
        if not name:
            name = f"{self._n_rows * self._n_cols}-well-plate"
        self._name = name

    def __del__(self):
        pass

    def __repr__(self):
        ostring = "WellPlate (plate specs)\n"
        ostring += f"Name:    {self.get_name()}\n"
        ostring += f"Dims:    {self.get_dims()}\n"
        return ostring

    def get_name(self):
        return self._name

    def num_rows(self):
        return self._n_rows

    def num_cols(self):
        return self._n_cols

    def get_dims(self):
        """Return number of rows, number of cols"""
        return self.num_rows(), self.num_cols()

    def num_wells(self):
        return self.num_rows() * self.num_cols()


# ----------------- Canned 96 well plate lable / index lists ----------------

PLATE96_ROWS = list("ABCDEFGH")
PLATE96_COLS = [str(i + 1) for i in range(12)]
PLATE96_N_ROWS = len(PLATE96_ROWS)
PLATE96_N_COLS = len(PLATE96_COLS)

# Mapping rows to 0-based coords; e.g. B >--> 1
PLATE96_ROW_TO_INT = dict(zip(PLATE96_ROWS, [i for i in range(PLATE96_N_ROWS)]))


def _plate96_well_index_maps():
    """Initialize mappings between well coords, lables and indexs"""
    # Mapping plate 0-base coords to well string; e.g. (1,3) >--> B4
    row_col_to_well = {}
    # Mapping plate 1-base int index to well string; e.g. 16 <--> B4
    wint_to_well = {}
    well_to_wint = {}

    wind = 0
    for i, row in enumerate(PLATE96_ROWS):
        for j, col in enumerate(PLATE96_COLS):
            well = f"{row}{col}"
            coords = (i, j)
            row_col_to_well[coords] = well
            # 1-based indexes
            wind += 1
            wint_to_well[wind] = well
            well_to_wint[well] = wind

    return row_col_to_well, wint_to_well, well_to_wint


# Create once as globals
(
    PLATE96_ROW_COL_TO_WELL,
    PLATE96_WINT_TO_WELL,
    PLATE96_WELL_TO_WINT,
) = _plate96_well_index_maps()


def wells_for_96plate(n_rows=0, as_well=True, as_int=False):
    """Get list of wells for 96 well plate, or fewer rows

    n_rows = number of rows (e.g. different kits have diff row count)
    as_well = flag to return wells, else index

    return list of strings
    """
    wells = []
    # How many plate rows?
    max_row = PLATE96_N_ROWS if n_rows == 0 else n_rows

    for r, row in enumerate(PLATE96_ROWS):
        if r >= max_row:
            break
        for col in range(PLATE96_N_COLS):
            wells.append(f"{row}{col+1}")

    # Possibly get indexes?
    if not as_well:
        wells = plate96_well_to_wind(wells, as_int=as_int)

    return wells


def rows_cols_for_96plate(n_rows=0, as_list=True):
    """Get lists of rows and col labels for 96 well plate, or fewer rows

    n_rows = number of rows (e.g. different kits have diff row count)
    as_list = flag to return lists, else counts

    return tuple of lists or ints (rows, cols)
    """
    # How many plate rows?
    max_row = PLATE96_N_ROWS if n_rows == 0 else n_rows

    rows = PLATE96_ROWS[:max_row]
    cols = PLATE96_COLS

    if not as_list:
        rows = len(rows)
        cols = len(cols)

    return rows, cols


# ----------------- Parse well specification to wells / indexs ----------------

# Well definition regex match patterns (For single plate)
# Single well, Colon:range, Dash-range; Get letter and number separate
WELL96_1_REGEX = re.compile(r"([ABCDEFGH])([1-9]|1[012])$")
WELL96_C_REGEX = re.compile(r"([ABCDEFGH])([1-9]|1[012]):([ABCDEFGH])([1-9]|1[012])$")
WELL96_D_REGEX = re.compile(r"([ABCDEFGH])([1-9]|1[012])-([ABCDEFGH])([1-9]|1[012])$")


def split_multi_plate96_wells(wspec):
    """Split well specification for multiple plates

    Return list of per-plate well specification strings
    """
    plate_wells = re.split("_+", wspec)
    return plate_wells


def join_multi_plate_well_str(wspec_list):
    """Joint multi-plate well specs; ['A1:D2', 'A1:B12'] >--> 'A1:D2__A1:B12'

    Return str
    """
    if isinstance(wspec_list, str):
        wspec_list = [wspec_list]
    return "__".join(wspec_list)


def parse_multi_plate96_wells(wspec, as_well=True, as_int=False):
    """Parse potentially multi-plate well specification into well indexes

    Wrapper for parse_96wells to handle multipe plate specs: e.g. "wells_wells"

    Return parse status, list of well specifications per plate, list of well
    indexes per plate, story if there are any parsing issues; Returns lists
    for wells and indexes even if there is just one plate

    Return tuple (status, list[str], list[list], story)
    """
    ok = True
    story = ""
    wlist_list = []
    # Split plate wells on (one or more) underscores
    plate_wells = split_multi_plate96_wells(wspec)
    for p, pwells in enumerate(plate_wells):
        ok, wlist, story = parse_plate96_wells(pwells, as_well=as_well, as_int=as_int)
        if not ok:
            # If multi-plate, add starting spec to story
            if len(plate_wells) > 1:
                story = f"{story}, plate {p+1} of {wspec}"
            break
        wlist_list.append(wlist)

    # Nothing if any issue
    if not ok:
        wlist_list = plate_wells = []

    return ok, plate_wells, wlist_list, story


def parse_plate96_wells(s, as_well=True, as_int=False):
    """Parse well specification string into well indexes

    s = well specification; e.g. 'A4,C5:D8' 'B4-B9', etc
    as_well = flag to return well str ('B1'), else 1-base index ('13')
    as_int = flag to return int index (if not well), else string

    Parses well specs into list of individual wells;
    How wells are denoted depends on flags; Default = well str

    Return tuple (status, list, story)
    """
    ok = True
    story = ""

    wells = np.arange(96, dtype=int).reshape(8, 12)
    sub_wells = []
    # Try will fail on non-match cases that are processed with checks
    try:
        blocks = s.upper().split(",")
        for b in blocks:
            if ":" in b:
                vals = _unpack_well_match(WELL96_C_REGEX.match(b), 4)
                s_row = PLATE96_ROW_TO_INT[vals[0]]
                s_col = vals[1] - 1
                e_row = PLATE96_ROW_TO_INT[vals[2]]
                e_col = vals[3]
                sub_wells += list(wells[s_row : e_row + 1, s_col:e_col].flatten())
            elif "-" in b:
                vals = _unpack_well_match(WELL96_D_REGEX.match(b), 4)
                s_row = PLATE96_ROW_TO_INT[vals[0]]
                s_col = vals[1] - 1
                e_row = PLATE96_ROW_TO_INT[vals[2]]
                e_col = vals[3] - 1
                sub_wells += list(
                    np.arange(wells[s_row, s_col], wells[e_row, e_col] + 1)
                )
            else:
                vals = _unpack_well_match(WELL96_1_REGEX.match(b), 2)
                s_row = PLATE96_ROW_TO_INT[vals[0]]
                s_col = vals[1] - 1
                sub_wells += [wells[s_row, s_col]]
        sub_wells = list(np.unique(sub_wells))
    except Exception:
        ok = False
        story = f"Problem well specification: '{b}' from '{s}'"
        pass

    # Possible conversion to various output formats if ok
    if ok:
        # Shift to 1-based indexes
        sub_wells = [w + 1 for w in sub_wells]
        if as_well:
            sub_wells = [PLATE96_WINT_TO_WELL[i] for i in sub_wells]
        else:
            if not as_int:
                sub_wells = [str(w) for w in sub_wells]
    else:
        sub_wells = []

    return (ok, sub_wells, story)


def _unpack_well_match(mat, num):
    """Unpack well regex match with two or four items [letter, well, [letter, well]]"""
    vals = []
    if num == 2:
        vals = [mat[1], int(mat[2])]
    elif num == 4:
        vals = [mat[1], int(mat[2]), mat[3], int(mat[4])]
    return vals


def plate96_cannonical_wspec(wspec, fatal=True):
    """Get cannonical well specification (via expand wells then contract)

    Return str
    """
    ok, well_list, _ = parse_plate96_wells(wspec)
    if ok:
        new_wspec = plate96_wind_to_well(well_list, minimize=True)
    else:
        if fatal:
            story = f"plate96_cannonical_wspec unparsable: |{wspec}|"
            raise ValueError(story)
    return new_wspec


def plate96_wind_to_well(wind, sort=True, minimize=False, as_list=False):
    """Get string from well indexes; e.g. [16, 17] >--> 'B4,B5'

    wind = well index or collection of well indexes (str or int)
    sort = flag to sort (simple, non-minimized) list
    minimize = flag to run minimization algorithm
    as_list = flag to return list of wells; Doesn't work with minimize

    If minimize
        Try to minimize well specification via block algorithm
    Else
        Generate simple comma,sep,list of wells: e.g 'A1,A2,A3,A4'

    return str or list
    """
    well = ""
    # If given a list, process each index
    if isinstance(wind, list) or isinstance(wind, set):
        # Make sure all int if given str (e.g. as well 'B1' or index '13')
        new_list = []
        for w in wind:
            if w in PLATE96_WELL_TO_WINT:
                new_list.append(PLATE96_WELL_TO_WINT[w])
            else:
                new_list.append(int(w))
        wind = new_list

        # Minimize spec string?
        if minimize:
            plate_df = get_96plate_df(val=0)
            plate_df = mark_96plate_wells(plate_df, wind, mark=1)
            well = samp_well_specs_for_96plate(plate_df, 1)
        # Just enumerate well index list
        else:
            well_list = []
            if sort:
                wind = sorted(list(wind))
            for i in wind:
                well_list.append(plate96_wind_to_well(i))
            if as_list:
                well = well_list
            else:
                well = ",".join(well_list)
    # Single well index
    else:
        # Make sure int
        wind = int(wind)
        well = PLATE96_WINT_TO_WELL[wind]
        if as_list:
            well = [well]
    return well


def plate96_well_to_wind(well, sort=False, as_int=False, as_list=True, na_val=-1):
    """Get well indexes from strings; e.g. 'B4,B5' >--> ['16', '17']

    sort = flat to sort
    as_int = flag to return int (else str)
    as_list = flat go return list; True if given a list
    na_val = value to use if well not found

    NOTE: For condensed well annotation (e.g. 'A1-B12'), use parse function()

    Return list or str or int
    """
    wind_list = []
    if isinstance(well, list) or isinstance(well, set):
        for w in well:
            wind_list.append(PLATE96_WELL_TO_WINT.get(w, na_val))
        if sort:
            wind_list.sort()
        if not as_int:
            wind_list = [str(w) for w in wind_list]
        wind = wind_list
    else:
        wind = PLATE96_WELL_TO_WINT.get(w, na_val)
        if not as_int:
            wind = str(wind)
        if as_list:
            wind = [wind]
    return wind


# ----------------- Plate matrix dataframe and reporting --------------------


def get_96plate_df(n_rows=0, val=0):
    """Get dataframe dimensioned as 96-well plate, or fewer rows

    n_rows = max number of rows for plate; default 0 = all
    val = value for matrix elements

    Return dataframe of int
    """
    row_labs, col_labs = rows_cols_for_96plate(n_rows=n_rows)
    plate_df = pd.DataFrame(
        np.zeros((len(row_labs), len(col_labs))), index=row_labs, columns=col_labs
    )
    # Possibly replace zero with value and cast to int
    plate_df = plate_df.replace(0, val).astype(int)
    return plate_df


def mark_96plate_wells(plate_df, wint_list, mark=1, non_mark=None, open_mark=None):
    """Mark plate dimensioned dataframe with well indexes

    plate_df = plate dataframe
    wind_list = list of 1-base int well indexes
    mark = value to mark sample wells
    non_mark = value to mark non-sample wells
    open_mark = optional value to flag so wells are only marked once; i.e. when open

    return dataframe
    """
    # Setting background?
    if non_mark is not None:
        plate_df[:] = non_mark

    # indexes one at a time
    for wind in wint_list:
        # 1-base so subtract 1
        r = (wind - 1) // PLATE96_N_COLS
        c = (wind - 1) % PLATE96_N_COLS
        row = PLATE96_ROWS[r]
        col = PLATE96_COLS[c]
        # Possibly ignore if already marked
        if open_mark is not None:
            if plate_df.at[row, col] != open_mark:
                continue
        plate_df.at[row, col] = mark

    return plate_df


def sample_marked_96plate_df(samp_list, n_rows=0):
    """Get dataframe plate marked with sample number

    samp_list = list of SpSamp objects or well-spec strings ['A1-B4','B5:C12'...]
    n_rows = max number of rows for plate; default 0 = all

    Return dataframe
    """
    # Get initial plate with all zero
    open_val = 0
    df = get_96plate_df(n_rows=n_rows, val=open_val)

    for i, samp in enumerate(samp_list):
        # Well spec
        if isinstance(samp, str):
            ok, wint_list, _ = parse_plate96_wells(samp, as_well=False, as_int=True)
            if not ok:
                return None
        # Assume sample object
        else:
            wint_list = samp.get_wind_list(as_int=True)

        # Mark dataframe; Only open positions get marked (first sample if overlap)
        df = mark_96plate_wells(df, wint_list, mark=i + 1, open_mark=open_val)
    return df


def report_plate_samples_str(spipe, pref="# Plate_samples    "):
    """Create text report of plate samples

    Return str
    """
    if spipe.get_mode("comb,mkref"):
        return ""

    str_rows = []

    str_rows.append("#" * 28 + "  Plate sample layout  " + "#" * 28)

    df = get_plate_sample_df(spipe)
    plate_str = utils.report_df_str(df, pref=pref, row_num=False, index_name="Wells")
    str_rows.append(plate_str)

    # Multi plates?
    if spipe.num_sample_plates() > 1:
        str_rows.append("#")
        str_rows.append("# Multi-plate samples; Only round1 plate shown")

    return "\n".join(str_rows)


def get_plate_sample_df(spipe, no_samp_val=0, meta=False):
    """Get dataframe with plate (round1) sample indexes per well

    no_samp_val = value to use for non-samples
    meta        = flag to include meta samples

    Return dataframe
    """
    well_list = []
    for i, samp in enumerate(spipe.get_samples(meta=meta)):
        # Ignore and all-well except if there are no others
        if samp.is_allwell() and (spipe.num_samples() > 1):
            continue
        well_list.append(samp.get_wspec(plate=1))

    # Get round 1 dims (bc_round is 1-based)
    n_rows, _ = spipe.get_bc_rows_cols(bc_round=1)
    df = sample_marked_96plate_df(well_list, n_rows=n_rows)

    return df


# ----------------- Plate (marked) to well spec processing -------------------


def all_samp_specs_for_96plate(wells_df):
    """Parse plate well dataframe into per-sample well specification strings

    wells_df = plate dimensioned dataframe with samples marked by int > 0

    Each sample (unique positive int in 'plate') gets a sample spec str (e.g. 'A5:D8')

    Return list of well spec strings
    """
    # Copy dataframe, as this gets set to zero to track well settings
    wells_df = wells_df.copy()
    wspec_list = []
    # Samples denoted by numbers in dataframe, 1,2.. max; Unset / consumed = zero
    samp_nums = sorted([x for x in set(wells_df.values.flatten()) if x > 0])
    for snum in samp_nums:
        well_spec = samp_well_specs_for_96plate(wells_df, snum)
        wspec_list.append(well_spec)

    return wspec_list


def samp_well_specs_for_96plate(wells_df, snum):
    """Parse plate well dataframe for single sample well specification string

    wells_df = plate dataframe with samples marked by int > 0
    snum = int number for sample to consider (i.e. which int to process)

    return well spec string
    """
    # Total wells for this sample
    snum_wells = num_snum_wells(wells_df, snum)
    well_def_list = []
    # Cook up 'blocks' of wells while any remain
    while snum_wells > 0:
        # Start coord for block
        st_coords = find_snum_well_start(wells_df, snum)
        if not st_coords:
            break
        # End coord and count; This modifies wells dataframe
        en_coords, num = mark_snum_well_block(wells_df, snum, st_coords, mark=0)
        if not en_coords:
            break
        snum_wells -= num
        # Save well block string
        wspec = well_block_str(st_coords, en_coords)
        well_def_list.append(wspec)

    # Final well spec for sample = comma,delim,list of well block strings
    well_spec = ",".join(well_def_list)
    return well_spec


def find_snum_well_start(wells_df, snum):
    """Find first well with given sample number

    Return coord tuple (row,col)
    """
    rows, cols = wells_df.shape
    for r in range(rows):
        for c in range(cols):
            if wells_df.iloc[r, c] == snum:
                return (r, c)
    return None


def count_snum_well_run(wells_df, snum, st_coords, down_col=False):
    """Find run of wells matching sample number, across or down from start (row,col)

    wells_df = plate dimensioned dataframe with samples marked by int
    snum = target sample int to count

    Return num
    """
    rows, cols = wells_df.shape
    st_row, st_col = st_coords
    if (not (0 <= st_row < rows)) or (not (0 <= st_col < cols)):
        return 0
    num = 0
    if down_col:
        for r in range(st_row, rows):
            if wells_df.iloc[r, st_col] == snum:
                num += 1
            else:
                break
    else:
        for c in range(st_col, cols):
            if wells_df.iloc[st_row, c] == snum:
                num += 1
            else:
                break
    return num


def mark_snum_well_block(wells_df, snum, st_coords, mark=0):
    """Mark out well block for sample number starting from plate coords

    wells_df = plate dimensioned dataframe with samples marked by int
    snum = target sample int to count

    Return end (row,col) tuple and number of wells marked
    """
    # Unpack bounds and starting row,col
    max_row, max_col = wells_df.shape
    st_row, st_col = st_coords
    # print(f"Dim {max_row}, {max_col}, snum={snum}, start = {st_coords}")
    # Init
    en_row = en_col = -1
    num = 0

    # Runs of snum down rows or cols
    col_run = count_snum_well_run(wells_df, snum, st_coords)
    row_run = count_snum_well_run(wells_df, snum, st_coords, down_col=True)
    # print(f" row_run {row_run}, col_run {col_run}")
    # Nothing either direction?
    if max(row_run, col_run) < 1:
        return None, 0

    # Go with the bigger of row or col to mark at a time
    if col_run > row_run:
        # print(f"Setting by row, {col_run} cols at a time")
        # Set by row, col_run at a time
        en_col = st_col + col_run
        row = st_row
        n = col_run
        while n == col_run:
            # print("Setting", row, ",", st_col, en_col)
            wells_df.iloc[row, st_col:en_col] = mark
            en_row = row
            row += 1
            if row >= max_row:
                break
            n = count_snum_well_run(wells_df, snum, (row, st_col), down_col=False)
            num += col_run
        # End col reduced to be inclusive (i.e. last == en_col)
        en_col -= 1
    else:
        # print(f"Setting by col, {row_run} rows at a time")
        # Set by col, row_run at a time
        en_row = st_row + row_run
        col = st_col
        n = row_run
        while n == row_run:
            # print("Setting", st_row, en_row, ",", col)
            wells_df.iloc[st_row:en_row, col] = mark
            en_col = col
            col += 1
            if col >= max_col:
                break
            n = count_snum_well_run(wells_df, snum, (st_row, col), down_col=True)
            num += row_run
        # End row reduced to be inclusive (i.e. last == en_row)
        en_row -= 1

    return (en_row, en_col), num


def num_snum_wells(wells_df, snum):
    """Return count of wells matching sample number

    wells_df = plate dimensioned dataframe with samples marked by int
    snum = target sample int to count
    """
    if snum < 0:
        num = (wells_df > 0).values.flatten().sum()
    else:
        num = (wells_df == snum).values.flatten().sum()
    return num


def well_block_str(st_coords, en_coords=None):
    """Get string to specify plate wells in block from st to en coords

    Return string
    """
    # Pair of coords, start (top left) to end (bottom right)
    if en_coords:
        st_row, st_col = st_coords
        en_row, en_col = en_coords
        # Start well alone via recur
        st_str = well_block_str(st_coords)
        # Same well so just one; recur
        if (st_row == en_row) and (st_col == en_col):
            wspec = st_str
        else:
            en_str = well_block_str(en_coords)
            # If one row, dash, else colon
            sep = "-" if st_row == en_row else ":"
            wspec = st_str + sep + en_str
    # Single coord
    else:
        wspec = PLATE96_ROW_COL_TO_WELL.get(st_coords, "")
    return wspec
