#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Sample loading table functions

import numpy as np
import pandas as pd


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import plates
    import utils
else:
    from splitpipe import plates
    from splitpipe import utils


### ------------------------------ SampleLoadingTable --------------------------------

# 96 well plate defs
PLATE96_ROWS = list("ABCDEFGH")
PLATE96_COLS = [str(i + 1) for i in range(12)]
PLATE96_N_ROWS = len(PLATE96_ROWS)
PLATE96_N_COLS = len(PLATE96_COLS)


def parse_sampleloadingtable(fname, keep_dfs=False, pref="# "):
    """Parse sample loading table

    Return dict, string
    """
    sltab_info = {}
    out_lines = []
    # Dataframes for sheets
    sltab_dict, story = get_sltab_excel_dfs(fname)
    if not sltab_dict:
        return None, story

    # List of samples and df with wells
    samp_list, story = parse_sltab_sample_df(sltab_dict["sample"])
    if not samp_list:
        # print("<< return none", story)
        return None, story
    if story:
        out_lines.append(story)
    # Plate well info, specification strings
    wells_df, story = parse_sltab_plate_df(sltab_dict["plate"])
    if wells_df is None:
        # print("<< return none", story)
        return None, story
    if story:
        out_lines.append(story)

    # Set well specification string for each sample
    wspec_list = plates.all_samp_specs_for_96plate(wells_df)
    for s_dict, well_str in zip(samp_list, wspec_list):
        s_dict["wells"] = well_str

    # Save info
    sltab_info["filename"] = fname
    sltab_info["samples"] = samp_list
    sltab_info["num_samples"] = len(samp_list)
    sltab_info["num_wells"] = plates.num_snum_wells(wells_df, -1)
    if keep_dfs:
        sltab_info["wells_df"] = wells_df
        sltab_info["plate_df"] = sltab_dict["plate"]
        sltab_info["sample_df"] = sltab_dict["sample"]

    # Successful parse story
    out_lines.append(f"{pref}Samples from SampleLoadingTable: {fname}")
    for s_dict in samp_list:
        story = f"{pref} {s_dict['number']}, {s_dict['name']}, {s_dict['wells']}"
        out_lines.append(story)
    out_lines.append(f"{pref} --------------------------------------------")
    story = f"# Table has {len(samp_list)} samples in {sltab_info['num_wells']} wells:"
    out_lines.append(story)
    story = utils.report_df_str(wells_df, pref="#  ", index_name="Wells", row_num=False)
    out_lines.append(story)
    out_lines.append(f"{pref} --------------------------------------------")

    out_str = "\n".join(out_lines)
    return sltab_info, out_str


def get_sltab_excel_dfs(fname):
    """Attempt to parse SampleLoadingTable excel to dataframes per sheet

    Get dict with dataframes 'sample' and 'plate'

    Return, string
    """
    sheets = {}
    story = ""
    try:
        # 2024-02; Cannot get default "openpyxl" lib to work with "filter" excel
        # raw_sheets = pd.read_excel(fname, sheet_name=None)
        raw_sheets = pd.read_excel(fname, sheet_name=None, engine="calamine")
    except Exception as e:
        story = f"Failed to parse excel: {e}"
        story += f"\nMaybe check and disable 'Filter' in {fname}"
        return None, story

    # Find sample table and plate configuration sheets
    samp_df = plate_df = None
    for sname in raw_sheets.keys():
        if sname.lower().startswith("sample"):
            samp_df = raw_sheets[sname]
        if sname.lower().startswith("plate"):
            plate_df = raw_sheets[sname]
    if (samp_df is None) or (plate_df is None):
        story = "Could not find sample and plate sheets"
        story += f"Sheet names: {raw_sheets.keys()}"
        return None, story

    sheets["sample"] = samp_df
    sheets["plate"] = plate_df
    return sheets, story


def df_cell_str_match(df, query, row=-1, col=-1, startswith=True, as_num=False):
    """Check if given query string matches dataframe cell(s)

    Return match, row, col
    """
    if isinstance(df, pd.Series):
        df = pd.DataFrame(df)

    hit_str = hit_row = hit_col = None
    found = False

    query = str(query).lower()
    # One cell
    if (row >= 0) and (col >= 0):
        cell_str = str(df.iloc[row, col]).lower().strip()
        if as_num:
            good, _ = utils.str_to_num(cell_str)
            if good:
                found = True
        elif startswith and cell_str.startswith(query):
            found = True
        elif cell_str == query:
            found = True
        if found:
            hit_row = row
            hit_col = col
            hit_str = str(df.iloc[row, col])
    # One row
    elif row >= 0:
        for col in range(len(df.columns)):
            cell_str, _, _ = df_cell_str_match(
                df, query, row=row, col=col, startswith=startswith, as_num=as_num
            )
            if cell_str is not None:
                hit_row = row
                hit_col = col
                hit_str = str(df.iloc[row, col])
                break
    # One col
    elif col >= 0:
        for row in range(len(df.index)):
            cell_str, _, _ = df_cell_str_match(
                df, query, row=row, col=col, startswith=startswith, as_num=as_num
            )
            if cell_str is not None:
                hit_row = row
                hit_col = col
                hit_str = str(df.iloc[row, col])
                break
    # Whole thing
    else:
        for row in range(len(df.index)):
            cell_str, _, col = df_cell_str_match(
                df, query, row=row, startswith=startswith, as_num=as_num
            )
            if cell_str is not None:
                hit_row = row
                hit_col = col
                hit_str = str(df.iloc[row, col])
                break

    # Match a number? Convert string to number
    if (hit_str is not None) and as_num:
        _, hit_str = utils.str_to_num(hit_str)

    return hit_str, hit_row, hit_col


def parse_sltab_sample_df(df, pref="# WARNING "):
    """Parse 'Sample Table' (excel) sheet dataframe

    Get per-sample info from SampleSheet: Number (index), name, number wells
    Any issues reported in returned string

    Return list[dict], string
    """
    samp_data = []
    story_lines = []
    num_samp = 0

    # Get number of samples row in first col then number of samples
    query = "Number of Samples"
    hit_str, row, col = df_cell_str_match(df, query, col=0)
    if not hit_str:
        story = f"{pref}Could not find '{query}' if first col of sample sheet"
        return None, story
    num_samp, row, col = df_cell_str_match(df, query, row=row, as_num=True)
    if not num_samp:
        story = f"{pref}Could not get number of samples from row {row} of sample sheet"
        return None, story

    # Only look at remaining rows
    df = df.iloc[row + 1 :]

    # Collect columns for number, name, and number of wells
    # Find row with sample info, first col
    query = "Sample #"
    hit_str, row, col = df_cell_str_match(df, query, col=0)
    if not hit_str:
        story = f"{pref}Could not find '{query}' if first col of sample sheet"
        return None, story
    col_keys = {"number": col}

    query = "Sample name"
    hit_str, _, col = df_cell_str_match(df, query, row=row)
    if not hit_str:
        story = f"{pref}Could not find '{query}' in sample sheet"
        return None, story
    col_keys["name"] = col

    # Number of wells may differ v1 or v2 versions of the excel file
    query = "# wells"
    hit_str, _, col = df_cell_str_match(df, query, row=row)
    if hit_str:
        col_keys["num_wells"] = col
    query = "Number of Wells"
    hit_str, _, col = df_cell_str_match(df, query, row=row)
    if hit_str:
        col_keys["num_wells"] = col
    if "num_wells" not in col_keys:
        story = f"{pref}Could not find number of wells in sample sheet"
        return None, story

    # Collec info for samples, only look at remaining rows
    df = df.iloc[row + 1 :]
    n = 0
    for row in range(len(df.index)):
        name = df.iloc[row, col_keys["name"]]
        num = df.iloc[row, col_keys["number"]]
        n_wells = df.iloc[row, col_keys["num_wells"]]
        # If bogus line, name may be "nan", num should be ok, n_wells zero
        ok, n_wells = utils.str_to_num(n_wells)
        if (not ok) or (n_wells < 1):
            break
        samp_data.append({"number": num, "name": name, "n_wells": n_wells})
        n += 1

    # Check correct sample count
    if n < num_samp:
        story = f"{pref}Found {n} valid sample data rows, expecting {num_samp}"
        return None, story

    story = "\n".join(story_lines)
    return samp_data, story


def parse_sltab_plate_df(df):
    """Parse 'Plate Configuration' grid from excel sheet dataframe

    Dataframe should be dimensioned like plate, sample number in each well
    Any issues reported in returned string

    Return dataframe, string
    """
    # Plate start corner
    query = "A"
    hit_str, row, col = df_cell_str_match(df, query, startswith=False)
    if not hit_str:
        story = f"Could not find plate 'A' marker; Row {row}, col {col}"
        return None, story

    # New dataframe from table
    # No checking excel; Just grab rows x cols first
    wells_df = df.iloc[row : row + PLATE96_N_ROWS, col + 1 : col + 1 + PLATE96_N_COLS]
    wells_df.index = PLATE96_ROWS
    wells_df.columns = PLATE96_COLS

    # Remove non-real wells; i.e. drop rows with all missing numbers
    ok_rows = []
    for row in wells_df.index:
        if not wells_df.loc[row].isnull().all():
            ok_rows.append(row)

    # Legit rows, missing = 0 and all as int
    wells = wells_df.loc[ok_rows, :].replace(np.nan, 0).astype(int)
    return wells, ""
