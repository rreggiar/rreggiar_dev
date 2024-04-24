#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import pandas as pd

from collections import Counter

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
    import utils
else:
    from splitpipe import bcutils
    from splitpipe import utils


# ---------------------------------------------------------------------------
class MapTarg:
    """Class for map-target info"""

    def __init__(self, spipe, name, fpath):
        """Initialize structure

        spipe = top-level splitpipe object
        name = name of target
        fpath = file path to sequence
        """
        self.spipe = spipe
        self.name = name
        # Check that supplied file is readable
        utils.check_infile(fpath, verb=True, toxic=True)
        self.fpath = fpath
        # Init empty fields
        self.kmer_set = None
        self.rc_kmer_set = None
        self.kmer_len = 0
        self.kmer_step = 0
        self.kmer_search_max = 0
        self.kmer_search_step = 0
        self.kmer_counts = Counter()
        self.rc_kmer_counts = Counter()

    def __del__(self):
        pass

    def __repr__(self):
        ostring = "MapTarg\n"
        ostring += f"Name: {self.name}\n"
        ostring += f"File: {self.fpath}\n"
        if self.kmer_set:
            ostring += f"kmer set sizes: {len(self.kmer_set)} {len(self.rc_kmer_set)}\n"
        ostring += f"kmer sets len {self.kmer_len} step {self.kmer_step}\n"
        ostring += (
            f"kmer search max {self.kmer_search_max} step {self.kmer_search_step}\n"
        )
        ostring += f"kmer counts {len(self.kmer_counts)} {len(self.rc_kmer_counts)}\n"
        return ostring

    def get_name(self):
        return self.name

    def get_path(self):
        return self.fpath

    def clear_counts(self):
        self.kmer_counts.clear()
        self.rc_kmer_counts.clear()

    def search_kmer_sets(self, seq, do_rc=False):
        """Search seq against kmer sets

        Return if kmer from seq is found in kmer set
        """
        found = False
        # Rev comp set or forward
        if do_rc:
            kmers = self.rc_kmer_set
        else:
            kmers = self.kmer_set
        # Scan into seq
        for i in range(0, self.kmer_search_max, self.kmer_search_step):
            try:
                sseq = seq[i : i + self.kmer_len]
                if sseq in kmers:
                    found = True
                    break
            except Exception:
                break
        return found


# ---------------------------------------------------------------------------

#### One time processing; Called in dge once globally (not per sample)


def proc_map_targ_counts(spipe):
    """Get map target counts (if any map targets are defined)

    Output saved to file

    Return status.
    """
    assert not spipe.get_mode("comb"), "proc_map_targ_counts not for comb mode"

    if spipe.num_map_targs() < 1:
        spipe.report_run_story2("No map targets to screen")
        return True

    # Dataframe with fastq sample; Don't retain for later (keep=False)
    ok, seqs_df = spipe.get_fastq_samp_df(0, keep=False)
    if not ok:
        spipe.set_problem("Failed to get fastq bc-corrected seq sample")
        return False
    # Lists of barcode indexes and seqs
    bci_list = list(seqs_df["head"].map(bcutils.fqrec_to_bc1_bci).astype(int))
    seq_list = list(seqs_df["seq"])

    story = f"Map target screening; {spipe.num_map_targs()} targets"
    spipe.report_run_story2(story)

    # List of all round 1 barcode indexes
    plate_bci_list = spipe.get_bc_bci_list(1)

    # Init counters
    # Collect total counts and lists for each map_targ for/rev
    mt_counts = {"count": dict(zip(plate_bci_list, [0] * len(plate_bci_list)))}
    c_cols = ["count"]
    for mapt in spipe.get_map_targ_list():
        name = mapt.get_name()
        f_key = name + "_count_for"
        r_key = name + "_count_rev"
        mt_counts[f_key] = dict(zip(plate_bci_list, [0] * len(plate_bci_list)))
        mt_counts[r_key] = dict(zip(plate_bci_list, [0] * len(plate_bci_list)))
        c_cols.append(f_key)
        c_cols.append(r_key)

    # Screen each map_targ
    for i, mapt in enumerate(spipe.get_map_targ_list()):
        name = mapt.get_name()
        story = f"Screening map_targ {name} ({i+1}) ({len(seq_list)} reads)"
        spipe.report_run_story2(story)
        f_key = name + "_count_for"
        r_key = name + "_count_rev"
        # Combinations of barcode index and sequence
        for bci, seq in zip(bci_list, seq_list):
            mt_counts["count"][bci] += 1
            if mapt.search_kmer_sets(seq):
                mt_counts[f_key][bci] += 1
            if mapt.search_kmer_sets(seq, do_rc=True):
                mt_counts[r_key][bci] += 1

    # Copy round1 bc data, bci = index, keep only well and seq type
    new_df = spipe.get_bc_df(1).copy()
    out_cols = ["well", "stype"]
    new_df = new_df[out_cols]
    # Dataframe with for and rev counts; Add combined 'any' counts
    c_df = pd.DataFrame.from_dict(mt_counts).astype(int)
    c_df.index = new_df.index
    out_cols.append("count")
    last_cols = []
    for mapt in spipe.get_map_targ_list():
        name = mapt.get_name()
        a_key = name + "_count_any"
        f_key = name + "_count_for"
        r_key = name + "_count_rev"
        c_df[a_key] = c_df[f_key] + c_df[r_key]
        out_cols.append(a_key)
        last_cols.append(f_key)
        last_cols.append(r_key)

    # Combine counts with well and type, then save with given col order
    new_df = pd.concat([new_df, c_df], axis=1)
    out_cols += last_cols

    spipe.write_df(new_df[out_cols], "PF_MAP_TARG_CT")
    return True
