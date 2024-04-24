#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Barcode funcitons


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
else:
    from splitpipe import bcutils


# Data filename templates
BC_DF_TEMP = "{path}/barcodes/bc_data_{bc_name}.csv"
BC_DIC_TEMP = "{path}/barcodes/bc_dict_{bc_name}.json"


# Default perfect barcode fraction factor (for bc correction)
DEF_CELL_THRESH = 0.92


# ---------------------------------------------------------------------------
class BCInfo:
    """Class for barcode set info"""

    def __init__(self, bc_rounds, path, amp_seq, cell_thresh=0):
        """Initialize structure

        bc_rounds = list as: [['bc1','n192_v4'], ['bc2','v1'], ['bc3','v1']]
        path =
        amp_seq = barcode-position-encoded read structure
        """

        self.bc_rounds = bc_rounds
        self.path = path
        # Threshold used for barcode correciton
        self.cell_thresh = cell_thresh if cell_thresh else DEF_CELL_THRESH

        # Amplicon info; Seq and start,end for each round (zero = polyN)
        self.amp_seq = amp_seq
        bc_starts, bc_ends = bcutils.get_amp_part_bounds(amp_seq)
        self.bc_starts = bc_starts
        self.bc_ends = bc_ends

        # Init lists with nothing (i.e. no round 0; rounds 1-3 have stuff)
        self.bc_dicts = [None]
        self.bc_data = [None]
        self.bc_seqs = [None]
        self.bc_names = [None]
        # field-to-field mapping
        # bci = unique index which maps 1:1 with unique seq
        # wind = well index which maps 1:1: with well string; zp = zero pad string
        # Multpiple bci/seq can map to single well/wind, so can't unambiguously
        #   map from well/wind to bci/seq)
        self.bc_seq_to_bci = [None]
        self.bc_seq_to_well = [None]
        self.bc_seq_to_type = [None]
        self.bc_seq_to_windzp = [None]
        self.bc_windzp_to_well = [None]
        self.bc_well_to_windzp = [None]
        self.bc_bci_to_well = [None]
        self.bc_bci_to_windzp = [None]
        self.bc_bci_to_type = [None]
        self.bc_bci_to_seq = [None]

        # Load barcode seqs and within-edit-distance dicts
        bc_round = 1
        for rlist in bc_rounds:
            bc_name = rlist[1]
            df_fname = BC_DF_TEMP.format(path=path, bc_name=bc_name)
            dict_fname = BC_DIC_TEMP.format(path=path, bc_name=bc_name)

            bc_df = bc_dict = None
            try:
                bc_df = bcutils.read_bc_data_csv(df_fname)
                bc_dict = bcutils.read_bc_dict(dict_fname)
            except Exception as e:
                story = f"Barcode {rlist} failed; file(s) {df_fname} {dict_fname}; {e}"
                raise Exception(story)

            # Add round to bc dataframe
            bc_df["bc_round"] = bc_round
            bc_round += 1

            self.bc_data.append(bc_df)
            self.bc_dicts.append(bc_dict)
            self.bc_seqs.append(list(bc_df["sequence"]))
            self.bc_names.append(bc_name)

            # Field-to-field mapping dicts
            bci_list = bcutils.get_bc_bci_list(bc_df)
            wind_list = bcutils.get_bc_wind_list(bc_df, unique=False)
            well_list = bcutils.get_bc_well_list(bc_df, unique=False)
            seq_list = list(bc_df["sequence"])
            type_list = list(bc_df["stype"])

            self.bc_seq_to_bci.append(dict(zip(seq_list, bci_list)))
            self.bc_seq_to_well.append(dict(zip(seq_list, well_list)))
            self.bc_seq_to_type.append(dict(zip(seq_list, type_list)))
            self.bc_seq_to_windzp.append(dict(zip(seq_list, wind_list)))
            self.bc_windzp_to_well.append(dict(zip(wind_list, well_list)))
            self.bc_well_to_windzp.append(dict(zip(well_list, wind_list)))
            self.bc_bci_to_seq.append(dict(zip(bci_list, seq_list)))
            self.bc_bci_to_well.append(dict(zip(bci_list, well_list)))
            self.bc_bci_to_windzp.append(dict(zip(bci_list, wind_list)))
            self.bc_bci_to_type.append(dict(zip(bci_list, type_list)))

    def __del__(self):
        pass

    def __repr__(self):
        ostring = "BCInfo (barcode set information)\n"
        for r, name in enumerate(self.get_bc_name(as_list=True)):
            if not name:
                continue
            df_story = self.get_bc_df_story(r)
            ostring += f"[{r}] {name: <8} {df_story}\n"
        return ostring

    def get_bc_df_story(self, which):
        """Get story for barcodes of given round"""
        df = self.get_bc_data(which)
        story = bcutils.get_bc_df_story(df)
        return story

    def num_rounds(self):
        """Return number of bc rounds"""
        # List is 1-based, so sub 1
        n = len(self.get_bc_name(as_list=True)) - 1
        return n

    def get_amp_seq(self):
        """Get bc-encoded read amplicon sequence"""
        return self.amp_seq

    def get_cell_thresh(self):
        """Get barcode correction cell threshold value"""
        return self.cell_thresh

    def get_bc_dict(self, which=1, as_list=False):
        """Get barcode round off-by-x-edit-dist dictionary for given bc round

        as_list = flag to return 1-based list of all dicts

        Return dict or list of dict
        """
        if as_list:
            bc_dict = self.bc_dicts
        else:
            bc_dict = self.bc_dicts[which]
        return bc_dict

    def get_bc_data(self, which):
        """Get barcode seq data dataframe for given bc round

        Return dataframe
        """
        df = self.bc_data[which]
        return df

    def get_bc_seqs(self, which, as_set=False):
        """Get barcode sequences for given bc round

        as_set = flag to return set of all names

        Return list or set
        """
        seqs = self.bc_seqs[which]
        if as_set:
            seqs = set(seqs)
        return seqs

    def get_bc_name(self, which=1, as_list=False):
        """Get barcode set name for given bc round

        as_list = flag to return 1-based list of all names

        Return str or list
        """
        if as_list:
            name = self.bc_names
        else:
            name = self.bc_names[which]
        return name

    def get_amp_starts_ends(self):
        """Get start and end coords for bc rounds in read amplicon seq

        Return two lists
        """
        starts = self.bc_starts
        ends = self.bc_ends
        return starts, ends

    def get_seq_to_bci(self, which):
        """Get mapping from seq to index for round

        Return dict
        """
        bc_map = self.bc_seq_to_bci[which]
        return bc_map

    def get_seq_to_well(self, which):
        """Get mapping from seq to well for round

        Return dict
        """
        bc_map = self.bc_seq_to_well[which]
        return bc_map

    def get_seq_to_type(self, which):
        """Get mapping from seq to (primer sequence) type for round

        Return dict
        """
        bc_map = self.bc_seq_to_type[which]
        return bc_map

    def get_seq_to_windzp(self, which):
        """Get mapping from seq to zero-pad windex for round

        Return dict
        """
        bc_map = self.bc_seq_to_windzp[which]
        return bc_map

    def get_windzp_to_well(self, which):
        """Get mapping from windex to well for round

        Return dict
        """
        bc_map = self.bc_windzp_to_well[which]
        return bc_map

    def get_well_to_windzp(self, which):
        """Get mapping from well to zero-pad windex for round

        Return dict
        """
        bc_map = self.bc_well_to_windzp[which]
        return bc_map

    def get_bci_to_seq(self, which):
        """Get mapping from index to sequence for round

        Return dict
        """
        bc_map = self.bc_bci_to_seq[which]
        return bc_map

    def get_bci_to_well(self, which):
        """Get mapping from index to well for round

        Return dict
        """
        bc_map = self.bc_bci_to_well[which]
        return bc_map

    def get_bci_to_windzp(self, which):
        """Get mapping from index to zero-pad windex for round

        Return dict
        """
        bc_map = self.bc_bci_to_windzp[which]
        return bc_map

    def get_bci_to_type(self, which):
        """Get mapping from index to (primer sequence) type for round

        Return dict
        """
        bc_map = self.bc_bci_to_type[which]
        return bc_map
