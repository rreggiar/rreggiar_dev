#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import re
from itertools import product
from natsort import natsort_keygen


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import plates
    import spclass
    import spoutdir as spout
    import utils
else:
    from splitpipe import plates
    from splitpipe import spclass
    from splitpipe import spoutdir as spout
    from splitpipe import utils


# ---------------------------------------------------------------------------
class SpSamp:
    """Class for 'sample' info

    Single mode = data from single analysis, single plate of wells
    Combine mode = multiple analysis, (dict of) well info for each sublibrary

    Simple = given well definition
    Meta = combined wells from other samples
    """

    def __init__(
        self,
        name,
        wspec="",
        plate_dims=None,
        sublibs=None,
        samples=None,
        def_plate_dims=(8, 12),
    ):
        """Initialize structure

        name = Name of target
        wspec = Well specification (string; e.g. A2-A6,B5:C10)
        plate_dims = Plate dimension info; n_rows, tuple(s); Required for simple
        sublibs = Collection of sublibrary structs (combine mode only)
        samples = Collection of sample structs (i.e. meta parts if meta sample)
        def_plate_dims = default 96-well = (8, 12)

        Need either wells or samples (single mode) or sublibs (combine mode)
        """
        if sublibs is None:
            sublibs = []
        assert isinstance(sublibs, list), f"Expected list got {type(sublibs)}"

        if samples is None:
            samples = []
        assert isinstance(samples, list), f"Expected list got {type(samples)}"

        self.name = name
        self._sublibs = sublibs
        self._meta_parts = samples
        self._am_meta = False
        self._status = "good"
        self._stat_step = "start"

        # Collection of info for sublib structs; samples, well def strings and well lists
        # Lists universal to can handle collections same if single or combine mode
        # All (top-level) lists will be same length as N sublibs (or len = 1 if not combine)
        self._all_sub_samples = []
        # Well list lists (e.g. [sublib][plate] = ['A1', 'A2', 'B1']) for each sublib, each plate
        self._all_well_list = []
        # Well string list (e.g. [sublib] = 'A1-B12') one string for each sublib
        self._all_wspec_list = []

        # WellPlate structs with plate dims; List for multi-plate
        self._wplates = []
        n_plates = 1

        # Combine mode. If sublibraries, collect well / sample info for each sublib
        if self._sublibs:
            # If dict, get list
            if isinstance(self._sublibs, dict):
                self._sublibs = list(self._sublibs.values())

            # Collect same-name sample info for each sublib
            is_meta = False
            for slib in self._sublibs:
                assert isinstance(slib, spout.SpOut), f"Expected SpOut got {type(slib)}"
                # Get well and sub_well, or nothing / empty
                ssamp = slib.get_samples(one=name)
                # print("SSS", ssamp)
                if ssamp:
                    if ssamp.num_meta_parts():
                        is_meta = True
                    # Well spec; May be multi-plate
                    ss_wspec = ssamp.get_wspec()
                    # Wells as list of list for multi-plate; Take directly for sublib 0
                    ss_well_list = ssamp._all_well_list[0]
                    # Plate dims same for all sublibs; plate=0 gets all plates
                    plate_dims = ssamp.get_plate_dims(plate=0)
                else:
                    ss_wspec = ""
                    # Assumes max number of plates = 2; Not know here
                    ss_well_list = [] * 2

                # Set values for this sublib
                self._all_sub_samples.append(ssamp)
                self._all_wspec_list.append(ss_wspec)
                self._all_well_list.append(ss_well_list)
            # If any parts collected, this is meta
            self._am_meta = is_meta

        # Single mode. Just one sublib
        else:
            # Meta sample so combine all parts
            if self._meta_parts:
                compatible_meta_samples(self._meta_parts, fatal=True)

                # Plate dims same for all parts; Save first, plate=0 gets all plates
                plate_dims = self._meta_parts[0].get_plate_dims(plate=0)
                # Combine wells for each plate separately; Also well specs
                well_list_list = []
                wspec_list = []
                for p in self._meta_parts[0].get_plate_indexes():
                    # Union of all wells for plate
                    well_set = set()
                    for samp in self._meta_parts:
                        well_set.update(samp.get_well_list(plate=p))
                    # Well list and minimized spec from this
                    well_list = list(well_set)
                    well_list_list.append(well_list)
                    wspec = plates.plate96_wind_to_well(well_list, minimize=True)
                    wspec_list.append(wspec)
                # Single well spec in case multi-plate
                wspec = plates.join_multi_plate_well_str(wspec_list)
                # Set sample list via function
                self.set_meta_part_samples(self._meta_parts)
            # Simple, non-meta
            else:
                # Parse wells; Returns status, list[wspecs], list[well_list], issue story
                ok, _, well_list_list, story = plates.parse_multi_plate96_wells(wspec)
                if not ok:
                    raise ValueError(story)
                # Plate dimensions; If given an int, set as number of rows
                assert plate_dims is not None, "Simple SpSamp init requires plate dims"
                if isinstance(plate_dims, int):
                    plate_dims = (plate_dims, 12)

            # Save as lists even for single mode
            self._all_sub_samples = [None]
            self._all_well_list = [well_list_list]
            self._all_wspec_list = [wspec]
            # In case multi-plate; Number of plates = len well_list_list
            n_plates = len(well_list_list)

        self._set_plate_dims(plate_dims, n_plates=n_plates)

    def __del__(self):
        pass

    def __repr__(self):
        ostring = "SpSamp (sample definition)\n"
        ostring += f"Name:    {self.name}\n"
        raw_wells = self.get_def_str()
        ostring += f"Wells:   {raw_wells}\n"

        # Plates as list of tuple(row,col)
        dim_list = [wp.get_dims() for wp in self._wplates]
        ostring += f"Plate:   {dim_list}\n"

        type_story = self.get_type_str()
        ostring += f"Type:    {type_story}\n"
        if self.is_meta():
            ostring += f"Samples: {self.get_meta_parts(as_name=True)}\n"

        status, stat_step = self.get_status()
        ostring += f"Status:  {status} ({stat_step})\n"

        return ostring

    def _set_plate_dims(self, plate_dims, def_plate_dims=(8, 12), n_plates=1):
        """Set up well plate objects for each plate (e.g. multi-plate)

        plate_dims = list of plate dimension tuples
        def_plate_dims = default plate dims
        n_plates = number of plates; If fewer than this given, use defaults

        Return nothing
        """
        # If no plate dims, set default, and if tuple, put into list
        # Need at least 1 plate; Extend to 2 in general case
        if plate_dims is None:
            plate_dims = []
        elif isinstance(plate_dims, tuple):
            plate_dims = [plate_dims]
        assert isinstance(plate_dims, list), f"Expected list got {type(plate_dims)}"

        # Need defaults?
        while len(plate_dims) < n_plates:
            plate_dims.append(def_plate_dims)

        # Plates dims into structs and save
        wplate_list = []
        for dim_rc in plate_dims:
            n_rows, n_cols = dim_rc
            wplate_list.append(plates.WellPlate(n_rows=n_rows, n_cols=n_cols))
        self._wplates = wplate_list

    def get_name(self):
        return self.name

    def get_def_str(self, expand_meta=True):
        """Get sample definition string

        Return str
        """
        well_story = ""
        if self.is_meta() and (not expand_meta):
            well_story += self.get_meta_parts(as_name=True)
        else:
            # All wells, across all plates
            wspec = self.get_wspec()
            n_wells = self.num_wells(all_plates=True)
            well_story += f"{wspec} ({n_wells})"
            # Possibly list for multi-plates
            w_parts = plates.split_multi_plate96_wells(wspec)
            if len(w_parts) > 1:
                for p, wspec in enumerate(w_parts):
                    n_wells = self.num_wells(plate=p + 1)
                    well_story += f"   plate_{p+1} {wspec} ({n_wells})"
        return well_story

    def get_type_str(self):
        comb = self.get_combine_type_str()
        meta = self.get_meta_type_str()
        story = f"{comb}, {meta}"
        return story

    def get_combine_type_str(self):
        n_subs = self.num_sublibs()
        story = f"combined {n_subs} sublibs" if (n_subs > 1) else "single"
        return story

    def get_meta_type_str(self):
        n_samps = self.num_meta_parts()
        story = f"meta {n_samps} samples" if self.is_meta() else "simple"
        return story

    def get_wspec(self, sublib=0, plate=0, cannonical=False):
        """Get wells short-hand specification; e.g. 'A1:D2' or 'A1:D2__A1:B12'

        sublib = sublib index
        plate = 1-based plate index; 0 = all
        cannonical = flag to make spec cannonical form

        Return str
        """
        wspec = self._all_wspec_list[sublib]
        # Possibly break into multi-plate parts?
        if (plate > 0) and (plate < self.num_plates()):
            wspec = plates.split_multi_plate96_wells(wspec)[plate - 1]
        # Cannonical may be different than given
        if cannonical:
            wspec = plates.plate96_cannonical_wspec(wspec)
        return wspec

    def get_well_list(self, sublib=0, plate=1, sep="_"):
        """Get wells as list; ['A1', 'A2', 'B1'...] or ['A1__C1', 'A2__C1'...]

        sublib = sublib index
        plate = 1-based plate index if multiplate; If 0, get all combinations
        sep = seperator string for multi-plate combos (i.e. 'A1_C4')

        Return list
        """
        # Multi-plate combination case, create all plate 1 + 2 combos
        if self.is_multi_plate() and (plate == 0):
            r1_list = self.get_well_list(sublib=sublib, plate=1)
            r2_list = self.get_well_list(sublib=sublib, plate=2)
            well_list = [r1 + sep + r2 for r1, r2 in list(product(r1_list, r2_list))]
            return well_list

        # All plates or normal case, given plate
        if plate == 0:
            plate = 1
        well_list = self._all_well_list[sublib][plate - 1]
        well_list = sorted(well_list, key=natsort_keygen())
        return well_list

    def get_wind_list(self, sublib=0, plate=1, zero_pad=True, as_int=False, sep="_"):
        """Get well indexes; 1-based string ints; ['01', '02', '13'...'

        sublib = sublib index
        plate = 1-based plate index if multiplate; If 0, get all combinations
        zero_pad = flag to pad index strings (e.g. so same width with leading zero)
        as_int = flag to return list of int (not for combs)
        sep = seperator string for multi-plate combos (i.e. '01_09')

        Return list
        """
        # Multi-plate combination case, create all plate 1 + 2 combos
        if self.is_multi_plate() and (plate == 0):
            r1_list = self.get_wind_list(sublib=sublib, plate=1, zero_pad=zero_pad)
            r2_list = self.get_wind_list(sublib=sublib, plate=2, zero_pad=zero_pad)
            wind_list = [r1 + sep + r2 for r1, r2 in list(product(r1_list, r2_list))]
            return wind_list

        # Get max well index to set how wide may need to pad numbers
        max_wind = 0
        # Plate 0 gets list of all plates
        dim_list = self.get_plate_dims(plate=0)
        for n_row, n_col in dim_list:
            if (n_row * n_col) > max_wind:
                max_wind = n_row * n_col

        if plate == 0:
            plate = 1
        well_list = self.get_well_list(sublib=sublib, plate=plate)
        wind_list = []
        for w in plates.plate96_well_to_wind(well_list, as_int=True):
            if zero_pad:
                if max_wind > 99:
                    w = f"{w:03d}"
                else:
                    w = f"{w:02d}"
            else:
                w = str(w)
            wind_list.append(w)

        # Return as int?
        if as_int:
            wind_list = [int(w) for w in wind_list]
        return wind_list

    def num_wells(
        self,
        sublib=0,
        plate=1,
        all_sublibs=False,
        all_plates=False,
        as_frac=False,
        sublib_mean=True,
    ):
        """Count of wells

        sublib = sublib index
        plate = 1-based plate index; Plate 0 = all plates
        all_sublibs = Flag to sum across all sublibs
        all_plates = Flag to sum across all plates
        as_frac = Flag to return fraction of plate wells
        sublib_mean = Flag to return mean count for sublibs (combine only)

        Return int
        """
        if plate == 0:
            all_plates = True
        # If sublib mean, need to sum all sublibs first
        if sublib_mean:
            all_sublibs = True

        # Multiple or sum
        if as_frac:
            answer = 1
        else:
            answer = 0

        # Each possible sublibs
        for s in range(self.num_sublibs()):
            # Not this one?
            if (not all_sublibs) and (s != sublib):
                continue
            # Each possible plate
            for p in self.get_plate_indexes():
                if (not all_plates) and (p != plate):
                    continue
                n = len(self.get_well_list(sublib=s, plate=p))
                if as_frac:
                    n_row, n_col = self.get_plate_dims(plate=p)
                    answer *= n / (n_row * n_col)
                else:
                    answer += n

        # Mean of sublibs? As int if whole number
        if sublib_mean and (self.num_sublibs() > 1):
            answer = utils.whole_float_to_int(answer / self.num_sublibs())

        return answer

    def get_plate_dims(self, plate=1, as_num=False):
        """Get plate row, col number

        plate = 1-based plate index; If 0, return list

        Return tuple, list of tuples
        """
        if plate == 0:
            answer = [wp.get_dims() for wp in self._wplates]
        else:
            answer = self._wplates[plate - 1].get_dims()
        return answer

    def get_plate_num_wells(self, plate=1):
        """Get plate dim based number of wells

        plate = 1-based plate index; If 0, return list

        Return int
        """
        if plate == 0:
            dim_list = self.get_plate_dims(plate=plate)
        else:
            dim_list = [self.get_plate_dims(plate=plate)]

        # Each potential plate
        answer = 0
        for n_row, n_col in dim_list:
            answer += n_row * n_col
        return answer

    def get_sublib_sample(self, sublib=0):
        """Get (set of) well indexes"""
        return self._all_sub_samples[sublib]

    def is_allwell(self, meta=True):
        """True if all-well sample

        meta = flag to include meta sample

        return Boolean
        """
        if self.is_meta() and (not meta):
            return False

        # All well has full count of wells; plate 0 is all plates
        n_plate_wells = 0
        for n_row, n_col in self.get_plate_dims(plate=0):
            n_plate_wells += n_row * n_col
        n_samp_wells = self.num_wells(plate=0)

        # All well count = max
        return bool(n_samp_wells == n_plate_wells)

    def is_allsample(self, samp_list):
        """True if all-sample sample

        samp_list = list or dict with sample names or SplitPipe obj that can yield this

        all-sample is made of all non-meta samples except all-well

        return Boolean
        """
        # Get list of all simple, non-all-well sample names
        if isinstance(samp_list, list):
            pass
        elif isinstance(samp_list, dict):
            samp_list = samp_list.keys()
        elif isinstance(samp_list, spclass.SplitPipe):
            samp_list = samp_list.get_samples(as_name=True, meta=False, allwell=False)
        else:
            story = f"is_allsample given unknown 'kit' type {type(samp_list)}"
            raise Exception(story)

        # Lists of sample names for this sample and given list
        this_samp_set = set(self.get_meta_parts(as_name=True))
        simp_samp_set = set(samp_list)

        # All sample if collection of parts is the same as all simple (non-meta) samples
        return bool(this_samp_set == simp_samp_set)

    def is_combine(self):
        """If current sample combined mode or single?

        return Boolean
        """
        return bool(self.num_sublibs() > 1)

    def num_sublibs(self, min_wells=0, plate=0):
        """Get number of sublibs that have wells

        min_wells = qualifer to count sublib or not

        Even is no sublibs for combine mode, return 1 for current single sublib

        Return int >= 1
        """
        n_slib = len(self.get_sublibs())
        if not n_slib:
            n_slib = 1

        n = 0
        for s in range(n_slib):
            if len(self.get_well_list(sublib=s, plate=plate)) >= min_wells:
                n += 1
        return n

    def get_sublibs(self):
        """List of sublibrary structs

        Return list
        """
        return self._sublibs

    def get_plate_indexes(self):
        """Return list of plate indexes"""
        return [p + 1 for p in range(self.num_plates())]

    def num_plates(self):
        """Number of plates"""
        return len(self._wplates)

    def is_multi_plate(self):
        """If current sample has wells across multipe plates

        return Boolean
        """
        return self.num_plates() > 1

    def set_meta_part_samples(self, samp_list):
        """Set meta sample samples"""
        # Make sure list of samples
        for samp in samp_list:
            assert isinstance(samp, SpSamp), f"Expected SpSamp got {type(samp)}"
        self._meta_parts = samp_list
        self._am_meta = True

    def is_meta(self):
        """If current sample is a meta sample of others

        return Boolean
        """
        return self._am_meta

    def is_simple(self):
        """If current sample is simple (not-meta)

        return Boolean
        """
        return not self.is_meta()

    def num_meta_parts(self, meta=True, simple=True):
        """Get number of (meta) samples

        Return int
        """
        parts = self.get_meta_parts(meta=meta, simple=simple)
        return len(parts)

    def get_meta_parts(self, as_name=False, meta=True, simple=True, unwrap=False):
        """List of meta sample parts

        as_name     = flag to return names rather than objects
        meta        = flag to include meta samples
        simple      = flag to include simple (non-meta) samples
        unwrap      = flag to unwrap to only get non-meta (simple) samples

        Return list
        """
        # Special case to unwrwap down to non-meta sample parts
        if unwrap:
            return unwrap_meta_parts(self, as_name=as_name)

        answer = []
        for part in self._meta_parts:
            keep = True
            # Meta or simple qualifiers
            if self.is_meta():
                if not meta:
                    keep = False
            else:
                if not simple:
                    keep = False
            if keep:
                if as_name:
                    part = part.get_name()
                answer.append(part)
        return answer

    def set_status(self, status, step):
        """Set sample status flag and processing-step

        status = str [good|bad|issue] or true/false
        step = str
        """
        if not isinstance(status, str):
            status = "good" if status else "bad"
        else:
            if status.lower()[0] == "g":
                status = "good"
            elif status.lower()[0] == "b":
                status = "bad"
            else:
                status = "issue"
        self._status = status
        self._stat_step = step

    def get_status(
        self, status=False, step=False, is_good=False, is_bad=False, is_issue=False
    ):
        """Return status flag and step

        status = flag to return status
        step = flag to return step

        Return string, bool or tuple(str, str)
        """
        # All is_* flags are about status
        if is_good or is_bad or is_issue:
            status = True
        # Status, step or tuple of both
        if status:
            if is_good:
                answer = True if self._status == "good" else False
            elif is_bad:
                answer = True if self._status == "bad" else False
            elif is_issue:
                answer = True if self._status == "issue" else False
            else:
                answer = self._status
        elif step:
            answer = self._stat_step
        else:
            answer = self._status, self._stat_step

        return answer


# ---------------------- Non-class functions ------------------------------


def compare_samples(f_samp, s_samp, name=True, wells=True, parts=True):
    """Compare two samples

    name    = Flag to compare name
    wells   = Flag to compare wells (and num plates)
    parts   = Flag to compare meta-sample parts

    Return tuple (status, story if false)
    """
    same = True
    story = []
    # Compare indicated parameters
    if name:
        if f_samp.get_name() != s_samp.get_name():
            same = False
            story.append(f"Names differ: {f_samp.get_name()} vs {s_samp.get_name()}")
    if wells:
        f_nplates = f_samp.num_plates()
        s_nplates = s_samp.num_plates()
        if f_nplates != s_nplates:
            same = False
            story.append(f"Well plate numbers differ: {f_nplates} vs {s_nplates}")
        else:
            for p in range(f_nplates):
                f_list = f_samp.get_well_list(plate=p)
                s_list = s_samp.get_well_list(plate=p)
                if f_list != s_list:
                    same = False
                    story.append(
                        f"Well lists differ: lens {len(f_list)} vs {len(s_list)}"
                    )
    if parts:
        f_parts = f_samp.get_meta_parts(as_name=True)
        s_parts = s_samp.get_meta_parts(as_name=True)
        if f_parts != s_parts:
            same = False
            story.append(f"Part lists differ: {f_parts} vs {s_parts}")
    return same, story


def sample_well_overlap(f_samp, s_samp):
    """Check for overlap of sample wells, including multi-plate case

    Return number, story
    """
    n = 0
    story = ""

    assert (
        f_samp.num_plates() == s_samp.num_plates()
    ), f"Differing number of plates: {f_samp.num_plates()} {s_samp.num_plates()}"

    # plate 0 = all plate combos in multi-plate case
    f_wells = set(f_samp.get_well_list(plate=0))
    s_wells = set(s_samp.get_well_list(plate=0))
    # Shared in common
    com_wells = f_wells & s_wells
    if com_wells:
        n = len(com_wells)
        what = "well combintations" if f_samp.is_multi_plate() else "wells"
        story = (
            f"Samples {f_samp.get_name()} and {s_samp.get_name()} overlap {n} {what}: "
        )
        story += ",".join(com_wells)

    return n, story


def unwrap_meta_parts(samp, as_name=False):
    """Unwrap sample to basis set of non-meta (simple) parts

    Return list
    """
    # Recursively collect list of samples
    answer = []
    if samp.is_meta():
        for part in samp.get_meta_parts():
            if part.is_meta():
                answer.extend(unwrap_meta_parts(part))
            else:
                answer.append(part)
    else:
        answer.append(samp)
    # Make sure samples are unique; Maybe return just name
    answer = sorted(list(set(answer)), key=lambda s: s.get_name())
    if as_name:
        answer = [s.get_name() for s in answer]
    return answer


def compatible_meta_samples(samples, fatal=True):
    """Check list of samples for compatability as meta sample parts

    return tuple (status, story)
    """
    ok = True
    story = ""

    # Make sure all samples are legit
    if ok:
        for i, samp in enumerate(samples):
            if type(samp) != SpSamp:
                ok = False
                story = f"Expected SpSamp got {type(samp)} for item {i+1}"
                break
    # Check if any sublibs all agree
    if ok:
        f_list = samples[0].get_sublibs()
        for i, samp in enumerate(samples[1:]):
            s_list = samp.get_sublibs()
            if s_list != f_list:
                ok = False
                story = "(Meta)sample sublib lists differ\n"
                story += f"1 {samples[0].get_name()} list: {f_list}"
                story += f"{i+1} {samp.get_name()} list: {s_list}"
                break
    if fatal and not ok:
        raise Exception(story)

    return ok, story


def sort_samp_list_sim_met_all(samples, sort_name=False):
    """Sort sample list so simple, then meta, then all(well or sample)

    samples = dict or list of sample objs
    sort_name = flag to sort on sample name

    Return list
    """
    # Dict with names as keys
    if isinstance(samples, dict):
        s_dict = samples
    else:
        s_dict = {s.get_name(): s for s in samples}

    # Names, maybe sort order
    name_list = s_dict.keys()
    if sort_name:
        name_list = sorted(name_list)
    # Lists of simple, meta, and all-x samples
    s_list = []
    m_list = []
    a_list = []
    for name in name_list:
        samp = s_dict[name]
        if name in ["all-well", "all-sample"]:
            a_list.append(name)
        elif samp.is_meta():
            m_list.append(name)
        else:
            s_list.append(name)

    # From names back to obj list
    answer = []
    name_list = s_list + m_list + a_list
    for name in name_list:
        answer.append(s_dict[name])

    return answer


def order_sample_list(slist, sort_name=False):
    """Order sample list to allow meta sample processing

    Meta samples have to reference previously defined samples

    slist       = sample list; May be objs or list of sample def lists
    sort_name   = flag to sort by name

    Return tuple (status, sorted list, story)
    """
    # List of names and dict with either defs or sample objs
    name_list = []
    in_dict = {}
    as_def = False
    for s in slist:
        if isinstance(s, SpSamp):
            name = s.get_name()
            val = s
        elif isinstance(s, list):
            # As sample definition string <name> <def>
            name, val = s[:2]
            as_def = True
        else:
            story = f"order_sample_list: unknown list item: {type(s)}"
            raise Exception(story)
        in_dict[name] = val
        name_list.append(name)
    if sort_name:
        name_list = sorted(name_list)

    # List of outputs and set to keep track of what's been assigned
    simple_list = []
    order_list = []
    order_set = set()
    ok = True
    story = ""

    # First pass = collect non meta samples (simple well defintions)
    for name in name_list:
        val = in_dict[name]
        if as_def:
            # Well definition; Returns tuple (status, well_list, well_list_list, story-if-problem)
            # If parse status is good, this is not a meta sample
            not_meta, _, _, _ = plates.parse_multi_plate96_wells(val, as_well=False)
        else:
            # val is a sample obj
            not_meta = val.is_simple()
        if not_meta:
            simple_list.append(name)
            if as_def:
                order_list.append([name, val])
            else:
                order_list.append(val)
            order_set.add(name)

    # Second pass for meta samples
    while ok and len(order_list) < len(in_dict):
        n_add = 0
        for name in name_list:
            # Only consider not-yet-ordered samples; If already have, skip
            if name in order_set:
                continue
            val = in_dict[name]

            # Meta definition as comma,sep,samples or wildcard*
            if as_def:
                def_str = val
            else:
                def_str = val.get_meta_parts(as_name=True)
            part_list = get_meta_wspec_parts(def_str, order_set)
            if len(part_list) < 2:
                ok = False
                story = f"Problem with sample def: '{name}' '{def_str}''"
                break

            # If all meta sample parts are already found, save current
            if set(part_list).issubset(order_set):
                n_add += 1
                order_set.add(name)
                if as_def:
                    order_list.append([name, def_str])
                else:
                    order_list.append(val)

        # If no additions, we're stuck (cyclic graph, etc?)
        if not n_add:
            ok = False
            story = "Problem parsing wells / resolving meta-parts for sample def(s);"
            for name in name_list:
                if name not in order_set:
                    val = in_dict[name]
                    story += f"\nSample: '{name}' '{val}'"

    return ok, order_list, story


def get_meta_wspec_parts(wspec, name_list):
    """Get list of names for meta well spec; Comma list or wildcard

    wspec = well specification (e.g. name list, wildcard for meta sample)
    name_list = source of names to match (e.g. non-meta base sample names)

    Return list of found/matching names
    """
    if isinstance(name_list, dict):
        name_list = name_list.keys()

    part_list = []
    # Comma sep list
    if "," in wspec:
        for part in wspec.split(","):
            # Wildcard part?
            if "*" in part:
                new_parts = get_meta_wspec_parts(part, name_list)
                if new_parts:
                    # Only new matches
                    for npart in new_parts:
                        if npart not in part_list:
                            part_list.append(npart)
            else:
                # Simple part has to be in list, else return empty list
                if (part in name_list) and (part not in part_list):
                    part_list.append(part)
    # Wildcard
    elif "*" in wspec:
        regex_pattern = wspec.replace("*", ".*")
        regex = re.compile(regex_pattern)
        part_list = [name for name in name_list if regex.match(name)]

    return part_list


def clean_sample_def_name(name, ok_num_start=True):
    """Clean given sample name chars and so not ambiguous with respect to wells

    Clean up restricted chars, etc (names are part of filepaths in output)
    If a sample is named as a well, change; e.g. 'd1' >--> 'sample_d1'

    Return name, story
    """
    well_list = plates.wells_for_96plate()

    new_name = name
    # Clean chars; "non-variable" chars go to underscore.
    # Allow dash, and name can start with number
    new_name = utils.clean_vname_chars(
        name, to_under=True, dash_ok=True, int_first=True
    )
    if not new_name:
        new_name = "sample"
    # Name like a well?
    if new_name.upper() in well_list:
        new_name = "sample_" + new_name
    # Naked number? (Note clean step changes floats like '1.3' >--> '1_3'
    ok, val = utils.str_to_num(new_name)
    if ok:
        new_name = "sample_" + new_name
    # Don't start with a number? May need prefix
    if not ok_num_start:
        ok, val = utils.str_to_num(new_name[0])
        if ok:
            new_name = "sample_" + new_name
    # Name change?
    story = ""
    if new_name != name:
        story = f"Cleaned sample name '{name}' >--> '{new_name}'"

    return new_name, story


def is_multi_plate_wspec(wspec, strict=False):
    """Check if well specifcation is multi-plate

    wspec = well specification (e.g. "A1:B6", "A1:B6__C1:D12") or list of these

    Return bool
    """
    # If list, bail on first multi (assume all are the same)
    if isinstance(wspec, list):
        for ws in wspec:
            if is_multi_plate_wspec(ws):
                return True
        return False

    # One spec; Parse returns list of separate well specs
    ok, plate_wells, _, story = plates.parse_multi_plate96_wells(wspec)
    if not ok:
        # If not strict (e.g. meta-samples), then don't fail
        if strict:
            story += f"\nWell spec: |{wspec}|"
            raise ValueError(story)

    # More than one plate?
    multi = True if len(plate_wells) > 1 else False
    return multi


def get_wspec_for_all_well(n_rows, n_plates=1):
    """Get well short hand spec for all-well sample given plate dims

    Return string (e.g. 'A1:D12')
    """
    assert isinstance(n_rows, int), f"Expected int got {type(n_rows)}"

    # Row and col lists (of string)
    rows, cols = plates.rows_cols_for_96plate(n_rows=n_rows)
    # First row, col to last row, col
    wspec = f"{rows[0]}{cols[0]}:{rows[-1]}{cols[-1]}"

    # Multiplate? Second plate has 8 rows
    if n_plates > 1:
        rows, cols = plates.rows_cols_for_96plate(n_rows=8)
        wspec = wspec + "__" + f"{rows[0]}{cols[0]}:{rows[-1]}{cols[-1]}"

    return wspec


def get_sample_specs(samp):
    """Get define settings for sample (e.g. for json)

    Return dict
    """
    ndict = {}
    ndict["name"] = samp.get_name()
    ndict["wspec"] = samp.get_wspec()

    # Saved per plate
    ndict["num_plates"] = samp.num_plates()
    plate_dim_list = []
    well_num_list = []
    well_list_list = []
    for p in samp.get_plate_indexes():
        plate_dim_list.append(list(samp.get_plate_dims(plate=p)))
        well_num_list.append(samp.num_wells(plate=p))
        well_list_list.append(samp.get_well_list(plate=p))
    ndict["plate_dims"] = plate_dim_list
    ndict["num_wells"] = well_num_list
    ndict["well_list"] = well_list_list

    # Status
    status, stat_step = samp.get_status()
    ndict["status"] = status
    ndict["stat_step"] = stat_step

    # Parts for meta sample
    ndict["meta"] = samp.is_meta()
    if samp.is_meta():
        samples = samp.get_meta_parts()
        samp_list = [s.get_name() for s in samples]
        ndict["sample_num"] = len(samp_list)
        ndict["sample_list"] = samp_list

    # Sublibs for combine mode
    ndict["combined"] = samp.is_combine()
    if samp.is_combine():
        sublibs = samp.get_sublibs()
        slib_list = []
        for i, slib in enumerate(sublibs):
            slib_dict = {
                "sublib_name": slib.get_name(),
                "sublib_wells": samp.get_wspec(sublib=i),
            }
            slib_list.append(slib_dict)
        ndict["sublib_num"] = len(slib_list)
        ndict["sublib_list"] = slib_list
    return ndict


def write_def(samp, ofname):
    """Write def / specs to file

    Returns nothing
    """
    o_dict = {}
    o_dict["header"] = utils.get_header_dict(desc="Sample specification")
    o_dict["sample"] = get_sample_specs(samp)
    utils.write_json(o_dict, ofname)
