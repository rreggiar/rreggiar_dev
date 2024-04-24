#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import os
import sys
import pandas as pd


LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
    import spsamp
else:
    from splitpipe import utils
    from splitpipe import spsamp


# ---------------------------------------------------------------------------
class SpOut:
    """Class for pipeline output dirs"""

    def __init__(self, fpath, name=""):
        """Initialize structure; Minimally need path"""
        self.fpath = ""
        self.name = ""
        self.kit = ""
        self.use_case = ""
        self.genome_dir = ""
        self._sample_dict = {}

        # Get directory file dict
        f_dict = spoutdir_paths(fpath)
        if not f_dict:
            story = f"Failed to parse top level dir: Path='{fpath}'"
            raise ValueError(story)

        # Set path to parsed top (may not be exactly given fpath)
        self.fpath = f_dict["DIR_TOP"]
        # Name is last part of top path unless given
        if not name:
            name = self.fpath.split("/")[-1]
        self.name = name

        # Init run process info with object values
        odict = {"object": {"path": fpath, "name": name}}
        # Get run process info dict defining everything
        pdict = spoutdir_runproc(f_dict["DIR_TOP"])
        self._runproc_info = {**odict, **pdict}

        # Instantiate defined samples if not a combine mode sublib
        if not pdict["run"]["mode"].startswith("com"):
            # Save meta sample definitions for second pass with simple samples
            meta_def_list = []
            # Samples as list of dicts under 'sample' top level dict
            for sdict in pdict["sample"]["samples"]:
                # Meta or simple? Create simple sample, save meta as tuple here
                if sdict["meta"]:
                    meta_def_list.append(
                        [
                            sdict["name"],
                            sdict["sample_list"],
                            sdict["plate_dims"],
                        ]
                    )
                else:
                    plate_dims = sdict["plate_dims"]
                    samp = spsamp.SpSamp(
                        sdict["name"], wspec=sdict["wspec"], plate_dims=plate_dims
                    )
                    samp.set_status(sdict["status"], sdict["stat_step"])
                    self._sample_dict[samp.get_name()] = samp

            if not self._sample_dict:
                story = f"Failed to define any samples: Path='{fpath}'"
                raise ValueError(story)

            # Now create any meta samples
            for m_name, m_name_list, plate_dims in meta_def_list:
                sample_list = []
                for sname in m_name_list:
                    sample_list.append(self._sample_dict[sname])
                samp = spsamp.SpSamp(m_name, samples=sample_list, plate_dims=plate_dims)
                samp.set_status(sdict["status"], sdict["stat_step"])
                self._sample_dict[samp.get_name()] = samp

        # Set for comparision of sublibs
        self.chemistry = pdict["chemistry"]["name"]
        self.kit = pdict["kit"]["kit_name"]
        self.use_case = pdict["run"]["use_case"]
        self.genome_dir = pdict.get("genome_dir", "")

    def __del__(self):
        pass

    def __repr__(self):
        ostring = "split-pipe output dir (SpOut obj)\n"
        ostring += f"Name:     {self.name}\n"
        ostring += f"Path:     {self.fpath}\n"
        ostring += f"Kit:      {self.kit}\n"
        ostring += f"Use_case: {self.use_case}\n"
        ostring += f"Samples count {len(self._sample_dict)}\n"
        ostring += f"Genome:   {self.genome_dir}"
        return ostring

    def get_name(self):
        return self.name

    def get_path(self):
        return self.fpath

    def get_runproc_info(self, key=None):
        """Get runproc dict or value of one key

        key = possible info key to extract; Can be multi.level via dots

        Return tuple (dict or key value, source)
        """
        pdict = self._runproc_info
        source = self.get_path()
        if key:
            # This handles nested.keys (e.g. 'run.use_case')
            pdict = utils.get_nested_dict_val(pdict, key)
        return (pdict, source)

    def get_samples(self, one="", as_name=False):
        """List or named sample structs (not input strings)"""
        answer = None
        s_dict = self._sample_dict
        if one:
            answer = s_dict.get(one)
        else:
            if as_name:
                answer = list(s_dict.keys())
            else:
                answer = list(s_dict.values())
        return answer

    def proc_sample_defs(self, f_list):
        """Create SpSamp objects from definition files in list"""
        # Init empty dict
        self._sample_dict = {}
        # Collection of samples from definitions
        for path in f_list:
            def_dict = utils.parse_def_dict(path, keys="name,wells,type")
            if not def_dict:
                story = f"PROBLEM: Failed to parse: '{path}'"
                print(story)
                self._sample_dict = {}
                break
            # Must be normal (not combined) type
            if def_dict["type"].lower() != "normal":
                story = f"PROBLEM: Sample {path} type not normal: {def_dict}"
                print(story)
                self._sample_dict = {}
                break
            # Create and save
            samp = spsamp.SpSamp(
                def_dict["name"],
                wspec=def_dict["wspec"],
                plate_dims=(def_dict["plate_n_rows"], def_dict["plate_n_cols"]),
            )
            self._sample_dict[samp.get_name()] = samp

    def report_outdir(self, ofile=None):
        """Report subdirectory details"""
        if not ofile:
            ofile = sys.stdout
        print("SpOut", file=ofile)
        print(f"Name: {self.name}", file=ofile)
        print(f"Path: {self.fpath}", file=ofile)
        for samp in self.get_samples():
            print(f"Sample: {samp.get_name()} {samp.get_wspec()}", file=ofile)


# -------------------------------------------------------------------------
#              Non-class functions
# -------------------------------------------------------------------------


def spout_equal_key_vals(odir1, odir2, rpkey):
    """Compare SpOuts on given runproc_info key

    Return tuple (bool, story)
    """
    assert isinstance(odir1, SpOut) and isinstance(odir2, SpOut)

    val1, _ = odir1.get_runproc_info(key=rpkey)
    val2, _ = odir2.get_runproc_info(key=rpkey)
    equal = False
    story = ""
    if val1 and val2:
        if val1 == val2:
            equal = True
        else:
            story = f"Different {rpkey}; '{val1}' vs val2 '{val2}'"
    else:
        story = f"Missing {rpkey}; val1 '{val1}' and/or val2 '{val2}'"
    return equal, story


def spoutdir_obj(top_dir, name=""):
    """Get output dir object (if top_dir is legit)

    Return SpOut
    """
    obj = None
    fpath = spoutdir_str(top_dir)
    if fpath:
        obj = SpOut(fpath, name=name)
    return obj


def spoutdir_str(top_dir):
    """Get output dir path str (if top_dir is legit)

    Return str
    """
    fpath = ""
    if isinstance(top_dir, str):
        if is_spoutdir(top_dir):
            fpath = top_dir
    elif isinstance(top_dir, SpOut):
        fpath = top_dir.get_path()
    else:
        raise ValueError(f"Unsupported top_dir type {type(top_dir)}")
    return fpath


def spoutdir_runproc(top_dir, key=None):
    """Get run process info for pipeline output dir

    top_dir = path or obj
    key = subfield for proc info; May be multi.level (e.g. 'run.use_case')

    Return dict (nested json)
    """
    # print(f">>> spout spoutdir_runproc")
    ndict = {}
    # Different potential sources
    if isinstance(top_dir, str):
        # print(f">>> spoutdir_runproc top_dir str |{top_dir}|")
        fdict = spoutdir_paths(top_dir)
        # print(f"+ fdict {fdict}")
        if fdict and fdict["PF_RUNPROC_DEF"]:
            ndict = utils.read_json(fdict["PF_RUNPROC_DEF"])
    elif isinstance(top_dir, SpOut):
        # print(">>> spoutdir_runproc top_dir SpOut")
        fdict, _ = top_dir.get_runproc_info()
    else:
        raise ValueError(f"Unsupported top_dir type {type(top_dir)}")
    # If have dict and key, get sub field
    if ndict and key:
        ndict = utils.get_nested_dict_val(ndict, key)
    # print(f"<<< spoutdir_runproc str {type(ndict)}")
    return ndict


def spoutdir_paths(fpath, dir_slash=False):
    """Check if given path is a (non-mkref) pipeline output dir

    Return dict with various file paths that are found
    """
    # Get cannonical starting top
    # print(f"# >>> spoutdir_paths |{fpath}|")
    top_dir = find_spoutdir(fpath, dir_slash=dir_slash)
    if not top_dir:
        # print(f"<<< spoutdir_paths; No top_dir")
        return {}

    o_dict = None
    run_proc_def = ""
    sample_dirs = {}

    # Top level will always have process
    proc_subdir = utils.check_infile(
        "process", path=top_dir, verb=False, dir_slash=dir_slash
    )
    # print(f"+ spoutdir_paths; top={top_dir} proc={proc_subdir}")
    # Run process def file
    fnames = utils.fnames_for_path(path=proc_subdir, depth=1, name="run_proc_def.json")
    if len(fnames) == 1:
        run_proc_def = fnames[0]

    # Sample sub_dirs all contain subdir named "report"; depth 2 <top_dir>/<samp>/report
    fnames = utils.fnames_for_path(
        path=top_dir, depth=2, get_subdirs=True, get_files=False, name="report"
    )
    rep_end = "/report"
    for report in fnames:
        if report.endswith(rep_end):
            samp_dir = report[: -len(rep_end)]
            samp_name = samp_dir.split("/")[-1]
            if dir_slash:
                samp_dir += "/"
            sample_dirs[samp_name] = samp_dir

    # Return dict with file paths
    o_dict = {
        "DIR_TOP": top_dir,
        "DIR_PROC": proc_subdir,
        "PF_RUNPROC_DEF": run_proc_def,
        "DIR_SAMP": sample_dirs,
    }
    # print(f"<<< spoutdir_paths", o_dict)
    return o_dict


def find_spoutdir(path, dir_slash=False):
    """Attempt to find top-level pipeline output dir from path

    Return path if found
    """
    top_dir = ""
    # Check paths starting from given
    parts = path.split("/")
    while parts:
        tpath = "/".join(parts)
        if is_spoutdir(tpath):
            top_dir = tpath
            break
        parts = parts[:-1]
    # Add trailing slash if found and flag
    if top_dir and dir_slash:
        top_dir += "/"
    return top_dir


def is_spoutdir(path):
    """Guess if path is pipeline output top_dir

    Return Boolean
    """
    # print(f">>> is_spoutdir |{path}|")
    is_top = False
    if os.path.isdir(path):
        _dirs = _files = []
        # Only want first iteration (dir sub search) then break
        #   Otherwise: for root, dirs, files in os.walk(path):
        for _, _dirs, _files in os.walk(path):
            break
        # both dirs and files
        is_top = check_spoutdir_subdirs(_dirs) and check_spoutdir_files(_files)
    # print(f">>> is_spoutdir {is_top}")
    return is_top


def check_spoutdir_subdirs(dirs, min_dirs=1):
    """Look for subdirs indicating pipeline output top_dir

    Return Boolean
    """
    found = 0
    for subdir in dirs:
        if subdir == "process":
            found += 1
    return found >= min_dirs


def check_spoutdir_files(flist, min_files=1):
    """Look for files indicating pipeline output top_dir

    Return Boolean
    """
    # Specific files
    named_list = "agg_samp_ana_summary.csv,all_summaries.zip".split(",")
    found = 0
    for fname in flist:
        if fname in named_list:
            found += 1
            continue
        if fname.startswith("split-pipe_"):
            if fname.endswith(".log") or fname.endswith(".log.gz"):
                found += 1
    return found >= min_files


def clear_env_files(top_dir, verb=True):
    """Remove run environment files

    Return number of files removed
    """
    f_lis = []
    f_lis.append(f"{top_dir}/process/run_proc_def.json")
    if verb:
        print(f"Checking to remove {len(f_lis)} files:")
    # Go through list explicitly to get status for each one
    tot = 0
    for fpath in f_lis:
        n = utils.check_and_rm_file(fpath)
        if n < 0:
            print(f"Problem removing file: {fpath}")
            break
        tot += n
    return tot


def load_top_dir_cell_bcs(path):
    """Load all cell barcode info from top dir

    Return tuple (dict, story); dict[sample] = set of barcodes
    """
    # Make sure top dir
    path = find_spoutdir(path)
    if not path:
        return (None, f"Not spipe top dir: {path}")

    # Dict and set to tally all unique barcodes (for reporting)
    new_dict = {}
    all_bc_set = set()
    # Collect barcodes from all (filtered) cell data files
    flist = utils.fnames_for_path(
        path=path, get_subdirs=False, get_files=True, name="cell_metadata.csv"
    )
    for cpath in flist:
        # Parts of path like this;
        #   /analysis/sublib2/samp4/DGE_targeted_filtered/cell_metadata.csv
        parts = cpath.split("/")
        dge_dir = parts[-2]
        samp_dir = parts[-3]

        # Don't want unfiltered
        if "unfiltered" in dge_dir:
            continue

        # Only barcode columns
        df = pd.read_csv(cpath, usecols=["bc_wells"])
        new_dict[samp_dir] = set(df["bc_wells"])
        all_bc_set |= new_dict[samp_dir]

    story = f"Loaded {len(all_bc_set)} cell barcodes in {len(new_dict)} samples"
    # Add final union set
    new_dict["all_bc_set"] = all_bc_set

    return (new_dict, story)


def get_sublib_name_dict(sub_dirs):
    """Get unique names for list of sublibraries

    Assumes listed paths already cleaned up (e.g. unique, expanded, no trailing '/')
    Names consider path parts, end toward root until unique, like:

        /x/y/z/sub1    >-->    z__sub1
        /x/y/z/sub2    >-->    sub2
        /x/y/w/sub1    >-->    y__w__sub1
        /x/v/w/sub1    >-->    v__w__sub1

    Return dict with [path] = name
    """
    # Init names to nothing
    sub_dict = {s: "" for s in sub_dirs}
    p = 1
    while True:
        # Unique and duplicated subsets of current dict
        u_dict, d_dict = utils.sep_dic_by_dup_vals(sub_dict)
        # No duplicate names = done
        if not d_dict:
            sub_dict = u_dict
            break
        # Update names in dup dict, including next path part; Add to unqiue dict
        for k in d_dict.keys():
            parts = k.split("/")
            u_dict[k] = "__".join([p for p in parts[-p:] if p])
        p += 1
        sub_dict = u_dict
    return sub_dict


def compare_spoutdirs(
    f_out,
    s_out,
    ver=True,
    use_case=True,
    chem=True,
    kit=True,
    genome=True,
    samples=False,
):
    """Compare first and second output dirs with variuos critera

    ver     = Flag to compare version
    use_case = Flag to compare use_case
    kit     = Flag to compare kit
    genome  = Flag to compare genome
    samples = Flag to compare samples

    Return tuple (status, story if false)
    """
    same = True
    story = []
    # Compare indicated keys
    if same and ver:
        ok, eqstr = spout_equal_key_vals(f_out, s_out, "header.version")
        if not ok:
            same = False
            story.append(eqstr)
    if same and use_case:
        ok, eqstr = spout_equal_key_vals(f_out, s_out, "run.use_case")
        if not ok:
            same = False
            story.append(eqstr)
    # Kit; Compare name and number
    if same and kit:
        comp_keys = "kit_name,kit_num".split(",")
        for ckey in comp_keys:
            key = f"kit.{ckey}"
            ok, eqstr = spout_equal_key_vals(f_out, s_out, key)
            if not ok:
                same = False
                story.append(eqstr)
    # Chemistry
    if same and chem:
        ok, eqstr = spout_equal_key_vals(f_out, s_out, "chemistry.name")
        if not ok:
            same = False
            story.append(eqstr)
    # Genome has various things that can be compared ...
    if same and genome:
        comp_keys = "genome_num,genome_names,genome_total_size".split(",")
        for ckey in comp_keys:
            key = f"genome.{ckey}"
            ok, eqstr = spout_equal_key_vals(f_out, s_out, key)
            if not ok:
                same = False
                story.append(eqstr)

    # Samples; First collection of names, then sample details
    if same and samples:
        f_samps = set(f_out.get_samples(as_name=True))
        s_samps = set(s_out.get_samples(as_name=True))
        if not (f_samps == s_samps):
            same = False
            dif_str = f"Sample lists differ: {sorted(f_samps)} -vs- {sorted(s_samps)}"
            story.append(dif_str)
    # Now details; Names are the same at this stage
    if same and samples:
        for samp_name in sorted(f_samps):
            f_samp = f_out.get_samples(one=samp_name)
            s_samp = s_out.get_samples(one=samp_name)
            if not spsamp.compare_samples(f_samp, s_samp):
                same = False
                dif_str = f"Sample {samp_name} details differ: {f_samp} -vs- {s_samp}"
                story.append(dif_str)

    # Return (status, story)
    return same, "\n".join(story)
