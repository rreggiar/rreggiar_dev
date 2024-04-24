#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

import os
import re
import psutil
import inspect
import random
import numpy as np
import pandas as pd
from collections import OrderedDict
from itertools import product
from natsort import natsort_keygen

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import bcutils
    import ccutils
    import crispr
    import genes
    import kits
    import maptarg
    import pb_const as const
    import pb_mem
    import plates
    import run_modes as modes
    import sltab
    import spsamp
    import spoutdir as spout
    import tcr
    import utils
else:
    from splitpipe import bcutils
    from splitpipe import ccutils
    from splitpipe import crispr
    from splitpipe import genes
    from splitpipe import kits
    from splitpipe import maptarg
    from splitpipe import pb_const as const
    from splitpipe import pb_mem
    from splitpipe import plates
    from splitpipe import run_modes as modes
    from splitpipe import sltab
    from splitpipe import spsamp
    from splitpipe import spoutdir as spout
    from splitpipe import tcr
    from splitpipe import utils


# ---------------------------------------------------------------------------
class SplitPipe:
    """Main singleton class for pipeline"""

    def __init__(self, version, args=None):
        """Initialize; Set all parameters to init values

        version = program version; i.e. to report out
        args = command line arguments

        Parameter precedence:
            Defaults < parameter file < command line

        Any issues / problems are set in the structure itself
        """
        print(f"# Init SplitPipe {version}")
        self._version = version

        # Defaults
        self._init_defaults()

        # Standardized args keys and get into dict
        args = const.clean_par_dict_keys(args)

        # Parameter file?
        if "parfile" in args:
            self._parse_params(args["parfile"])
        # Any parse problems so far?
        if self.any_problem(report=True):
            return

        # Set command line args into parameters
        self._save_command_line()
        if args:
            self._add_comargs(args)
            self._reset_temp_comarg_params()
        self._pars = const.clean_par_dict_keys(self._pars)

        # File paths
        self._pkg_path = os.path.dirname(__file__)
        self._init_filepath_dict()

        # Init / pre-setup functions for specific things; Order may matter
        # Any issues are reported as problems which are caught in eval loop
        call_lis = [
            # Expand any filename parameters
            "_init_fnames",
            # If sublib list given, load and expand these paths
            "_init_sublib_list",
            "_init_combine",
            # Runtime process info from PF_RUNPROC_DEF files; From sublib if combine
            "_init_runproc_info",
            "_init_mode",
            # Set working info / any missing parameters; These don't fail
            "_init_logfname",
            "_init_sys_info",
            "_init_nthreads",
            "_init_rand_seed",
            "_init_in_fastqs",
            # These can fail (i.e. set problem)
            "_init_parent",
            # Chemistry, kit, samples, genome, etc may come from parent
            "_init_chemistry",
            "_init_kit",
            "_init_barcodes",
            "_init_sample_list",  # After know kit
            "_init_genome",
            "_init_tcr",
            "_init_targen",
            "_init_fb_crispr",
            "_init_gene_urls",
            "_init_dge_vars",
            "_init_map_targs",
            "_init_mkref",
            "_init_run_name",
            # Set use case and run steps after everything else
            "_init_use_case",
            "_init_run_steps",
        ]
        # Call each function and check for problems
        for call in call_lis:
            # print(call)
            eval(f"self.{call}")()
            if self.any_problem(report=True):
                break

    def __del__(self):
        self.close_down()

    def __repr__(self):
        use_case = self.get_use_case()
        mode = self.get_mode()
        # Current mode only if anything is set
        cur_mode = self.get_cur_mode()
        if cur_mode:
            cur_mode = f"\t({cur_mode})"
        chem = self.get_chemistry()
        kit_s, kit_n, _ = self.get_kit()
        ostring = "SplitPipe\n"
        ostring += f"Use_case:  {use_case}\n"
        ostring += f"Mode:      {mode}{cur_mode}\n"
        if not self.get_use_case("mkref"):
            ostring += f"Chemistry: {chem}\n"
            ostring += f"Kit:       {kit_s}\t({kit_n} wells)"
        return ostring

    def _init_defaults(self):
        """Set all fields to initial / default values"""
        print("# Initializing defaults")
        # State values; Collection of init lines, Problems, calls
        self._init_story = []
        self._command_line = []
        self._prob_story = ""
        self._prob_count = 0
        self._time0 = self._s_time = utils.get_time(as_str=False)
        self._last_calls = []
        # runtime process info dict and source path
        self._runproc_info = {}
        self._runproc_source = ""
        self._top_dir = None
        # working mode and list of run steps
        self._cur_mode = ""
        self._run_steps = []

        # Par and log file
        self._parfname = ""
        self._logfname = ""
        self.log_FH = None

        # Initialize empty sample, map_target, sublib, object collections ...
        self._sample_dict = OrderedDict()
        self._mtarg_dict = OrderedDict()
        self._sublib_dict = OrderedDict()
        # Excutables
        self._execs = {}
        # Barcode, genome, target enrichment, parent dir info structs
        self._bc_info = None
        self._gene_info = None
        self._guide_info = None
        self._targen_info = None
        self._parent_info = None

        # Working data objs
        self._tscp_assign = None
        self._stat_df = None
        self._stat_df_sample = self._stat_df_top_dir = ""
        self._fq0_samp = None
        self._fq1_samp = None
        self._fq2_samp = None
        self._vbc_count = None
        self._2vbc_count = None
        # Kit score details
        self._kit_score_info = None
        # Sample plate; i.e. based on kit, multi-plate
        self._wplates = None
        # Cell barcodes (dict of sets per sample; e.g. from parent dir)
        self._cell_bc = None
        # gene url format string collection for reporting
        self._gene_url_dict = {}

        # Default parameters defined in const module
        self._pars = const.get_spipe_default_pars()
        # Add temp command-line-only params for parsing
        temp_pars = const.get_spipe_temp_par_dict()
        for k in temp_pars.keys():
            self._pars[k] = None

    def _parse_params(self, parfname):
        """Parse parameter file (if any)

        All raw <key> <value> [<val>] lists kept in par_pfraw
        Values also copied into (i.e. set in) main par collection

        Sets problem if encountered

        Return status
        """
        self._parfname = parfname
        self._par_pfraw = {}
        if parfname:
            print(f"# Parsing parameters from {parfname}")
            line_n = 0
            try:
                # Lists for sample, map_targ, gene url, barcode definitions
                sample_list = []
                mtarg_list = []
                gurl_list = []
                bc_rset_list = []
                gfasta_list = []
                fasta_slice = []
                # First pass read in, check and save
                INFILE = open(parfname)
                for line in INFILE:
                    line_n += 1
                    # Strip off anything after comment, then split to <key> <vals...>
                    line = line.split("#")[0].strip()
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    # Keys are lowercase and any period replace with underscore
                    key = const.clean_param_key(parts[0])
                    if key not in self._pars:
                        story = f"Unrecoginzed parameter key: '{key}'\n"
                        story += f"Problem line {line_n}: '{line}'"
                        self.set_problem(story)
                        return False
                    # Check number of values is ok
                    vals = parts[1:]
                    n_val = const.param_num_vals(key)
                    if len(vals) != n_val:
                        story = f"Parameter '{key}' should have {n_val} values; Found '{vals}'\n"
                        story += f"Problem line {line_n}: '{line}'"
                        self.set_problem(story)
                        return False

                    # List are special; Append, else just collect
                    if key == "sample":
                        sample_list.append(vals)
                    elif key == "map_targ":
                        mtarg_list.append(vals)
                    elif key == "rep_gene_url":
                        gurl_list.append(vals)
                    elif key == "bc_round_set":
                        bc_rset_list.append(vals)
                    elif key == "gfasta":
                        gfasta_list.append(vals)
                    elif key == "fastq_samp_slice":
                        fasta_slice = [int(v) for v in vals]
                    else:
                        self._par_pfraw[key] = vals
                INFILE.close()

                # Copy raw parameter values into real par collection
                for key, vals in self._par_pfraw.items():
                    self._set_par_val(key, vals[0])

                # Handle any collected sample or map_targ; Append anything new
                if sample_list:
                    self.get_par_val("sample", as_list=True).extend(sample_list)
                if mtarg_list:
                    self.get_par_val("map_targ", as_list=True).extend(mtarg_list)
                # Add geno urls
                if gurl_list:
                    self.get_par_val("rep_gene_url", as_list=True).extend(gurl_list)
                # Fasta slice list; Just add as list
                if fasta_slice:
                    self._set_par_val("fastq_samp_slice", fasta_slice)
                # Barcode sets; Just append here, Fix name duplicates in barcode setup
                if bc_rset_list:
                    self.get_par_val("bc_round_set", as_list=True).extend(bc_rset_list)
                if gfasta_list:
                    self.get_par_val("gfasta", as_list=True).extend(gfasta_list)
            except Exception:
                story = f"Problem opening / parsing parameter file: '{parfname}'"
                self.set_problem(story)
                return False

        return True

    def _reset_temp_comarg_params(self):
        """Reset any temp command line parameters to final form;

        Command line args may be different (e.g. shorter) than internal params.
        Logic may also be inverted; i.e. --flag  >--> not flag_long_formg
        This function copies shortened keys to full length keys and deletes short keys

        Some flags may be command line only; e.g. --kit_list, --chem_list, --bc_list
        If the 'long' version is missing, then simply delete these

        Returns nothing
        """
        temp_pars = const.get_spipe_temp_par_dict()
        for t_par, full_par in temp_pars.items():
            # Only set if got something (temp comarg) from command line
            if self._pars[t_par]:
                val = self._pars[t_par]
                # full_par is tuple
                long_par, logic_flag = full_par
                # No long version = just delete
                if logic_flag == "none":
                    long_par = None
                # Invert logic
                if logic_flag == "not":
                    val = not val
                if long_par:
                    self._pars[long_par] = val
            del self._pars[t_par]

    def _save_command_line(self, sep=" "):
        """Get system args into command line strings"""
        # Save both wrapped and non-wrapped
        cl1 = utils.get_sys_command_line(wrap=True)
        cl2 = utils.get_sys_command_line(wrap=False)
        self._command_line = [cl1, cl2]

    def _add_comargs(self, args):
        """Add passed arguments to parameter collection

        args = dict; keys have to be real or temp parameters
        """
        # Temp (com line) to internal parameter dict
        temp_par_dict = const.get_spipe_temp_par_dict()
        # Passed parameters have to be real; Actual or temp
        for par, val in args.items():
            if par in self._pars:
                key = par
            elif par in temp_par_dict:
                key = temp_par_dict[par]
            else:
                self.set_problem(f"Unrecoginzed parameter given: '{par}'")
                return False

            # Only set if some value or int (may be zero)
            #   Else if set all keys, will overwrite defaults "", None, False, etc
            if val or isinstance(val, int):
                # print(f"PPP yes setting {key} {val} {type(val)}")
                self._set_par_val(key, val)
            # else:
            #    print(f"PPP NOT setting {key} {val} {type(val)}")

    def _set_par_val(self, key, val, story=None):
        """Set parameter value based on key"""
        # If given any story, val formatted as string with appended story
        if story:
            val = f"{val}\t{story}"
        # Standardize key
        key = const.clean_param_key(key)
        self._pars[key] = val
        # print("VVV", key, self._pars[key])

    def _init_filepath_dict(self):
        """Set up file path format string dict

        Strings are used by 'filepath(key, sample)' function
        """
        self._filepath_dict = const.get_spipe_filepath_dict()

    def _init_fnames(self):
        """Change filename parameters to full explicit paths

        Sublibrary output dirs not checked here
        Doesn't create any files/dirs; Just names and maybe check existance
        """
        # Directories
        self._expand_one_par_fnames("output_dir", dir_slash=True, check_in=False)
        self._expand_one_par_fnames("input_dir", dir_slash=True)
        self._expand_one_par_fnames("genome_dir", dir_slash=True)
        self._expand_one_par_fnames("parent_dir", dir_slash=True)
        # Single files
        self._expand_one_par_fnames("fq1", dir_slash=False)
        self._expand_one_par_fnames("fq2", dir_slash=False)
        self._expand_one_par_fnames("samp_list", dir_slash=False)
        self._expand_one_par_fnames("samp_sltab", dir_slash=False)
        self._expand_one_par_fnames("sublib_list", dir_slash=False)
        self._expand_one_par_fnames("targeted_list", dir_slash=False)
        self._expand_one_par_fnames("dge_cell_list", dir_slash=False)
        self._expand_one_par_fnames("crsp_guides", dir_slash=False)
        # Params that are lists of filenames (files or subdirs)
        self._expand_one_par_fnames("fasta", dir_slash=False)
        self._expand_one_par_fnames("genes", dir_slash=False)
        # gfasta is special; List of lists
        gfasta = self.get_par_val("gfasta", as_list=True)
        if gfasta:
            new_gfasta = []
            for genome, fname in gfasta:
                fpath = utils.fname_full_path(fname, dir_slash=False, exists=False)
                new_gfasta.append([genome, fpath])
            self._set_par_val("gfasta", new_gfasta)

    def _expand_one_par_fnames(self, key, dir_slash=False, check_in=True):
        """Expand to full path file(s) associated with given parameter key

        key         = parameter key (not file key)
        dir_slash   = flag to add trailing slash to dir
        check_in    = flag to check if file exists (i.e. legit input)

        Potentially changed par value (str or list of str) is reset

        Returns status
        """
        # Get param as vanilla; no casting stuff
        fname = self.get_par_val(key, as_raw=True)
        if fname:
            # Subdir doesn't have to exist (here, yet)
            fpath = utils.fname_full_path(fname, dir_slash=dir_slash, exists=False)
            # Check existence?
            if check_in:
                if not utils.check_infile(fpath):
                    self.set_problem(f"Bad filepath: {fname}")
                    return False
            self._set_par_val(key, fpath)
        return True

    def _init_runproc_info(self, top_dir=None, comb_sub=True):
        """Init runtime process info"""
        # For mode 'all' don't load; Everything must be specified de novo
        if self.get_mode("all"):
            self.report_run_story(
                "Mode 'all'; Not using any prior run process definitions"
            )
            return

        fpath = self.filepath("PF_RUNPROC_DEF", None, top_dir=top_dir)
        ok, ndict = utils.check_read_json(fpath)
        if ok:
            self._runproc_info = ndict
            self._runproc_source = fpath
            self.report_run_story(f"Loaded run process definitions {fpath}")
        else:
            # If combined, try again with sublib
            if comb_sub and self.is_combine():
                slib_paths = self.get_sublibs(from_par=True)
                # Get runproc from first sublib; No more recurring
                self._init_runproc_info(top_dir=slib_paths[0], comb_sub=False)

    def _init_chemistry(self):
        """Possibly set chemistry from params or definition file"""
        if self.get_mode("mkref"):
            return

        chem = self.get_chemistry()
        if not chem:
            # Chemistry from runproc returns tuples (value, source)
            p_chem, source = self.get_runproc_info(key="chemistry.name")
            if p_chem:
                self._set_par_val("chemistry", p_chem)
                story = f"Set chemistry {p_chem}, {source} (init chemistry)"
                self.report_run_story(story)
        # Must have something, and set chem-specific params
        chem = self.get_chemistry()
        if chem:
            self._set_chem_defaults()
        else:
            self.set_problem("No chemistry set")

    def _init_kit(self):
        """Possibly set kit from params or definition file"""
        chem = self.get_chemistry()

        if self.get_mode("mkref"):
            # kit 1 is custom / special; No need to specify chemistry
            p_kit_s, _ = kits.parse_kit(1)
            self._set_kit(p_kit_s, chem, source="automatic")
            self.report_run_story(f"Kit set to {p_kit_s} (init kit; mkref mode)")
            return

        # Current parameter settings
        p_kit_s, p_kit_n, p_kit_source = self.get_kit()
        # Try to get info from file
        f_kit_s = f_kit_source = ""
        f_kit_dict, source = self.get_runproc_info(key="kit")
        if f_kit_dict and ("kit_name" in f_kit_dict):
            f_kit_s = f_kit_dict["kit_name"]
            # f_kit_source = f_kit_dict['kit_source']
            f_kit_source = source
            self.report_run_story(f"Got kit {f_kit_s}, from {f_kit_source}")

        # If only one (file or param) is set, set both to that one
        if f_kit_s and (not p_kit_s):
            p_kit_s = f_kit_s
        elif p_kit_s and (not f_kit_s):
            f_kit_s = p_kit_s
            f_kit_source = "given"

        # Conflicting kits?
        if p_kit_s != f_kit_s:
            kit = self.get_par_val("kit", as_str=True)
            story = (
                f"Conflicting kit specs: file '{f_kit_s}'; given '{kit}' ({p_kit_s})"
            )
            self.set_problem(story)
            return

        # Get standard name if not empty (p_kit == s_kit at this point)
        if p_kit_s:
            kit_s, k_story = kits.parse_kit(p_kit_s, chem=chem)
            # Non-legit?
            if p_kit_s and (not kit_s):
                story = f"Failed to parse kit: '{p_kit_s}'; {k_story})"
                self.set_problem(story)
            else:
                self._set_kit(kit_s, chem, source=f_kit_source)
                self.report_run_story(f"Kit set to {kit_s}, {f_kit_source} (init kit)")

    def _set_kit(self, kit, chem, source="", src_story=""):
        """Set kit value / source

        Save parsed name, source.
        """
        kit_s, story = kits.parse_kit(kit, chem=chem)
        if kit_s:
            kit_n = kits.kit_num_wells(kit_s, chem=chem)
            if kit_n > 1:
                story = f"({kit_n} wells)"
            else:
                story = f"({kit_n})"
            self._set_par_val("kit", kit_s, story)
        if source:
            self._set_par_val("kit_source", source, story=src_story)

    def _init_barcodes(self):
        """Possibly set barcodes from run proc defs"""
        bc_rounds, source = self.get_runproc_info(key="barcode.bc_round_set")
        if bc_rounds:
            bc_names = bcutils.get_bc_round_list(bc_rounds, names=True)
            # Only want rounds and names; Nothing more (like description)
            bc_clean = []
            for i, name in enumerate(bc_names):
                bc_clean.append([i + 1, name])
            self._set_par_val("bc_round_set", bc_clean, "")
            self.report_run_story(f"Barcode rounds set to {bc_names} (init barcodes)")

    def _set_cur_mode(self, mode):
        """Set current (run time step) mode flag"""
        if mode:
            assert self.is_legit_mode(mode), f"Non-legit mode: {mode}"
            self._cur_mode = self.parse_mode(mode)

    def get_cur_mode(self):
        """Get current (run time step) mode"""
        return self._cur_mode

    def _init_mode(self):
        """Set run mode"""
        # Get whatever is set / or not
        raw_mode = self.get_par_val("mode", as_str=True)
        if not raw_mode:
            mode = "<Not set>"
        else:
            mode = self.parse_mode(raw_mode)
        # If good, save and set other things
        if self.is_legit_mode(mode):
            self._set_par_val("mode", mode)
            # If combine mode, set fresh too
            if self.get_mode("comb"):
                self.report_run_story("Combine mode; Setting reuse False")
                self._set_par_val("reuse", False)
        else:
            self.set_problem(f"Bad mode specified: '{raw_mode}' >--> '{mode}'")

    def _init_logfname(self):
        """Set logfile name; If none, cook up from version"""
        if not self._logfname:
            logfname = self.get_version().replace(" ", "_").replace(".", "_") + ".log"
            outdir = self.get_par_val("output_dir", as_path=True)
            if outdir:
                if outdir[-1] != "/":
                    outdir = outdir + "/"
                logfname = outdir + logfname
        self._logfname = logfname

    def _init_sys_info(self):
        """Get system info; RAM, CPUs"""
        self._sysram = psutil.virtual_memory().total
        self._ncpu = psutil.cpu_count()

    def _init_nthreads(self):
        """Set nthreads to legit value (i.e. bound 1-ncpu)"""
        # Get as-is here in case fraction; Cast to int below
        n = self.get_par_val("nthreads", as_raw=True)
        n_max = self.get_par_val("nthreads", as_int=True)
        # Save given for reporting
        given_n = n
        ncpu = self.get_sys_cpu()
        # Not set?
        if not n:
            if n_max:
                n = min(n_max, ncpu)
            else:
                n = ncpu
        # Set
        else:
            # Fractional? Get any fractional part
            ok, num = utils.str_to_num(n, as_float=True)
            if ok and num:
                frac = n - int(n)
                if frac > 0:
                    n = round(frac * ncpu)

        # Final number 1 <= int <= CPU
        n = int(min(max(n, 1), ncpu))
        # Limit to max?
        lim_story = ""
        max_t = self.get_par_val("max_threads", as_int=True)
        if max_t:
            if n > max_t:
                n = max_t
                lim_story = f" limited to {max_t}"
        story = f"Setting threads to {n} (ncpu={ncpu}, initial nthreads={given_n}{lim_story}) (init threads)"
        self.report_run_story(story)
        # Automatic?
        story = "" if given_n else "(automatic)"
        self._set_par_val("nthreads", n, story)

    def _init_rand_seed(self):
        """Set random seed(s)"""
        rseed = self.get_par_val("rseed", as_int=True)
        if not rseed:
            self.report_run_story("No random seed given; Creating one")
            # Four digits from time
            t = utils.get_time(as_str=False)
            rseed = int(str(hash(t))[-4:])
            self._set_par_val("rseed", rseed, "(automatic)")
        # Both random and numpy have different random generators
        # https://stackoverflow.com/questions/11526975/set-random-seed-programwide-in-python/11527011
        random.seed(rseed)
        np.random.seed(rseed)

    def _init_in_fastqs(self):
        """Set input fastq files (if any)"""
        # From parameter first
        fq1 = self.get_par_val("fq1", as_path=True)
        fq2 = self.get_par_val("fq2", as_path=True)
        if fq1 and (not fq2):
            fq2 = utils.get_fq2_from_fq1(fq1)
            if fq2:
                self.report_run_story(f"No fq2; Setting from fq1: {fq2} (init fastqs)")
                self._set_par_val("fq2", fq2, "(automatic)")

        # If not set, load from run process
        # Only need to check fq1; If loading from run process, use both fq1 and fq2
        if not fq1:
            p_fq1, _ = self.get_runproc_info(key="fastq.fq1_name")
            p_fq2, _ = self.get_runproc_info(key="fastq.fq2_name")
            if p_fq1:
                self._set_par_val("fq1", p_fq1)
                self.report_run_story(f"Set fastq fq1 {p_fq1} (init fastqs)")
                self._set_par_val("fq2", p_fq2)
                self.report_run_story(f"Set fastq fq2 {p_fq2} (init fastqs)")

    def _init_sample_list(self):
        """If given files for sample definition, parse and set those

        List sorted for any meta samples

        Only puts sample definitions into parameter string; Nothing is processed
        """
        if self.get_mode("mkref"):
            return
        if self.get_mode("comb"):
            # Set plate dims from sublib info
            sub_dirs = self.get_sublibs(from_par=True)
            self._set_wellplates(sublib_list=sub_dirs)
            self.report_run_story("Combine mode; Samples will be set from sublibs")
            return

        new_list = []
        start_list = self.get_par_val("sample", as_list=True)
        what_str = also_str = n_story = ""

        # Possible sources of sample defs
        runproc_dict, source = self.get_runproc_info(key="sample")
        parent_dir = self.get_par_val("parent_dir", as_path=True)
        samp_list = self.get_par_val("samp_list", as_path=True)
        samp_sltab = self.get_par_val("samp_sltab", as_path=True)
        source_dict = {"samp_list": samp_list, "samp_sltab": samp_sltab}
        # Run process info
        if runproc_dict:
            what_str = "Run procss definition info"
            also_str = self._samp_also_source_story(source_dict)
            new_list, n_story = self._init_sample_from_runproc(runproc_dict)
        # Parent dir
        elif parent_dir:
            what_str = f"Parent dir {parent_dir}"
            also_str = self._samp_also_source_story(source_dict)
            new_list, n_story = self._init_sample_from_top_dir(parent_dir)
        # Simple list
        elif samp_list:
            what_str = f"Sample list {samp_list}"
            # Not using samp_sltab inputs, so get this story
            source_dict = {"samp_sltab": samp_sltab}
            also_str = self._samp_also_source_story(source_dict)
            new_list, n_story = self._init_sample_file_list(samp_list)
        # SampleLoadingTable
        elif samp_sltab:
            what_str = f"SampleLoadingTable {samp_sltab}"
            # Not using samp_list inputs, so get this story
            source_dict = {"samp_list": samp_list}
            also_str = self._samp_also_source_story(source_dict)
            new_list, n_story = self._init_sample_sltab_list(samp_sltab)
        # Check if problems (and bail early)
        if self.any_problem():
            self.set_problem(f"Sample init failed for {what_str}")
            return

        # Warn if skipped some inputs
        if also_str:
            self.report_warning(also_str)

        if new_list and start_list:
            story = f"Replacing ({len(start_list)}) given sample definitions"
            self.report_warning(story)

        # Nothing new
        if not new_list:
            new_list = start_list

        # Make sure names are ok
        name_list = []
        samp_def_list = []
        for i, samp_def in enumerate(new_list):
            name = samp_def[0]
            new_name, story = spsamp.clean_sample_def_name(name)
            # If name change, update def and report
            if new_name != name:
                samp_def[0] = new_name
                self.report_warning(story)
            # Make sure unique
            if new_name in name_list:
                story = f"Sample [{i+1}] name '{new_name}' ('{name}') is not unique"
                self.set_problem(story)
            name_list.append(new_name)
            samp_def_list.append(samp_def)

        # Sort to handle meta samples (also just alphabetic)
        ok, samp_def_list, story = spsamp.order_sample_list(samp_def_list)
        if not ok:
            self.set_problem(f"Sample init failed to parse / order")
            self.set_problem(story)
            return

        # Plate dims; Maybe multi plate
        self._set_wellplates(samp_def_list=samp_def_list)

        # Save def list and tell story
        self._set_par_val("sample", samp_def_list)
        if samp_def_list:
            self.report_run_story(f"Set {len(samp_def_list)} samples, {n_story}")
        else:
            self.report_run_story("No samples set")

    def _set_wellplates(self, samp_def_list=None, sublib_list=None):
        """Set well plate (dimension) info

        samp_def_list = list of sample well specs
        sublib_list = list of sublib dirs (from which to get sample_def_list)
        """
        if sublib_list is not None:
            # Cook up sample def list for first sublib; All should be same (re dims)
            sub_dir = sublib_list[0]
            pdict, source = self.get_runproc_info(key="sample.samples", top_dir=sub_dir)
            samp_def_list = []
            for sdict in pdict:
                samp_def_list.append(sdict["wspec"])
        else:
            # Sample def list; Assume [name, well-spec] per sample
            new_list = []
            for samp_def in samp_def_list:
                new_list.append(samp_def[1])
            samp_def_list = new_list
            #    ok, _, well_list_list, story = plates.parse_multi_plate96_wells(wspec)

        # Round 1 kit-specific plate dims
        kit, chem = self.get_kit_chem()
        k_n_rows, k_n_cols, _ = kits.kit_num_rows_cols_plates(kit, chem)

        # List of WellPlate structs
        self._wplates = [plates.WellPlate(n_rows=k_n_rows, n_cols=k_n_cols)]
        # Multi plate?
        if spsamp.is_multi_plate_wspec(samp_def_list):
            # Default plate is 96 well
            self._wplates.append(plates.WellPlate())

    def _samp_also_source_story(self, source_dict):
        """Create string with any non-null sources in dict

        Return str
        """
        also_str = ""
        also_list = []
        for k, v in source_dict.items():
            if v:
                also_list.append(k)
        if also_list:
            also_str = "Not using sample specs from " + ",".join(also_list)
        return also_str

    def _init_sample_from_runproc(self, pdict):
        """Parse samples from runproc info dict

        Return tuple (samples[[name,well],...], story)
        """
        samp_defs = []
        for sdict in pdict["samples"]:
            # Meta samples have sample list
            if "sample_list" in sdict:
                wspec = ",".join(sdict["sample_list"])
            else:
                wspec = sdict["wspec"]
            samp_defs.append([sdict["name"], wspec])
        story = "from runproc_info (init sample_runproc)"
        return samp_defs, story

    def _init_sample_from_top_dir(self, fpath):
        """Parse samples from given path

        Return tuple (samples[[name,well],...], story)
        """
        samp_defs = []
        story = ""
        pdict, source = self.get_runproc_info(key="sample", top_dir=fpath)
        if pdict:
            samp_defs, story = self._init_sample_from_runproc(pdict)
        return samp_defs, story

    def _init_sample_file_list(self, fname):
        """Parse sample definitions from simple file

        Return tuple (samples[[name,well],...], story)
        """
        samp_defs = []
        story = ""
        # Load file lines; Skip blanks. Skip comments is default
        f_lines = utils.lines_from_fname(fname, all_lines=False)
        for line in f_lines:
            parts = line.split()
            # Def should have <name> <wells> [<fraction>]
            if 2 <= len(parts) <= 3:
                samp_defs.append(parts)
            else:
                story = f"Sample spec not as <name> <wells>\n{line}"
                self.set_problem(story)
                samp_defs = []
                break
        if not samp_defs:
            story = f"Failed initializing samples from list {fname}"
            self.set_problem(story)
            samp_defs = []
        else:
            story = f"from list {fname} (init sample_list)"
        return samp_defs, story

    def _init_sample_sltab_list(self, fname):
        """Parse sample definitions from SampleLoadingTable excel file

        Return tuple (samples[[name,well],...], story)
        """
        self.report_run_story(f"Parsing SampleLoadingTable {fname}")
        samp_defs = []
        story = ""
        try:
            sltab_info, story = sltab.parse_sampleloadingtable(fname)
            self.report_run_story(story)
            if sltab_info:
                # List of per-sample info
                for samp in sltab_info["samples"]:
                    parts = [samp["name"], samp["wells"]]
                    samp_defs.append(parts)
                story = f"from list {fname} (init sample_sltab)"
            else:
                self.set_problem("Bad SampleLoadingTable")
        except Exception as e:
            self.set_problem(f"Bad SampleLoadingTable: {e}")
        return samp_defs, story

    def _init_sublib_list(self):
        """If given sublib_list file, parse paths to sublibrary list

        Only puts listed subdirs into parameter string; Nothing is processed
        """
        # If given sublib list file, add these to collection
        fname = self.get_par_val("sublib_list", as_path=True)
        if fname:
            # Load possible list of paths; Don't want blanks. Skip comments is default
            slib_paths = utils.lines_from_fname(fname, all_lines=False)
            if not slib_paths:
                self.set_problem(f"Failed to read sublib_lis: '{fname}'")
                return
            self._set_par_val("sublibraries", slib_paths)

    def _init_combine(self):
        """Check sublibraries are good if combine mode"""
        if not self.get_mode("comb"):
            return

        sub_dirs = self.get_sublibs(from_par=True)
        # Combine needs at least two subdirs
        if len(sub_dirs) < 2:
            self.set_problem(
                f"Combine mode requires >= 2 sublibraries; Only {len(sub_dirs)} given"
            )
            return

        # Make sure subdirs exist
        ok = utils.check_infile(sub_dirs, verb=True)
        if not ok:
            bad_list = utils.bad_infile_list(sub_dirs)
            self.set_problem(f"Sublibrarie(s) are bad: {bad_list}")
            return

        # Make sure all appear to be normal output dirs
        for subdir in sub_dirs:
            o_dict = spout.spoutdir_paths(subdir)
            if not o_dict:
                story = f"Sublibrary is not pipeline output dir: {subdir}"
                self.set_problem(story)
                return

        # Use case
        use_case, source = self.get_runproc_info(key="run.use_case")
        if not use_case:
            self.set_problem(f"Failed to get use_case from sublibrarie(s); {source}")
        self._set_par_val("use_case", use_case)
        story = f"Use_case behavior from sublibs; Set to {use_case} (init combine)"
        self.report_run_story(story)

    def _init_genome(self):
        """Possibly set genome stuff from definition"""
        if self.get_mode("mkref"):
            return

        # Current parameter
        p_dir = self.get_par_val("genome_dir", as_path=True)
        genome_dir = ""

        # Genome dir from file
        pdict, source = self.get_runproc_info(key="genome")
        if pdict and ("genome_dir" in pdict):
            # Check if any given genome_dir exists and matches
            genome_dir = pdict["genome_dir"]
            if p_dir:
                geno_dirs = [genome_dir, p_dir]
                if not utils.check_infile(geno_dirs, verb=False):
                    bad_list = utils.bad_infile_list(geno_dirs)
                    self.set_problem(f"Genome dir(s) bad: {bad_list}")
                elif not os.path.samefile(p_dir, genome_dir):
                    story = (
                        f"Given {p_dir} -vs- loaded {genome_dir} genome paths differ"
                    )
                    self.set_problem(story)

        # Genome param not set? Set from file if available
        if (not p_dir) and genome_dir:
            self._set_par_val("genome_dir", genome_dir)
            self.report_run_story(f"Genome set to {genome_dir} (init genome)")

    def _init_tcr(self):
        """Possibly set TCR status from file"""
        # Use case from file
        pdict, source = self.get_runproc_info(key="run")
        if pdict and ("use_case" in pdict):
            use_case = pdict["use_case"]
            if use_case.startswith("tcr"):
                self._set_par_val("tcr_analysis", True)
                self.report_run_story("TCR analysis true (init tcr)")

    def _init_parent(self):
        """Possibly set parent from definition

        Also may set other things from parent
        """
        if self.get_mode("mkref"):
            return

        # Possibly set from runproc if not given
        parent_dir = self.get_par_val("parent_dir", as_path=True)
        pdir, source = self.get_runproc_info(key="parent.path")
        if pdir and (not parent_dir):
            self._set_par_val("parent_dir", pdir, "")
            story = f"Set parent_dir via runproc_info {pdir} (init parent)"
            self.report_run_story(story)

        # Chemistry and kit from parent (if have one)
        if parent_dir:
            # Chemistry
            p_chem, _ = self.get_runproc_info(key="chemistry.name", top_dir=parent_dir)
            if not p_chem:
                self.set_problem(f"No chemistry found for parent {parent_dir}")
                return
            chem = self.get_par_val("chemistry", as_str=True)
            if chem:
                if not kits.same_chemistry(chem, p_chem):
                    self.report_warning(
                        f"Replaced given chemistry {chem} with {p_chem}"
                    )

            self._set_par_val("chemistry", p_chem, story="(from parent)")
            story = f"Set chemistry {p_chem} via parent {parent_dir} (init parent)"
            self.report_run_story(story)

            # Kit
            p_kit, _ = self.get_runproc_info(key="kit.kit_name", top_dir=parent_dir)
            if not p_kit:
                self.set_problem(f"No kit found for parent {parent_dir}")
                return
            kit = self.get_par_val("kit", as_str=True)
            if kit:
                if not kits.same_kit(kit, p_kit):
                    self.report_warning(f"Replaced given kit {kit} with {p_kit}")

            self._set_par_val("kit", p_kit, story="(from parent)")
            story = f"Set kit {p_kit} via parent {parent_dir} (init parent)"
            self.report_run_story(story)

    def _init_targen(self):
        """Possibly set target enrichment stuff from definition"""
        if self.get_mode("mkref"):
            return

        # Target path not given
        targeted_list = self.get_par_val("targeted_list", as_path=True)
        if not targeted_list:
            f_tpath, source = self.get_runproc_info(key="target.target_source")
            if f_tpath:
                self._set_par_val("targeted_list", f_tpath)
                self.report_run_story(
                    f"Target enrichment list set {f_tpath} (init targen)"
                )

        # Need if combine and targeted and not set
        targeted_list = self.get_par_val("targeted_list", as_path=True)
        if self.is_combine() and self.is_targeted() and (not targeted_list):
            # Source from first sublib; May not exist yet, so from_par=True
            top_dir = self.get_sublibs(from_par=True)[0]
            fpath = self.filepath("PF_TARGET_GENES", None, top_dir=top_dir)
            self._set_par_val("targeted_list", fpath)
            self.report_run_story(f"Target enrichment list set {fpath} (init targen)")

    def _init_fb_crispr(self):
        """Possibly set guides from definition"""
        if self.get_mode("mkref"):
            return
        if not self.is_focal_crispr():
            return

        # If combine, set guides from first sublib
        guide_fname = self.get_par_val("crsp_guides", as_path=True)
        if self.is_combine() and (not guide_fname):
            # May not exist yet, so from_par=True
            top_dir = self.get_sublibs(from_par=True)[0]
            p_fpath = self.filepath("PF_ALL_GUIDES", None, top_dir=top_dir)
            # Make sure file exists
            fpath = utils.check_infile(p_fpath, verb=False)
            if fpath:
                self._set_par_val("crsp_guides", fpath)
                story = f"Set guides from {fpath} (init targen)"
                self.report_run_story(story)
            else:
                story = f"Guides file is bad: {p_fpath}"
                self.set_problem(story)

    def _init_gene_urls(self):
        """Set up format strings for gene urls"""
        # If anything, should be list of [name, url] lists
        for gurl in self.get_par_val("rep_gene_url", as_list=True):
            self._gene_url_dict[gurl[0]] = gurl[1]

    def _init_dge_vars(self):
        """Set transcript / cell count cutoff vars"""
        # Make sure X-fold factor 1 or more
        gxfold = self.get_par_val("dge_cell_gxfold", as_float=True)
        if gxfold < 1:
            self.report_run_story(f"Cell target X-fold {gxfold} limited to 1")
            gxfold = 1
            self._set_par_val("dge_cell_gxfold", gxfold, story="(Limited to 1)")

        cell_targ = self.get_par_val("dge_cell_given", as_int=True)
        cell_est = self.get_par_val("dge_cell_est", as_int=True)
        # Cell target sets min and max; Estimate only max
        if cell_targ:
            cell_min = round(float(cell_targ) / gxfold)
            story = f"(Calculated from {cell_targ:.0f} and {gxfold:3.2f}-fold)"
            self._set_par_val("dge_cell_min", cell_min, story=story)
        if cell_targ or cell_est:
            cell_targ = float(cell_targ) if cell_targ else float(cell_est)
            cell_max = round(cell_targ * gxfold)
            story = f"(Calculated from {cell_targ:.0f} and {gxfold:3.2f}-fold)"
            self._set_par_val("dge_cell_max", cell_max, story=story)

    def _init_map_targs(self):
        """Init possible map_targ values; Should come after genome_dir is set"""
        path = self._get_map_targ_geno_path()
        if path:
            mt_list = self.get_par_val("map_targ", as_list=True)
            # Match fastas that may end *.fas or *.fa
            fastas = utils.fnames_for_path(path, regex=r".*\.fas?$")
            for fas in sorted(fastas):
                # Name = file basename without extension
                name = os.path.splitext(os.path.basename(fas))[0]
                # Add [name, path] list to map_targ list
                mt_list.append([name, fas])
            # Anything set?
            if fastas:
                self._set_par_val("map_targ", mt_list)
                story = f"Added {len(fastas)} map_targs from {path}"
                self.report_run_story(story)

    def _init_mkref(self):
        """Init for mkref"""
        if not self.get_mode("mkref"):
            return

        # CRISPR guides
        # At this stage, use case may not be set so check flag and guides directly
        if self.get_par_val("crsp_analysis", as_bool=True) and self.get_par_val(
            "crsp_guides"
        ):
            return

        # Normal case input; Each genome has <name> <genes> <fasta>
        # Numbers must match
        genome_names = self.get_par_val("genome_name", as_list=True)
        gtf_files = self.get_par_val("genes", as_list=True)
        fasta_files = self.get_par_val("fasta", as_list=True)
        nn = len(genome_names)
        ng = len(gtf_files)
        nf = len(fasta_files)
        if not (nn == ng == nf):
            story = f"Genome name {nn}, gtf {ng}, fasta {nf} numbers must match"
            self.set_problem(story)
            return

        # gfasta case; Each entry here has <name> <fasta>
        gfasta = self.get_par_val("gfasta", as_list=True)
        ngf = len(gfasta)
        # Need normal or gfasta to go forward
        if (nn + ngf) < 1:
            story = "Need [genome-name, genes, fasta] or [gfasta] inputs for mkref"
            self.set_problem(story)
            return

        # Check genome names; Can have dash but not underscore
        ok_story = "(Only alphanumeric and dash)"
        new_gnames = []
        for i, name in enumerate(genome_names):
            n_name = utils.clean_vname_chars(name, to_dash=True, to_under=False)
            if not n_name:
                story = f"Genome name ({i+1}) disallowed format: '{name}' {ok_story}"
                self.set_problem(story)
                return
            elif n_name != name:
                self.report_warning(f"Replaced genome_name {name} with {n_name}")
            new_gnames.append(n_name)
        self._set_par_val("genome_name", new_gnames)

        new_gfasta = []
        for i, (name, fname) in enumerate(gfasta):
            n_name = utils.clean_vname_chars(name, to_dash=True, to_under=False)
            if not n_name:
                story = f"Gfasta name ({i+1}) disallowed format: '{name}' {ok_story}"
                self.set_problem(story)
                return
            elif n_name != name:
                story = f"Replaced gfasta genome name '{name}' with '{n_name}'"
                self.report_warning(story)
            new_gfasta.append([n_name, fname])
        self._set_par_val("gfasta", new_gfasta)

    def _init_run_name(self):
        """Set run name if not set"""
        name = self.get_par_val("run_name", as_str=True)
        if not name:
            name = self.get_par_val("output_dir", as_path=True, dir_slash=False)
            if not name:
                self.set_problem("No output_dir, so cannot get a name")
                return
            name = name.split("/")[-1]
            self._set_par_val("run_name", name, "(automatic)")
            self.report_run_story(f"Set run_name {name} (init run_name)")

    def _init_use_case(self):
        """Set use_case var based on other settings

        Only set case; Doesn't validate that inputs work / are good
        """
        # Special cases
        if self.get_mode("mkref"):
            use_case = self.get_mode()
        else:
            # Will get from sublib if combine mode
            s_use_case, source = self.get_runproc_info(key="run.use_case")
            # Deduce from params if no starting value
            if not s_use_case:
                use_case = modes.decide_use_case(self)
            else:
                use_case = s_use_case

        self._set_par_val("use_case", use_case)
        story = f"Use_case behavior set to {use_case} (init use_case)"
        self.report_run_story(story)

    def _init_run_steps(self):
        """Get list of run steps for current state"""
        full_run_steps = modes.raw_run_step_list(self)
        if self.is_mkref():
            run_steps = full_run_steps
        else:
            # First step is mode, unless 'all'
            first_step = self.get_mode()
            if first_step == "all":
                first_step = full_run_steps[0]

            # Legit starting step?
            if first_step not in full_run_steps:
                use_case = self.get_use_case()
                story = f"Mode '{first_step}' not in use case {use_case} steps {full_run_steps}"
                self.set_problem(story)
                return

            # Truncate starting step?
            start_ind = full_run_steps.index(first_step)
            run_steps = full_run_steps[start_ind:]

            # Truncate end step?
            until_step = self.get_par_val("until_step", as_str=True)
            if self.get_par_val("one_step", as_bool=True):
                run_steps = run_steps[:1]
            elif until_step:
                # Clean up raw parameter
                u_step = modes.parse_mode(until_step)
                if u_step not in run_steps:
                    story = "Problem with specified mode / run steps\n"
                    story += f"Raw step list: {full_run_steps}\n"
                    story += f"Mode '{self.get_mode()}' (start) and until_step '{until_step}' (end) does't work"
                    self.set_problem(story)
                    return
                # Get end with cleaned up mode step. Include until step
                end_ind = run_steps.index(u_step) + 1
                run_steps = run_steps[:end_ind]

        self._run_steps = run_steps

    def _get_map_targ_geno_path(self, mkref=False):
        """Return path go map_targ genome dir (doesn't need to exist)"""
        path = ""
        # If mkref mode, genome_dir is output
        # Get top dirs 'as_path' with trailing slash
        if mkref:
            genome_dir = self.get_par_val("output_dir", as_path=True, dir_slash=True)
        else:
            genome_dir = self.get_par_val("genome_dir", as_path=True, dir_slash=True)
        # Path = top plus any subdir name
        if genome_dir:
            subdir = self.get_par_val(
                "map_targ_genome_dir_path", as_str=True, dir_slash=True
            )
            if subdir:
                path = genome_dir + subdir
        return path

    def _set_up_dirs(self):
        """Check / reset path names
        Set up top-level output dirs; Check / create
        """
        # Must always have output_dir
        top_dir = self.get_par_val("output_dir", as_path=True)
        if not top_dir:
            self.set_problem("No output directory given")
            return False

        # Check output dirs; Create if necessary
        dir_list = [self.filepath("DIR_TOP", None)]
        dir_list.append(self.filepath("DIR_PROC", None))
        dir_list.append(self.filepath("DIR_SCRATCH", None))
        # mkref is special; May have map_to subdir
        if self.get_mode("mkref"):
            # Not if fb/crispr
            if not self.is_focal_crispr():
                path = self._get_map_targ_geno_path(mkref=True)
                if path:
                    dir_list.append(path)
        # focal / crispr processing dir
        if self.is_focal_crispr():
            dir_list.append(self.filepath("DIR_PROC_CRSP", None))

        return utils.make_subdir(dir_list, toxic=False)

    def _set_up_run_env(self):
        """Check run environment (i.e. conda) is good, save info to file"""
        ok = True
        timeout = self.get_par_val("start_timeout", as_int=True)
        self.report_run_story(f"Checking run env (conda); timeout {timeout}")
        if timeout < 1:
            self.report_run_story("Conda list info not checked")
        else:
            # Call 'conda list' both checks and yields outputs
            ver = self.get_version()
            fpath = self.filepath("PF_CONDA_INFO", None)
            utils.file_header(ofile=fpath, story="conda listing", ver=ver)
            # Call, redircting to file
            com_args = "conda list"
            ok, output, ex_story = utils.call_exec(
                com_args, timeout=timeout, verb=False
            )
            if ok:
                with open(fpath, "a") as OUTFILE:
                    print(output, file=OUTFILE)
            else:
                # Log error
                self.set_problem("Subprocess failed: " + com_args)
                self.set_problem(ex_story)
                # Also suggest feedback
                utils.print_story("May need to install / init anaconda")
        return ok

    def _set_up_libs(self):
        """Initialize library-related stuff (e.g. verbosity, logging)"""
        if self.get_mode("mkref"):
            return True

        # Scanpy
        verbosity = self.get_par_val("ana_scanpy_verbosity", as_raw=True)
        utils.init_scanpy_params(verbosity=verbosity, verb=True)
        return True

    def _set_up_execs(self):
        """Check executables / environment is good

        Report save version info / output to files
        """
        timeout = self.get_par_val("start_timeout", as_int=True)
        self.report_run_story(f"Checking called executables; timeout {timeout}")
        if timeout < 1:
            exe_dict = const.get_spipe_exec_dict(call=False)
        else:
            exe_dict = const.get_spipe_exec_dict(call=True, timeout=timeout)

        # Version info repored and saved to file
        ver = self.get_version()
        fpath = self.filepath("PF_VERSION_INFO", None)
        utils.file_header(ofile=fpath, story="called program info", ver=ver)
        # Open as append to follow header
        with open(fpath, "a") as OUTFILE:
            # Unpack each exec; Results = tuple (status, path, version, story)
            for name, xinfo in exe_dict.items():
                # False status means not found / issue
                if not xinfo.status:
                    nf_story = f"'{name}' command not found; May need to install"
                    # Only warn if optional
                    if xinfo.optional:
                        self.report_warning(nf_story)
                        self.report_warning(f"'{name}' is optional; Will skip use")
                    else:
                        self.set_problem(nf_story)
                        self.set_problem(xinfo.story)
                        return False

                # Save lowercase name to obj
                self._execs[name.lower()] = xinfo

                # Report and save to file
                status = "Installed" if xinfo.status else "Not_found"
                self.report_run_story(f"Version info: {xinfo.ver} {xinfo.path}")
                print(name, status, xinfo.ver, xinfo.path, sep="\t", file=OUTFILE)

        return True

    def _set_up_log(self):
        """Open log, handling existing ones"""
        # Get log filename, then bump any existing versions
        fname = self.get_log_filename()
        b = utils.bump_logs(fname)
        story = f"(bumped {b} prior version(s))" if b else ""
        self.report_run_story(f"Opening log {fname} {story}")
        self.log_FH = open(self.get_log_filename(), "w")
        self.report_run_params(ofname=self.log_FH)
        # Mem trace log
        pb_mem.init_mem_profile(self)

    def _set_up_samples(self):
        """Set up internal sample-related stuff"""
        # Don't need if mkref or combine mode
        if self.get_mode("comb,mkref"):
            mode = self.get_mode()
            self.report_run_story(
                f"Set up samples; {mode} mode; Handled elsewhere; Ignoring"
            )
            return True

        self.report_run_story("Setting up samples")
        # Number of rows and cols in (96 well) plate depend on kit, chem
        kit, chem = self.get_kit_chem()
        k_n_rows, k_n_cols, _ = kits.kit_num_rows_cols_plates(kit, chem)
        plate_dims = (k_n_rows, k_n_cols)
        valid_wells = set(self.get_bc_well_list(1))
        # List of sample specification/definitions
        # List is ordered so meta samples are last and ref only prior definitions
        samp_def_list = self.get_par_val("sample", as_list=True)
        samp_dict = {}
        for samp_def in samp_def_list:
            # See if meta sample; Specification is list of prior-defined sample names
            sname, wspec = samp_def[:2]
            m_samples = []
            for part in spsamp.get_meta_wspec_parts(wspec, samp_dict):
                if part not in samp_dict:
                    m_samples = []
                    break
                m_samples.append(samp_dict[part])

            # Meta sample
            if m_samples:
                samp = spsamp.SpSamp(sname, samples=m_samples, plate_dims=plate_dims)
            else:
                # Check that wells are legit for kit
                ok, _, well_list_list, story = plates.parse_multi_plate96_wells(wspec)
                if not ok:
                    self.set_problem(f"Sample {sname} wells |{wspec}| failed parsing")
                    self.set_problem(story)
                    return False
                # Only look at first plate wells
                well_set = set(well_list_list[0])
                bad_set = well_set - valid_wells
                if bad_set:
                    story = f"Sample '{sname}' wells |{wspec}| does't work for kit {kit}; {bad_set}"
                    self.set_problem(story)
                    return False

                samp = spsamp.SpSamp(sname, wspec=wspec, plate_dims=plate_dims)
                # Make sure wells are legit for current plate
                samp_wells = set(samp.get_well_list(plate=1))
                if len(samp_wells - valid_wells) > 0:
                    story = f"Sample specification has wells outside plate: {samp_def}"
                    self.set_problem(story)
                    bc_info = self.get_bc_info()
                    story = bc_info.get_bc_df_story(1)
                    self.set_problem(story)
                    return False
            samp_dict[sname] = samp

        # Check / warn if any well overlaps
        samp_list = list(samp_dict.values())
        for i, samp1 in enumerate(samp_list):
            # Ignore all-well and meta samples
            if samp1.is_allwell() or samp1.is_meta():
                continue
            for j in range(i + 1, len(samp_list)):
                samp2 = samp_list[j]
                if samp2.is_allwell() or samp2.is_meta():
                    continue
                # Overlap well count and story
                n, story = spsamp.sample_well_overlap(samp1, samp2)
                if n:
                    # Not allowed?
                    if self.get_par_val("samp_wells_unique", as_bool=True):
                        self.set_problem(story)
                        return False
                    self.report_run_story2(story)

        # Save collection of samples thus far; get_samples needs this
        self._sample_dict = samp_dict

        # Need all-well?
        need_allwell = self.get_par_val("make_allwell", as_bool=True)
        # If no given samples, need allwell
        if not samp_dict:
            need_allwell = True
        # If already have an allwell sample, don't need new one
        if self.get_samples(non_allwell=False, meta=False):
            need_allwell = False
        if need_allwell:
            n_plates = self.num_sample_plates()
            n_rows, _ = self.get_bc_rows_cols(from_kit=True)
            wspec = spsamp.get_wspec_for_all_well(n_rows, n_plates=n_plates)
            sname = self.get_allwell_sample(as_name=True)
            # If this is the only sample, report
            if not samp_dict:
                self.report_run_story2(
                    f"No samples; Single 'all-well' sample named '{sname}'"
                )
            samp_def = [sname, wspec]
            samp_def_list.append(samp_def)
            samp = spsamp.SpSamp(sname, wspec=wspec, plate_dims=plate_dims)

            samp_dict[sname] = samp
            self.report_run_story(f"Added '{sname}' sample; {samp_def}")

        # all-sample sample?
        if self.get_par_val("make_allsample", as_bool=True):
            # More than one non-allwell and non-meta samples, then make all-sample
            samp_list = self.get_samples(allwell=False, meta=False)
            if len(samp_list) > 1:
                sname = "all-sample"
                parts = ",".join([s.get_name() for s in samp_list])
                samp_def = [sname, parts]
                samp_def_list.append(samp_def)
                samp = spsamp.SpSamp(sname, samples=samp_list, plate_dims=plate_dims)
                samp_dict[sname] = samp
                self.report_run_story(f"Added '{sname}' sample; {samp_def}")
                # Check if doesn't cover all wells
                if not samp.is_allwell():
                    # Each plate
                    p_story = ""
                    for plate in samp.get_plate_indexes():
                        n_row, n_col = samp.get_plate_dims(plate=plate)
                        tot_wells = n_row * n_col
                        n_wells = samp.num_wells(plate=plate)
                        if n_wells < tot_wells:
                            p_story += (
                                f"; plate {plate}: {n_wells} of {tot_wells} wells"
                            )

                    story = f"'all-sample' sample covers < all wells{p_story}"
                    self.report_warning(story)

        # Save final collection of samples
        self._sample_dict = samp_dict

        # Add well counts (just for reporting)
        new_sample_par = []
        for samp in self.get_samples():
            specs = [samp.get_name(), samp.get_def_str()]
            if samp.is_meta():
                specs.append(samp.get_meta_parts(as_name=True))
            new_sample_par.append(specs)
        self._set_par_val("sample", new_sample_par)

        self.report_run_story(f"Set up {self.num_samples()} samples")
        return True

    def get_allwell_sample(self, as_name=True):
        """Get "the" all-well sample or name

        Return str or obj
        """
        samp = None
        # Attempt to get list of only all-well samples
        samp_list = self.get_samples(non_allwell=False)
        if samp_list:
            samp = samp_list[0]

        # Default name "all-sample" even if no actual samples
        if as_name:
            if samp:
                samp = samp.get_name()
            else:
                samp = "all-sample"
        return samp

    def _set_up_samp_dirs(self):
        """Make sure all sample subdirs exist / are created"""
        ok = True
        # Base case; Always have these
        dir_keys = "DIR_TOP,DIR_PROC,DIR_SAMP,DIR_S_REP".split(",")
        # Saving figures?
        if self.get_par_val("rep_save_figs", as_bool=True):
            # Not for TCR
            if not self.is_tcr():
                dir_keys.append("DIR_S_FIG")
        # Specific DGE outputs depend on use case
        dir_keys += const.use_case_out_dir_key(self, as_list=True)
        # Make sure all dirs exist for all samples
        for samp in self.get_samples():
            for fkey in dir_keys:
                # Explicitly try to make subdir; If fail don't crash, but catch here
                path = self.filepath(fkey, samp, mkdir=False)
                ok = utils.make_subdir(path, toxic=False)
                if not ok:
                    self.set_problem(
                        f"Bad dir for sample |{samp.get_name()}| path |{path}|"
                    )
                    return False
        return ok

    def _get_sample_specs(self):
        """Get specification settings for samples (e.g. for json)

        Return dict
        """
        ndict = {}
        samp_list = self.get_samples()
        sdic_list = []
        ndict["number_samples"] = len(samp_list)
        for samp in samp_list:
            sdic_list.append(spsamp.get_sample_specs(samp))
        ndict["samples"] = sdic_list
        return ndict

    def _set_up_map_targs(self):
        """Set up map target structs; Create, init and save"""
        # Not if combine or mkref mode
        if self.get_mode("comb,mkref"):
            return True

        ok = True
        genome_dir = self.get_par_val("genome_dir", as_path=True)
        mtarg_list = self.get_par_val("map_targ", as_list=True)
        # Create, init and collect any map target objects
        for name, mpath in mtarg_list:
            # Make sure map-target path is good
            # First as is, else then with genonome path
            in_path = utils.check_infile(mpath, verb=False)
            if not in_path:
                in_path = utils.check_infile(mpath, path=genome_dir, verb=False)
            if not in_path:
                self.set_problem(f"map_targ path: |{genome_dir}| file |{mpath}|")
                ok = False
                break

            mapt = maptarg.MapTarg(self, name, in_path)
            mapt.kmer_len = self.get_par_val("map_targ_kmer_len", as_int=True)
            mapt.kmer_step = self.get_par_val("map_targ_kmer_step", as_int=True)
            mapt.kmer_set, mapt.rc_kmer_set = utils.fasta_kmer_sets(
                in_path, kmer_len=mapt.kmer_len, kmer_step=mapt.kmer_step
            )
            # Also keep search parameters
            mapt.kmer_search_max = self.get_par_val(
                "map_targ_kmer_search_max", as_int=True
            )
            mapt.kmer_search_step = self.get_par_val(
                "map_targ_kmer_search_step", as_int=True
            )

            self._mtarg_dict[name] = mapt
        return ok

    def _set_up_mkref(self):
        """Set up for mkref mode ..."""
        # Only mkref mode
        if not self.get_mode("mkref"):
            return True

        path = self.filepath("DIR_TOP", None, mkdir=False)
        ok = utils.make_subdir(path, toxic=False)
        return ok

    def _set_up_kit(self):
        """Set up / score / check kit specific parameters"""
        if self.get_mode("mkref"):
            return True

        # Record of scoring
        score_info = {}
        proc_story = given_story = ""
        all_scores = None

        ok = True
        # Get current kit, chem info
        kit_s, kit_n, kit_source = self.get_kit()
        chem = self.get_chemistry()
        # Score unless combine, focal, skipping, or already scored
        if self.get_mode("comb"):
            proc_story = "Combine mode; No kit scoring"
            self.report_run_story(proc_story)
        elif self.is_focal_crispr():
            proc_story = "Use case focal / crispr; No kit scoring"
            self.report_run_story(proc_story)
        elif self.get_par_val("kit_score_skip", as_bool=True):
            proc_story = "Skipping, no kit scoring"
            self.report_warning("Skipping kit scoring")
        elif kit_source.startswith("scor"):
            proc_story = "Kit scoring already done; Skipping"
            self.report_run_story(proc_story)
        elif kit_source.startswith("loaded runproc_info"):
            proc_story = "Kit scoring from runproc_info; Skipping"
            self.report_run_story(proc_story)
        elif kits.is_custom_kit(kit_s):
            proc_story = "Custom kit; No scoring; Skipping"
            self.report_run_story(proc_story)
        else:
            proc_story = "Score kit barcodes"
            given_story = kit_s if kit_s else "None"

            kit_guess, kit_score, all_scores = kits.guess_kit_from_bc(self)
            # Check score against warning and min
            warn_score = self.get_par_val("kit_score_warn", as_float=True)
            min_score = self.get_par_val("kit_score_min", as_float=True)
            if min_score <= kit_score < warn_score:
                self.report_warning(f"Best kit score is low: {round(kit_score, 3)}")
            if kit_score < min_score:
                story = "Kit score {kit_score} < min {min_score}; Do not trust"
                story += "\n  Can bypass kit score restriction with --kit_score_skip\n"
                self.set_problem(story)
                return False

            # Compare given to best score
            if kit_s:
                if kit_guess == kit_s:
                    self.report_run_story2(f"Given and guessed kit match: {kit_s}")
                else:
                    story = f"Given kit {kit_s} doesn't match guess {kit_guess}; Do not trust"
                    story += (
                        "\n  Can bypass kit score restriction with --kit_score_skip\n"
                    )
                    self.set_problem(story)
                    return False

            # Take the best if no kit given
            else:
                kit_s, _ = kits.parse_kit(kit_guess, chem=chem)
                kit_source = f"score {round(kit_score, 3)}"

        score_info["processing"] = proc_story
        score_info["given_kit"] = given_story
        score_info["kit_bc_scores"] = all_scores
        # Save kit scoring info; Retrieved later, e.g. for runproc info saving
        self._kit_score_info = score_info

        # Make sure have kit
        if kit_s:
            # Set standardized kit info
            self._set_kit(kit_s, chem, source=kit_source)
            # Set kit-specific defaults
            self._set_kit_defaults()
        else:
            self.set_problem("Invalid / no kit")
            # Set last call to report kits
            self.add_last_call("kits.report_valid_kits()")
            ok = False

        return ok

    def _get_run_specs(self):
        """Get define settings for run (e.g. for json)

        Return dict
        """
        ndict = {}
        ndict["use_case"] = self.get_use_case()
        ndict["mode"] = self.get_mode()
        ndict["run_name"] = self.get_par_val("run_name", as_str=True)
        ndict["output_dir"] = self.get_par_val("output_dir", as_path=True)
        return ndict

    def _get_chemistry_specs(self):
        """Get define settings for chemistry (e.g. for json)

        Return dict
        """
        chem = self.get_chemistry()
        return {"name": chem}

    def _get_kit_specs(self, score_list=None):
        """Get define settings for kit (e.g. for json)

        Return dict
        """
        kit_s, kit_n, kit_source = self.get_kit()
        new_dict = {"kit_name": kit_s, "kit_num": kit_n, "kit_source": kit_source}
        new_dict["kit_bc_scoring"] = self._kit_score_info
        return new_dict

    def _set_chem_defaults(self):
        """Set chem-specific default values"""
        chem = self.get_chemistry()
        par_dict = const.get_chem_default_pars(chem)
        # R2 read amp structure
        if not self.get_par_val("bc_amp_seq", as_str=True):
            self._set_par_val("bc_amp_seq", par_dict["bc_amp_seq"], "")
            story = f"Set bc_amp_seq for chemistry {chem} (set chem defaults)"
            self.report_run_story(story)

    def _set_kit_defaults(self):
        """Set kit-specific default values"""
        kit, chem = self.get_kit_chem()
        par_dict = const.get_kit_default_pars(kit, chem)

        # Nothing if custom
        if kits.is_custom_kit(kit):
            story = "Kit custom; Not setting default params, bc, etc (set kit defaults)"
            self.report_run_story(story)
        else:
            # Since parameter is a list, don't set story (not supported)
            self._set_par_val("bc_round_set", par_dict["bc_round_set"], "")
            story = f"Set barcodes for kit {kit} chemistry {chem} (set kit defaults)"
            self.report_run_story(story)
            # Set cell max if not already set
            if not self.get_par_val("dge_cell_max", as_int=True):
                self._set_par_val("dge_cell_max", par_dict["dge_cell_max"], "")
                story = f"Set dge_cell_max for kit {kit} chemistry {chem} (set kit defaults)"
                self.report_run_story(story)

    def _set_up_barcodes(self):
        """Set up barcode related structures

        Load seqs (pandas series) and similar-edit-dics
        Create seq >--> index dicts

        """
        # Don't need if mkref mode
        if self.get_mode("mkref"):
            return True

        # Barcode specification; List[num-bc-round] of lists[bc-key, bc-name,..]
        bc_rounds = self.get_par_val("bc_round_set", as_list=True)
        # If nothing set, set from chem and kit
        if not bc_rounds:
            chem = self.get_chemistry()
            kit_s, _, _ = self.get_kit()
            bc_rounds = kits.kit_bc_set_list(kit_s, chem, as_tup=True)
            self._set_par_val("bc_round_set", bc_rounds, "")
            self.report_run_story(f"Set barcodes for kit {kit_s} chemistry {chem}")
            self.report_run_story2(f"Barcode rounds {bc_rounds}")

        # Get bc round list again (in case just set)
        bc_rounds = self.get_par_val("bc_round_set", as_list=True)
        path = self.get_pkg_path()
        amp_seq = self.get_par_val("bc_amp_seq", as_str=True)
        cell_thresh = self.get_par_val("pre_perfect_bc_cell_thresh", as_float=True)
        bc_info, story = bcutils.get_bc_info(
            bc_rounds, path, amp_seq, cell_thresh=cell_thresh
        )
        if not bc_info:
            self.set_problem(story)
            self.set_problem("Failed to load barcode data")
            return False

        # Make sure we have all barcodes set
        all_rounds = set([1, 2, 3])
        given_rounds = set(bcutils.get_bc_round_list(bc_rounds, as_int=True))
        if given_rounds < all_rounds:
            self.set_problem(
                f"Need barcodes {sorted(all_rounds)} specified; Found {bc_rounds}"
            )
            return False

        # Update barcode parameters for reporting
        # Extend lists from [bcX, bc_name] to [bcX, bc_name, story]
        for r, rlist in enumerate(bc_rounds):
            # Index are 1-based
            story = bc_info.get_bc_df_story(r + 1)
            rlist.append(story)

        # Attach
        self._bc_info = bc_info

        return True

    def _set_up_genome(self):
        """Genome working files (if dir given)"""
        # Don't need if mkref, TCR
        if (
            self.get_mode("mkref")
            or self.is_tcr()
            or self.is_focal_crispr(map_direct=True)
        ):
            return True

        # Legit path?
        genome_dir = self.get_par_val("genome_dir", as_path=True)
        if not utils.check_infile(genome_dir):
            self.set_problem(f"Genome dir not found: '{genome_dir}'")
            return False

        fname = self.filepath("GENO_GENE_INFO", None)
        if not utils.check_infile(genome_dir):
            self.set_problem(f"Genome data file not found: '{fname}'")
            self.set_problem("Reference genome maybe needs to be rebuilt (mkref)")
            return False

        # Check min version?
        min_ver = self.get_par_val("align_min_genome_version", as_str=True)
        if min_ver:
            # Attempt to get genome version; Return value, source (if found)
            g_ver, _ = self.get_runproc_info(key="header.ver_number", mkref=True)
            if not utils.ver_at_least_ver(g_ver, min_ver):
                self.set_problem(
                    f"Genome {genome_dir} version {g_ver} < min required {min_ver}"
                )
                self.set_problem("Reference genome needs to be rebuilt (mkref)")
                return False
            self.report_run_story(
                f"Genome {genome_dir} version {g_ver} >= min {min_ver}"
            )

        # Attempt to load genome info; Takes a bit of time, so feed back
        self.report_run_story("Setting up genome files")
        ok = True

        try:
            # Load gene info
            gene_info = genes.load_gene_info(fname)
            self._set_gene_info(gene_info)
        except Exception:
            self.set_problem(f"Failed to set up genome; Dir: {genome_dir}")
            ok = False

        # Copy all_gene file (always, in case genome changed)
        fpath = self.filepath("PF_ALL_GENES", None)
        if ok:
            in_path = self.filepath("GENO_GENE_DATA", None)
            story = utils.copy_file(in_path, fpath)
            if not story:
                self.set_problem(f"Failed to copy gene data {in_path} to {fpath}")
                self.set_problem("Reference genome maybe needs to be rebuilt (mkref)")
                ok = False
            else:
                self.report_run_story(story)
        return ok

    def _set_gene_info(self, gene_info):
        """Set passed gene info into structure"""
        assert isinstance(gene_info, dict)

        gene_info["genome_dir"] = self.get_par_val("genome_dir", as_path=True)

        # Dataframe with all gene info
        gene_id_to_name = gene_info["gene_id_to_name"]
        gene_id_to_genome = gene_info["gene_id_to_genome"]
        all_genes_df = pd.DataFrame()
        all_genes_df["gene_id"] = sorted(list(gene_info["gene_id_to_name"].keys()))
        all_genes_df["gene_name"] = all_genes_df["gene_id"].apply(
            lambda s: gene_id_to_name[s]
        )
        all_genes_df["genome"] = all_genes_df["gene_id"].apply(
            lambda s: gene_id_to_genome[s]
        )
        gene_info["all_genes"] = all_genes_df
        gene_info["genome_list"] = sorted(list(all_genes_df["genome"].unique()))

        # Dict mapping gene id to all_genes index
        gene_id_to_idx = {v: k for k, v in dict(all_genes_df["gene_id"]).items()}
        gene_info["gene_id_to_idx"] = gene_id_to_idx

        self._gene_info = gene_info

    def _get_genome_specs(self):
        """Get specification settings for genome (e.g. for json)

        Return dict
        """
        gene_info = self.get_gene_info()
        if not gene_info:
            gene_info = self.get_guide_info()
            if not gene_info:
                return {}
        ndict = genes.get_genome_specs(gene_info)
        return ndict

    def _get_fastq_specs(self):
        """Get specification settings for fastq files

        Return dict
        """
        fq1 = self.get_par_val("fq1", as_path=True)
        fq2 = self.get_par_val("fq2", as_path=True)

        ndict = {}
        if fq1:
            ndict["fq1_name"] = fq1
            ndict["fq1_data"] = utils.report_file_info_str(fq1, name=False)
        if fq2:
            ndict["fq2_name"] = fq2
            ndict["fq2_data"] = utils.report_file_info_str(fq2, name=False)
        return ndict

    def _set_up_parent(self):
        """Parent dir"""
        # Don't need if mkref mode (Maybe more modes?)
        if self.get_mode("mkref"):
            return True

        # No list = nothing to setup
        path = self.get_par_val("parent_dir", as_path=True)
        if not path:
            return True

        # Try to get obj
        obj = spout.spoutdir_obj(path)
        if not obj:
            story = f"Bad parent dir {path}"
            self.set_problem(story)
            return False

        self.report_run_story(f"Processing parent dir {path}")
        # Get cell barcodes
        ok = False
        try:
            self.report_run_story2(f"Reading barcodes in parent {path}")
            bc_dict, story = spout.load_top_dir_cell_bcs(path)
            if bc_dict is None:
                pass
            else:
                self.report_run_story(story)
                # Check number of barcodes
                if len(bc_dict["all_bc_set"]) < 1:
                    story = f"Too few cell barcodes found: {len(bc_dict['all_bc_set'])}"
                    self.set_problem(story)
                    return False

                # Save info and set barcodes
                self._parent_info = {
                    "obj": obj,
                    "path": path,
                }
                self._cell_bc = bc_dict
                ok = True
        except Exception as e:
            story = f"Problem getting parent barcodes: {e} (exception)"
            self.set_problem(story)
        return ok

    def get_genome_dir(self):
        """Get current genome_dir path"""
        # For mkref mode, output_dir is what is usually genome_dir
        if self.get_mode("mkref"):
            genome_dir = self.get_par_val("output_dir", as_path=True)
        else:
            genome_dir = self.get_par_val("genome_dir", as_path=True)
        return genome_dir

    def get_parent_info(self, key=None):
        """Get parent info dict or value of one key

        Return dict or key value
        """
        pdict = self._parent_info
        if pdict and key:
            pdict = pdict.get(key)
        return pdict

    def have_parent(self, par_only=True):
        """Do we have a parent dir?"""
        have = False
        if self._parent_info is not None:
            have = True
        if par_only and self.get_par_val("parent_dir", as_str=True):
            have = True
        return have

    def _get_parent_specs(self):
        """Get define settings for parent (e.g. for json)

        Return dict
        """
        ndict = {}
        parent_info = self.get_parent_info()
        if parent_info is not None:
            ndict["path"] = parent_info["path"]
        return ndict

    def get_cell_bc(self, samp):
        """Get cell barcode set for sample

        Return set of well index triples
        """
        if self.have_cell_bc():
            if isinstance(samp, spsamp.SpSamp):
                samp = samp.get_name()
            return self._cell_bc[samp]

    def have_cell_bc(self):
        return bool(self._cell_bc is not None)

    def _set_up_dge_cell_list(self):
        """Set up DGE cell list

        Needs to be after samples are set up

        Return status
        """
        if self.get_mode("mkref"):
            return True

        # No list = nothing to setup
        fname = self.get_par_val("dge_cell_list", as_path=True)
        if not fname:
            return True

        if self.have_cell_bc():
            story = "Cannot have both parent_dir and cell_list"
            self.set_problem(story)
            return False

        ok = False
        try:
            self.report_run_story2(f"Reading cell list barcodes {fname}")
            bc_df = utils.read_csv(fname, head_if="bc_wells")
            if bc_df is None:
                pass
            else:
                story = f"Loaded {len(bc_df)} rows from {fname}"
                self.report_run_story(story)
                if self._set_cell_bc(bc_df):
                    ok = True
        except Exception as e:
            story = f"Bad cell list barcodes: {e} (exception)"
            self.set_problem(story)
        return ok

    def _set_cell_bc(self, bc_df):
        """Set cell list for use"""
        new_dict = {}
        # Create sets for each sample except all-well
        samp_set_dict = {}
        for samp in self.get_samples():
            if samp.is_allwell():
                continue
            # Create set, save for sample and fill with r1 well indexes
            samp_set = set()
            new_dict[samp.get_name()] = samp_set
            for r1 in samp.get_wind_list():
                samp_set_dict[r1] = samp_set

        # Set of all barcodes; Should be as str e.g. '25_60_03'
        all_bc_set = set(bc_df["bc_wells"].values)
        aw_name = self.get_allwell_sample(as_name=True)
        new_dict[aw_name] = all_bc_set
        for bc in all_bc_set:
            r1wind = int(bc.split("_")[0])
            samp_set_dict[r1wind].add(bc)

        self._cell_bc = new_dict
        return True

    def _set_up_targen(self):
        """Target enrichment"""
        # Don't need if mkref mode (Maybe more modes?)
        if self.get_mode("mkref"):
            return True

        # No list = nothing to setup
        fname = self.get_par_val("targeted_list", as_path=True)
        if not fname:
            return True

        # Get dataframe with target genes
        ok = False
        try:
            self.report_run_story2(f"Reading target enrichment list file {fname}")
            tgen_df = utils.read_csv(fname, head_if="gene_id,gene_name")
            if tgen_df is None:
                pass
            else:
                story = f"Loaded {len(tgen_df)} rows from {fname}"
                self.report_run_story(story)
                if self._set_targen_info(tgen_df, source=fname):
                    ok = True
        except Exception as e:
            story = f"Problem setting up target enrichment: {e} (exception)"
            self.set_problem(story)
        return ok

    def _set_targen_info(self, tgen_df, source=""):
        """Set target enrichment info from target dataframe"""
        gene_info = self.get_gene_info()
        all_genes = gene_info["all_genes"]
        # Get target column for matching
        if "gene_id" in tgen_df.columns.values:
            ag_col = "gene_id"
        elif "gene_name" in tgen_df.columns.values:
            ag_col = "gene_name"
        else:
            # Use values in first column to attempt match in all_genes
            tgen_list = list(tgen_df.iloc[:, 0])
            ag_col = utils.guess_df_col_by_values(all_genes, tgen_list)

        # Get list of targets from dataframe
        tgen_list = list(tgen_df[ag_col])
        story = f"Got {len(tgen_list)} target enrichment genes (column {ag_col})"
        self.report_run_story(story)

        # Get set; For matching but also to check unique
        tgen_set = set(tgen_list)
        if len(tgen_set) != len(tgen_list):
            story = f"Target enrichment genes not unique: {len(tgen_list)} vs {len(tgen_set)}"
            self.set_problem(story)
            return False

        # Get all_gene subset with target genes; Make sure all inputs found
        tgenes_df = all_genes[all_genes[ag_col].isin(tgen_set)]
        tgen_got_set = set(tgenes_df[ag_col])
        missing_set = tgen_set - tgen_got_set
        if missing_set:
            story = f"Target enrichment genes missing: Got {len(tgen_got_set)} vs {len(tgen_set)}"
            self.set_problem(story)
            story = f"{len(missing_set)} missing genes: {missing_set}"
            self.set_problem(story)
            return False

        # Reset target_genes index 0 - N (no longer maps to all_genes)
        tgenes_df.index = range(len(tgenes_df))
        gene_id_to_idx = {v: k for k, v in dict(tgenes_df["gene_id"]).items()}
        # Save things; Main thing gene_id in set
        gene_id_list = list(tgenes_df["gene_id"])
        targen_info = {
            "gene_id_list": gene_id_list,
            "target_number": len(gene_id_list),
            "in_targ_source": source,
            "in_allgene_col": ag_col,
            "in_targ_list": tgen_list,
            "target_genes": tgenes_df,
            "gene_id_to_idx": gene_id_to_idx,
        }
        self._targen_info = targen_info
        return True

    def _get_targen_specs(self):
        """Get define settings for target enrichment (e.g. for json)

        Return dict
        """
        ndict = {}
        targen_info = self.get_targen_info()
        if targen_info is not None:
            ndict["target_source"] = targen_info["in_targ_source"]
            ndict["target_number"] = targen_info["target_number"]
            ndict["in_allgene_col"] = targen_info["in_allgene_col"]
        return ndict

    def _write_targen_data(self, force=False):
        """Write target enrichement data to output file"""
        targen_info = self.get_targen_info()
        if targen_info is not None:
            # Only write data if option and file doesn't exist
            fpath = self.filepath("PF_TARGET_GENES", None)
            if force or (not utils.check_infile(fpath, verb=False)):
                self.write_df(
                    targen_info["target_genes"], "PF_TARGET_GENES", index=False
                )

    def _set_up_fb_crispr(self):
        """focal barcoding / crispr"""
        if not self.is_focal_crispr():
            return True
        # if self.get_mode("comb,mkref"):
        if self.get_mode("mkref"):
            return True

        if not self.get_par_val("parent_dir", as_path=True):
            self.set_problem("Focal barcoding / crispr requires parent_dir")
            return False

        guide_info = crispr.load_guide_info(self)
        if not guide_info:
            self.set_problem("Focal barcoding / crispr guide issue(s)")
            return False
        self._set_guide_info(guide_info)

        # Pepare 'genome' for use
        if not crispr.prep_guide_genome(self):
            self.set_problem("Focal barcoding / crispr requires guides")
            return False

        return True

    def _set_guide_info(self, guide_info):
        """Set guide info structure"""
        self._guide_info = guide_info

    def _write_fb_crispr_data(self, force=False):
        """Write fp crispr data to output file"""
        guide_info = self.get_guide_info()
        if guide_info:
            # Only write data if option or file doesn't exist
            fpath = self.filepath("PF_ALL_GUIDES", None)
            if force or (not utils.check_infile(fpath, verb=False)):
                self.write_df(guide_info["guide_df"], "PF_ALL_GUIDES", index=False)

    def _get_fb_crispr_specs(self):
        """Get settings for focal / crsipr (e.g. for json)

        Return dict
        """
        new_dict = {}
        guide_info = self.get_guide_info()
        if guide_info:
            new_dict = crispr.get_guide_info_specs(guide_info)
        return new_dict

    def _set_up_tcr(self):
        """focal barcoding / crispr"""
        if not self.is_tcr():
            return True
        if self.get_mode("comb,mkref"):
            return True

        ok = True
        # Make sure database installed
        if tcr.have_igblast_dbs(self, verb=True):
            self.report_run_story("TCR database(s) found")
        else:
            self.set_problem("TCR database(s) not found; Need to install")
            ok = False

        return ok

    def _get_barcode_specs(self):
        """Get define settings for combine mode (e.g. for json)

        Return dict
        """
        ndict = {}
        bc_rounds = self.get_par_val("bc_round_set", as_list=True)
        ndict["bc_round_set"] = bc_rounds
        return ndict

    def _set_up_combine(self):
        """Set up for combine mode ...

        Sets (SpOut) sublibrary structs
        Also sets up collection of (combine-mode SpSamp) samples

        """
        # Only combine mode
        if not self.get_mode("comb"):
            return True

        # Check number of sublibs (kit specific)
        if not kits.kit_ok_combine_sublib_count(self):
            return False

        self.report_run_story("Setting up for combine mode")
        # Get dict with unqiue names for each subdir
        sub_dirs = self.get_sublibs(from_par=True)
        sub_dict = spout.get_sublib_name_dict(sub_dirs)
        # Collect and set sublib objects
        for sdir_path, name in sub_dict.items():
            # Get obj
            slib = spout.spoutdir_obj(sdir_path, name=name)
            if not slib:
                self.set_problem(f"Sublibrary dir not valid: '{sdir_path}'")
                return False
            # Set into dict
            self._sublib_dict[name] = slib

        # Make sure all sublib outputs are 'equal' for combine purposes
        sublibs = self.get_sublibs()
        # Compare first to rest
        f_out = sublibs[0]
        # TCR has no genome, so no compare
        # Also currently samples must be identical across sublibs for TCR
        comp_genome = False if self.is_tcr() else True
        comp_samp = True if self.is_tcr() else False

        for s_out in sublibs[1:]:
            same, story = spout.compare_spoutdirs(
                f_out, s_out, genome=comp_genome, samples=comp_samp
            )
            if not same:
                self.set_problem("Sublibraries incompatible for combine")
                self.set_problem(f"First sublib: '{f_out.get_path()}'")
                self.set_problem(f"Other sublib: '{s_out.get_path()}'")
                self.set_problem(story)
                return False

        # Collect set of all unique sample names
        samp_name_set = set()
        for s_out in sublibs:
            for sname in s_out.get_samples(as_name=True):
                samp_name_set.add(sname)

        # Compare collected sample names vs any given sample names
        given_samples = self.get_par_val("sample", as_list=True)
        if given_samples:
            story = f"Ignoring {len(given_samples)} given sample defs"
            self.report_warning(story)
        self.report_run_story(f"Sublibs have {len(samp_name_set)} unique sample names")

        # Number of rows and cols in (96 well) plate depend on kit, chem
        kit, chem = self.get_kit_chem()
        k_n_rows, k_n_cols, _ = kits.kit_num_rows_cols_plates(kit, chem)
        plate_dims = (k_n_rows, k_n_cols)

        # Create samples from names
        samp_dict = {}
        for name in sorted(samp_name_set):
            # combine mode sample has only collection of subdirs; No wells or meta-parts
            samp = spsamp.SpSamp(name, sublibs=sublibs, plate_dims=plate_dims)
            samp_dict[name] = samp
        self.report_run_story(f"Samples set up, {len(sublibs)} sublibs deep")

        # Meta samples have recipe but need to have actual sample parts set
        for samp in samp_dict.values():
            if samp.is_simple():
                continue
            # Collect parts
            part_names = samp.get_meta_parts(as_name=True)
            samp_list = [samp_dict[p] for p in part_names]
            samp.set_meta_part_samples(samp_list)

        # Set struct and parameter story; Meta after simple
        new_sample_par = []
        for name, samp in samp_dict.items():
            self._sample_dict[name] = samp
            # Save parameter (just for reporting)
            if samp.is_meta():
                recipe = samp.get_meta_parts(as_name=True)
            else:
                recipe = samp.get_wspec()
            n_subs = samp.num_sublibs(min_wells=1)
            story = f"(Combine mode; {n_subs} of {len(sub_dirs)} sublibraries)"
            new_sample_par.append([name, recipe, story])
        self._set_par_val("sample", new_sample_par)
        return True

    def _get_combine_specs(self):
        """Get define settings for combine mode (e.g. for json)

        Return dict
        """
        ndict = {}
        if self.get_mode("comb"):
            sublibs = self.get_sublibs()
            ndict["sublib_number"] = len(sublibs)
            for i, slib in enumerate(sublibs):
                ndict[f"sublib_s{i+1}_name"] = slib.get_name()
                ndict[f"sublib_s{i+1}_path"] = slib.get_path()
        return ndict

    # ------------------------ public methods -------------------------------
    def get_command_line(self, wrap=True):
        """Return command line string"""
        if wrap:
            comline = self._command_line[0]
        else:
            comline = self._command_line[1]
        return comline

    def get_version(self, numonly=False):
        """Return version string

        numonly = flag to extract only number; e.g. "abc v.1.2" >--> "1.2"
        """
        ver = self._version
        if numonly:
            ver = re.search(r"\s*([\d.]+)", ver).group(1)
        return ver

    def parse_mode(self, mode):
        """Attempt to parse given mode into cannonical string"""
        mode_str = modes.parse_mode(mode)
        return mode_str

    def is_legit_mode(self, mode=None):
        """Check if currently set or given mode is legit"""
        if not mode:
            mode = self.get_mode()
        legit = bool(self.parse_mode(mode))
        return legit

    def get_mode(self, mode=None):
        """Get mode -OR- check mode

        If mode given, compare this to current value and return True/False
        If mode is a list of modes, return True/False if current mode in list
        When no mode given, return current mode value
        """
        rval = self.get_par_val("mode", as_str=True)
        # Parse to standardized form
        rval = self.parse_mode(rval)
        if mode:
            if isinstance(mode, str):
                if "," in mode:
                    mode = mode.split(",")
            if isinstance(mode, list):
                rval = any([self.get_mode(m) for m in mode])
            else:
                mode = self.parse_mode(mode)
                rval = self.is_legit_mode(mode=mode) and (mode == rval)
        return rval

    def get_use_case(self, use_case=None, top_dir=None):
        """Get use_case -OR- check it if passed"""
        # Given top_dir, get value from that
        if top_dir:
            rval, source = self.get_runproc_info(key="run.use_case", top_dir=top_dir)
        else:
            rval = self.get_par_val("use_case", as_str=True)
        # If given use_case, check this
        if use_case:
            if isinstance(use_case, str):
                if "," in use_case:
                    use_case = use_case.split(",")
            if isinstance(use_case, list):
                rval = any([self.get_use_case(m, top_dir=top_dir) for m in use_case])
            else:
                rval = use_case == rval
        return rval

    def get_run_steps(self, as_call=False, as_name=False, as_tup=False):
        """Get list of run steps

        as_call = return list of call functions
        as_name = return list of step names
        as_tup = return list with tuples (call, name)
        """
        step_list = self._run_steps
        # More than simple list of keys?
        if as_tup:
            as_call = as_name = True
        if as_call or as_name:
            step_list = modes.run_step_to_call_name(
                step_list, get_call=as_call, get_name=as_name
            )
        return step_list

    def get_runproc_info(self, key=None, top_dir=None, comb_sub=True, mkref=False):
        """Get run process dict or value of one key

        key = possible info key to extract; Can be multi.level via dots
        top_dir = output dir source; Can be path or SpOut struct
        comb_sub = Flag to look at subdir if combine mode
        mkref = Flag to get info from genome dir mkref file

        Return tuple (dict or key value, source of data)
        """
        # print(f"get_runproc_info >>> key={key} top={top_dir} comb_sub={comb_sub}")
        pdict = None
        source = ""
        from_sub = False
        # Mkref = genome dir file
        if mkref:
            genome_dir = self.get_par_val("genome_dir", as_path=True)
            if genome_dir:
                fpath = self.filepath("GENO_MKREF_DEF", None, genome_dir=genome_dir)
                pdict = utils.read_json(fpath)
                source = genome_dir
        else:
            # Given path
            if top_dir:
                # print(f"get_runproc_info ++ top_dir={top_dir}")
                pdict = {}
                if isinstance(top_dir, spout.SpOut):
                    pdict = top_dir.get_runproc_info()
                    source = top_dir.get_path()
                    # print(f"get_runproc_info + SpOut")
                elif isinstance(top_dir, str):
                    pdict = spout.spoutdir_runproc(top_dir)
                    source = top_dir
                    # print(f"get_runproc_info + str")
            # No path given
            else:
                if comb_sub and self.is_combine() and (not mkref):
                    # print(f"get_runproc_info ++ no top_dir, combine")
                    # Source from first sublib; May not exist yet, so from_par=True
                    top_dir = self.get_sublibs(from_par=True)[0]
                    pdict, source = self.get_runproc_info(
                        key=key, top_dir=top_dir, comb_sub=False
                    )
                    from_sub = True
                else:
                    # print(f"get_runproc_info ++ no top_dir, not combine")
                    pdict = self._runproc_info
                    source = "loaded runproc_info"

        # If nothing, set to empty dict
        # if not pdict:
        #    pdict = {}

        # Keys to limit output; Not if from sublib call
        # print(f"get_runproc_info ++ pdict type={type(pdict)} key={key}")
        if pdict and key and (not from_sub):
            # If we've got pdict from sublib via recur, it may already be a string
            if isinstance(pdict, dict):
                # This handles nested.keys (e.g. 'run.use_case')
                pdict = utils.get_nested_dict_val(pdict, key)

        # print("get_runproc_info ++ pdict", pdict)
        # what = 'dict' if isinstance(pdict, dict) else pdict
        # print(f"get_runproc_info <<< pdict={what}   source={source}")
        # print()

        return (pdict, source)

    def is_legit_kit(self, kit=None):
        """Check if currently set or given kit is legit"""
        if kit:
            chem = self.get_chemistry()
            kit_s, _ = kits.parse_kit(kit, chem=chem)
        else:
            kit_s, kit_n, _ = self.get_kit()
        # Ok if parse returned string
        legit = bool(kit_s)
        return legit

    def get_kit(self, g_kit=None):
        """Get kit int well number, str kit name, and source string

        If not given a kit to check, return current kit values (str, int, source)

        return tuple or bool
        """
        p_kit = self.get_par_val("kit", as_str=True)
        chem = self.get_chemistry()
        kit_s, story = kits.parse_kit(p_kit, chem=chem)

        # If given kit, check against this
        if g_kit:
            rval = False
            if kit_s:
                # Attempt to parse if int
                if isinstance(g_kit, int):
                    g_kit, story = kits.parse_kit(g_kit, chem=chem)
                # Split str with commas, else attempt to parse
                if isinstance(g_kit, str):
                    if "," in g_kit:
                        g_kit = g_kit.split(",")
                    else:
                        g_kit, story = kits.parse_kit(g_kit, chem=chem)
                # Match anything in list, else directly
                if isinstance(g_kit, list):
                    rval = any([self.get_kit(k) for k in g_kit])
                else:
                    rval = g_kit == kit_s
        # Get tuple of parts
        else:
            kit_n = 0
            if kit_s:
                kit_n = kits.kit_num_wells(kit_s, chem=chem)
            source = self.get_par_val("kit_source", as_raw=True)
            rval = (kit_s, kit_n, source)

        return rval

    def get_chemistry(self, as_raw=False, as_num=False):
        """Get chemistry

        as_raw = flag to return raw current value -vs- parsed
        """
        chem = self.get_par_val("chemistry", as_str=True)
        if not as_raw:
            chem = kits.parse_chemistry(chem, as_num=as_num)
        return chem

    def get_kit_chem(self, as_num=False):
        """Get kit and chemistry versions

        as_num = flag to return as numbers -vs- as string

        Return tuple kit, chem
        """
        kit_s, kit_n, _ = self.get_kit()
        chem = self.get_chemistry(as_num=as_num)
        if as_num:
            kit_s = kit_n
        return kit_s, chem

    def get_header_dict(self, desc="", fkey=""):
        """Get dict with header info (for json files)"""
        ver = self.get_version()
        vernum = self.get_version(numonly=True)
        ndict = utils.get_header_dict(ver=ver, vernum=vernum, desc=desc, fkey=fkey)
        return ndict

    def get_sys_mem(self, human=True):
        """Return system memory
        human = flag to convert memory to human-readable string
        """
        mem = self._sysram
        if human:
            mem = utils.report_bytes_str(mem)
        return mem

    def get_sys_cpu(self):
        """Return system num cpu"""
        return self._ncpu

    def get_pkg_path(self, dir_slash=False):
        """Return path to (python) package"""
        path = utils.fname_full_path(self._pkg_path, dir_slash=dir_slash)
        return path

    def get_par_filename(self):
        """Return parameter file name"""
        return self._parfname

    def get_log_filename(self):
        """Return log file name"""
        return self._logfname

    def get_log_fh(self):
        """Return log file handle"""
        return self.log_FH

    def get_par_dict(self):
        """Return parameter dictionary"""
        return self._pars

    def get_par_val(
        self,
        key,
        as_int=False,
        as_float=False,
        as_bool=False,
        as_str=False,
        as_path=False,
        dir_slash=False,
        as_list=False,
        as_raw=False,
    ):
        """Return parameter value based on key

        as_int = flag to return int
        as_float = flag to return float
        as_bool = flag to parse param as bool; Not just str exists; True, False, etc
        as_str = flag to return str
        as_path = flag to return as file path
        dir_slash = flag to make path end with slash (e.g. for dir/)
        as_raw = flag to return raw dict value, whatever that is

        If any type flags set (i.e. as_?), only first token returned (ignore any story)
        """
        c_par = const.clean_param_key(key)
        val = self._pars[c_par]
        # Raw value?
        if as_raw:
            return val
        # If val is string, get only first word
        if val and isinstance(val, str):
            val = val.split()[0]
        # Per type
        if as_int:
            val = int(val) if val else 0
        elif as_float:
            val = float(val) if val else 0.0
        elif as_bool:
            # Parse returns None if val is empty or None
            _, val = utils.str_to_bool(val)
            val = bool(val)
        elif as_path:
            if val and len(val) > 0:
                val = utils.fname_full_path(val, dir_slash=dir_slash, exists=False)
        elif as_str:
            val = str(val) if val else ""
        elif as_list:
            if val:
                assert isinstance(
                    val, list
                ), f"Expected list for key {key} got {type(val)}"
            else:
                val = []
        return val

    def get_par_story(self, key):
        """Get any story associated with a parameter key

        Return string
        """
        val = self.get_par_val(key, as_raw=True)
        story = ""
        if val:
            try:
                # Story is any text after a first space
                story = val.split(" ", 1)[1:]
            except Exception:
                pass
        return story

    def get_par_slice(self, key, lower=False):
        """Get slice parameters first : last

        key = slice parameter key
        lower = flag for reported story to be lowercase

        Return tuple (first, last, story)
        """
        samp_first, samp_last = self.get_par_val(key, as_list=True)
        story = utils.report_range_str(samp_first, samp_last, lower=lower)
        return (samp_first, samp_last, story)

    def filepath(self, fkey, sample_name, top_dir="", genome_dir="", mkdir=True):
        """Return filepath for dir/file indicated by fkey and sample_name

        fkey = Filename keyword
        samp_name = Sample name -OR- SpSamp object
        top_dir = Top level (i.e. output) dir
        genome_dir = Genome dir; Can be real or keyword to get from parameter
        mkdir = Flag to make directory if doesn't exist

        Return path
        """
        # Top dir not given; use current output_dir
        if not top_dir:
            top_dir = self.get_par_val("output_dir", as_path=True)
            # or current dir
            if not top_dir:
                top_dir = "."
        if top_dir.endswith("/"):
            top_dir = top_dir[:-1]

        # Genome_dir not given, use current
        if not genome_dir:
            genome_dir = self.get_par_val("genome_dir", as_path=True)
        else:
            # Given may be a parameter key rather than real path
            if genome_dir in self.get_par_dict().keys():
                genome_dir = self.get_par_val(genome_dir, as_path=True)
        if not genome_dir:
            genome_dir = ""
        if genome_dir.endswith("/"):
            genome_dir = genome_dir[:-1]

        # If given SpSamp, get string name
        if sample_name:
            if isinstance(sample_name, spsamp.SpSamp):
                sample_name = sample_name.get_name()
        if not sample_name:
            sample_name = ""

        # Get actual formatted string
        fpath = self._filepath_dict[fkey].format(
            top_dir=top_dir, sample_name=sample_name, genome_dir=genome_dir
        )
        # Directories are special
        if fkey.startswith("DIR_"):
            # Make sure has trailing slash
            if not fpath.endswith("/"):
                fpath = fpath + "/"
            # Making sure it exists?
            if mkdir:
                # This will try to create if needed; Crash on failure
                utils.make_subdir(fpath, toxic=True)
        return fpath

    def fileparts(self, fkey, sample_name, top_dir=""):
        """Return tuple of file path, name for file indicated by fkey and sample_name"""
        fpath = self.filepath(fkey, sample_name, top_dir=top_dir, mkdir=False)
        name = os.path.basename(fpath)
        path = os.path.dirname(fpath)
        return path, name

    def filename(self, fkey, sample_name, top_dir=""):
        """Return only file name (no path) for given key"""
        return self.fileparts(fkey, sample_name, top_dir=top_dir)[1]

    def get_sublibs(self, one="", as_name=False, from_par=False, raw_path=True):
        """List of sublibraries

        Get list or single sublib object or path

        from_par = flag to return path from parameters only
        fullpath = flag to get / check the fullpath of parameter-only path

        Return single obj or list of obj or str
        """
        answer = None
        # From parameters?
        if from_par:
            slibs = self.get_par_val("sublibraries", as_list=True)
            pref = self.get_par_val("sublib_pref", as_path=True, dir_slash=True)
            suff = self.get_par_val("sublib_suff", as_str=True)
            if one:
                if one in slibs:
                    # Full path from given plus pref and suff, then clean path has to exist
                    answer = pref + one + suff
                    if not raw_path:
                        answer = utils.fname_full_path(answer)
            else:
                answer = []
                for sub_lib in slibs:
                    fpath = self.get_sublibs(one=sub_lib, from_par=True)
                    answer.append(fpath)
        # Else object collection
        else:
            s_dict = self._sublib_dict
            # Single object
            if one:
                if one in s_dict:
                    answer = s_dict[one]
            # Whole list of names or objects
            else:
                if as_name:
                    answer = list(s_dict.keys())
                else:
                    answer = list(s_dict.values())

        return answer

    def num_sublibs(self, from_par=True):
        """Number of sublibraries"""
        return len(self.get_sublibs(from_par=from_par))

    def is_combine(self, par_only=True):
        """Combine mode"""
        # Combine mode should have libraries, but initially just param
        is_comb = False
        if self.num_sublibs() > 1:
            is_comb = True
        if par_only and self.get_par_val("sublibraries", as_str=True):
            is_comb = True
        return is_comb

    def is_mkref(self):
        """Mkref mode"""
        return self.get_mode("mkref")

    def get_samples(
        self,
        one="",
        as_name=False,
        allwell=True,
        non_allwell=True,
        meta=True,
        simple=True,
        sort=False,
        sim_met_all=True,
    ):
        """Get sample structs (not input sample spec strings)

        one         = name of single sample to get
        as_name     = flag to return names rather than objects
        allwell     = flag to qualify all-well samples (may be more than 'all-well')
        non_allwell = flag to qualify non-all-well samples
        meta        = flag to qualify meta samples
        simple      = flag to qualify simple (non-meta) samples
        sort        = flag to sort; Alphabetic and meta-reference (so meta after parts)
        sim_met_all = flag to sort so simple then meta then all(sample or well)

        Return one sample or list
        """
        answer = None
        # Single object
        if one:
            answer = self._sample_dict.get(one)
        # Whole list of objects
        else:
            samp_list = list(self._sample_dict.values())
            # Sort on name?
            if sort:
                samp_list = sorted(samp_list, key=lambda x: x.get_name())
            # Sort so simple then meta then all-X ?
            if sim_met_all:
                samp_list = spsamp.sort_samp_list_sim_met_all(samp_list)

            # Now qualify on flags
            answer = []
            for samp in samp_list:
                keep = True
                # Meta or simple
                if samp.is_meta():
                    if not meta:
                        keep = False
                else:
                    if not simple:
                        keep = False
                # All well or not
                if samp.is_allwell():
                    if not allwell:
                        keep = False
                else:
                    if not non_allwell:
                        keep = False
                if keep:
                    answer.append(samp)

            # As name only?
            if as_name:
                answer = [s.get_name() for s in answer]
        return answer

    def num_samples(
        self, allwell=True, non_allwell=True, meta=True, simple=True, is_good=False
    ):
        """Number of sublibraries

        First flags qualify samples; See get_samples()
        is_good = flag to only count samples with good status

        Return int
        """
        n = 0
        for samp in self.get_samples(
            allwell=allwell, non_allwell=non_allwell, meta=meta, simple=simple
        ):
            if is_good:
                if not samp.get_status(is_good=True):
                    continue
            n += 1
        return n

    def num_good_samples(self):
        """Get number of good samples

        Return num_good, num_toal, story
        """
        n_tot = self.num_samples()
        n_good = self.num_samples(is_good=True)
        story = f"Sample status: {n_tot} total, {n_good} good"
        return n_good, n_tot, story

    def get_plate_dims(self, plate=1):
        """Get plate row, col number

        plate = 1-based plate index; If 0, return list

        Return tuple or list of tuples
        """
        if plate == 0:
            answer = [wp.get_dims() for wp in self._wplates]
        else:
            answer = self._wplates[plate - 1].get_dims()
        return answer

    def num_sample_plates(self):
        """Number of plates for sample defs"""
        return len(self.get_plate_dims(plate=0))

    def get_map_targ_list(self):
        """List of map target structs"""
        return list(self._mtarg_dict.values())

    def num_map_targs(self):
        """Number of defined map targets"""
        return len(self.get_map_targ_list())

    def get_gene_url(self, genome, gene_id):
        """Get URL for given genome gene"""
        # Get format string
        try:
            fstr = self._gene_url_dict[genome]
        except Exception:
            fstr = self.get_par_val("rep_gene_def_url", as_str=True)
        # Get complete url
        try:
            url = fstr.format(gene_id=gene_id)
        except Exception:
            url = ""
        return url

    def force_fresh_files(self):
        """True if forcing fresh files"""
        return not self.get_par_val("reuse", as_bool=True)

    def keep_temp_files(self):
        """True if temp files should be kept"""
        return self.get_par_val("keep_temps", as_bool=True)

    def keep_going(self):
        """True if should keep going even if problems"""
        return self.get_par_val("keep_going", as_bool=True)

    def resume_mode(self):
        """True if not one step mode"""
        return not self.get_par_val("one_step", as_bool=True)

    def compress_temp_files(self):
        """True if temp / intermediate files should be compressed"""
        return not self.get_par_val("no_compress", as_bool=True)

    def set_up_env(self):
        """Do setup-up steps, initialize and check of run-time environ

        If all good to go, return True
        """
        if not self.is_legit_mode():
            self.set_problem("No valid mode specified")

        # Don't check kit yet; May need to be guessed

        # Check problems; Could have multiple sources above(?)
        if self.any_problem():
            self.set_problem("Cannot set up env with problems found")
            return False

        # Settings just for reporting
        use_case = self.get_use_case()
        mode = self.get_mode()
        chem = self.get_chemistry()
        kit_s, kit_n, story = self.get_kit()
        if not kit_s:
            kit_s = "not set"
        story = f"Setting up run-time env; Use_case {use_case}, Mode {mode}, Kit {kit_s} Chem {chem}"
        self.report_run_story(story)

        ok = True
        # Set up functions to call; Order may matter
        call_lis = [
            "_set_up_dirs",
            "_set_up_run_env",
            "_set_up_libs",
            "_set_up_execs",
            "_set_up_kit",
            "_set_up_barcodes",
            "_set_up_fb_crispr",
            "_set_up_tcr",
            "_set_up_genome",
            "_set_up_parent",
            "_set_up_targen",
            "_set_up_combine",
            "_set_up_samples",
            "_set_up_samp_dirs",
            "_set_up_dge_cell_list",
            "_set_up_map_targs",
            "_set_up_mkref",
        ]
        for call in call_lis:
            # print(f"# setup {call}")
            if not eval(f"self.{call}")():
                self.set_problem(f"Set up function '{call}' failed")
                ok = False
                break
        # Open log file (last part of set_up_env; only create if all good up to here)
        if ok:
            self._set_up_log()
        # Return set_up_env status
        return ok

    def write_env_files(self, force=True):
        """Write env files

        force = Flag to write runproc info even if file already exists
        """
        # Don't save anything for mkref here (it has it's own save);
        if self.get_mode("mkref"):
            return True

        if self.any_problem():
            if force:
                print(f"Problems found; Still saving env with force={force}")
            else:
                print("Problems found; Not saving env")
                return False

        # Run process info file
        fkey = "PF_RUNPROC_DEF"
        fname = self.filepath(fkey, None)
        # Possibly bumping file versions?
        b = utils.bump_logs(fname)
        story = f"(bumped {b} prior version(s))" if b else ""
        o_dict = self._get_current_runproc_specs()
        utils.write_json(o_dict, fname)
        self.report_run_story(f"Saved env info {fname} {story}")
        # Samples
        for samp in self.get_samples():
            ofname = self.filepath("SFR_SAMP_DEF", samp)
            spsamp.write_def(samp, ofname=ofname)

        # Other data
        bcutils.write_bc_user_data(self)
        self._write_targen_data()
        self._write_fb_crispr_data()

    def _get_current_runproc_specs(self, fkey=""):
        """Get current state into runproc info dict

        Return nested dict
        """
        # Collect dict of info for json
        o_dict = {}
        o_dict["header"] = self.get_header_dict(desc="run process info", fkey=fkey)
        # Defs for specific things
        call_funcs = {
            "run": "_get_run_specs",
            "chemistry": "_get_chemistry_specs",
            "kit": "_get_kit_specs",
            "fastq": "_get_fastq_specs",
            "genome": "_get_genome_specs",
            "sample": "_get_sample_specs",
            "parent": "_get_parent_specs",
            "target": "_get_targen_specs",
            "fb_crispr": "_get_fb_crispr_specs",
            "barcode": "_get_barcode_specs",
            "combine": "_get_combine_specs",
        }
        # Call each to get dict with info
        for key, call in call_funcs.items():
            o_dict[key] = eval(f"self.{call}")()
        return o_dict

    def close_down(self):
        """All done, finish up"""
        # Any last things to call
        self.run_last_call()

        # Final story and time
        print("# Closing down SplitPipe")
        utils.close_out_fh(self.get_log_fh())
        # Final time
        t = utils.get_time(stime=self._time0, as_str=True)
        print(f"# Total time {t}")
        return True

    def set_problem(self, story, report=False):
        """Set / record some problem state described by story"""
        self._prob_count += 1
        # Get info on calling function
        cframe = inspect.stack()[1]
        call_str = (
            f"Function {cframe.function}, file {cframe.filename} line {cframe.lineno}"
        )

        # Append new line if alreay have issues and add call info
        if self._prob_story:
            self._prob_story += "\n"
        self._prob_story += f"Problem: {call_str}\n"

        # Add prefix and keep story
        if isinstance(story, str):
            story = [story]
        for st_line in story:
            self._prob_story += f"Problem: {st_line}\n"

        if report:
            self.report_run_story(story, verb=0, pref="")

    def clear_problems(self, report=True):
        """Clear any problem story. Not problem count"""
        # Report anything already set
        self.any_problem(report=report)
        self._prob_story = ""
        if report:
            self.report_run_story("Cleared any noted problem(s) ... Can keep going")

    def any_problem(self, report=False):
        """Return any problem(s) in story string"""
        if report and self._prob_story:
            self.report_run_story(self._prob_story, verb=0, pref="")
        return self._prob_story

    def no_problems(self):
        """Return True / False if no / any problems are set"""
        return bool(not self.any_problem())

    def problem_count(self):
        """Return cumulative count of set problems"""
        return self._prob_count

    def get_status_story(self, pref="Status:"):
        """Get story of current status

        Return str
        """
        # Problems are counted even for 'keep_going' so could be many
        n_prob = self.problem_count()
        if self.no_problems():
            status = "OK with some issue(s)..." if n_prob else "successfully"
            story = f"{pref}{status} ({n_prob} problems)"
        else:
            story = f"{pref}with some issue(s)... {n_prob} problems"
            story += self.any_problem()
        return story

    def add_last_call(self, func, force=False):
        """Add given function call to list"""
        if (func not in self._last_calls) or force:
            self._last_calls.append(func)

    def run_last_call(self):
        """Call any functions in the last call list (and reset)"""
        for call in self._last_calls:
            eval(call)
        self._last_calls = []

    def report_proc_step(
        self, run_step, status=None, set_mode=True, to_log=False, newline=True
    ):
        """Report processing step; i.e. In/Out high-level step

        status = process status; None = start, True/False = end
        set_mode = flag to set current mode if starting step
        to_log = flag to output to log file
        """
        # Output file; Log or get_out_fh(ofname=None) = stdout
        if to_log:
            OUTFILE = self.get_log_fh()
        else:
            OUTFILE, _ = utils.get_out_fh(ofname=None)

        # output path for disk space reporting
        path = self.filepath("DIR_TOP", None, mkdir=False)
        # Start of process
        if status is None:
            self._s_time = utils.get_time(as_str=False)
            print(f"# >> Running {run_step}", file=OUTFILE)
            self.report_run_mem(story=run_step, to_log=to_log)
            self.report_run_disk(path, story=run_step, to_log=to_log)
            # Mem trace log; Only once, so not if otherwise output to log
            if not to_log:
                story = f"{run_step} start"
                pb_mem.report_to_mem_profile(self, story)
            # Set current mode?
            if set_mode:
                self._set_cur_mode(run_step)
        # End of process
        else:
            self.report_run_mem(story=run_step, to_log=to_log)
            self.report_run_disk(path, story=run_step, to_log=to_log)
            t = utils.get_time(stime=self._s_time, as_str=True)
            print(
                f"# << Running {run_step}; Status {status}; Duration {t}", file=OUTFILE
            )
            if newline:
                print("", file=OUTFILE)
        OUTFILE.flush()

        # Mirror in log file Single recursive call
        if (not to_log) and (self.get_log_fh()):
            self.report_proc_step(
                run_step, status=status, set_mode=set_mode, to_log=True, newline=newline
            )

    def report_run_story(self, story, verb=1, pref="#", to_log=False, no_log=False):
        """Report run-time story ... Intermediate status, feedback, etc.

        Default behavior is to write to stdout and to log via recursion
        If to_log, write to log and don't recur; If no_log, no recursion to log

        to_log = flag to output to log file / init str
        """
        if isinstance(story, list):
            for line in story:
                self.report_run_story(
                    line, verb=verb, pref=pref, to_log=to_log, no_log=no_log
                )
            return

        assert isinstance(story, str), f"Expected str got {type(story)}"

        # Output file; Log or get_out_fh(ofname=None) = stdout
        if to_log:
            OUTFILE = self.get_log_fh()
        else:
            OUTFILE, _ = utils.get_out_fh(ofname=None)
        # Cook up prefix
        if verb > 1:
            vpref = f"{pref}{verb}"
        else:
            vpref = pref
        if vpref:
            vpref = vpref + " "
        o_str = vpref + story
        # If have outfile, print and flush buffer, else just collect
        if OUTFILE:
            print(vpref + story, file=OUTFILE)
            OUTFILE.flush()
        else:
            self._init_story.append(o_str)

        # Mirror in log file; Single recursive call and no flag
        if (not to_log) and (not no_log):
            self.report_run_story(story, verb=verb, pref=pref, to_log=True)

    def report_warning(self, story, pref="# WARNING:", to_log=False):
        self.report_run_story(story, verb=1, pref=pref, to_log=to_log)

    # Wrappers for various verbosity level output
    def report_run_story2(self, story, pref="#", to_log=False, no_log=False):
        self.report_run_story(story, verb=2, pref=pref, to_log=to_log, no_log=no_log)

    def report_run_story3(self, story, pref="#", to_log=False):
        self.report_run_story(story, verb=3, pref=pref, to_log=to_log)

    def report_run_mem(self, story="", no_space=True, to_log=False, verb=2):
        rep_str = "MemStatus:"
        if story:
            if no_space:
                story = story.replace(" ", "_")
            rep_str += f" {story}"
        rep_str += utils.report_mem_str(pref=" ")
        self.report_run_story(rep_str, to_log=to_log, verb=verb)

    def report_run_disk(self, path="", story="", no_space=True, to_log=False, verb=2):
        rep_str = "DskStatus:"
        if story:
            if no_space:
                story = story.replace(" ", "_")
            rep_str += f" {story}"
        if not path:
            path = "/"
        rep_str += utils.report_disk_str(path, pref=" ")
        self.report_run_story(rep_str, to_log=to_log, verb=verb)

    def report_run_params(
        self,
        ofname=None,
        mode="a",
        unset="<Not set>",
        comline=True,
        dump_init=True,
    ):
        """Report run-time parameters

        ofname = output filename / file handle.
            If ofname is a filename or file handle, write to that, else stdout
        mode = file open mode
        unset = string to use for unset parameters

        """
        # Output goes?
        OUTFILE, new_open = utils.get_out_fh(ofname=ofname, mode=mode)

        # System story; Put in dict and dump to keep same format as params
        sdic = OrderedDict()
        sdic["Program"] = self.get_version()
        sdic["Run date"] = utils.get_time(as_str=True)
        node, machine, sysname, uname = utils.get_sys_login_info()
        sdic["Run on"] = f"{node} ({machine} {sysname}) by {uname}"
        sdic["System"] = f"{self.get_sys_mem()} RAM and {self.get_sys_cpu()} CPUs"
        sdic["PID"] = f"{os.getpid()}"
        print("#" * 80, file=OUTFILE)
        dict_str = utils.report_dict_str(sdic, k_fmt="# {k:<25} ", unwrap="Item")
        print(dict_str, file=OUTFILE)
        print("#" * 80, file=OUTFILE)

        # Command line
        if comline:
            print("#", file=OUTFILE)
            print(
                "#                ------ >>> Command line <<< -----------------",
                file=OUTFILE,
            )
            print("#", file=OUTFILE)
            print(self.get_command_line(), file=OUTFILE)

        # Any init info?
        if dump_init and self._init_story:
            print("#", file=OUTFILE)
            print(
                "#                ------ >>> Initialization messages <<< ------",
                file=OUTFILE,
            )
            print("#", file=OUTFILE)
            print("\n".join(self._init_story), file=OUTFILE)

        print("#", file=OUTFILE)
        # Collect values in dict so can replace unset vals if needed
        pdic = OrderedDict()
        pdic["use_case"] = self.get_use_case()
        pdic["mode"] = self.get_mode()
        chem = self.get_chemistry()
        pdic["chemistry"] = chem
        kit_s, kit_n, kit_source = self.get_kit()
        pdic["kit"] = kit_s
        pdic["kit_source"] = kit_source
        pdic["parfile"] = self.get_par_filename()
        pdic["logfile"] = self.get_log_filename()
        for k, v in self.get_par_dict().items():
            pdic[k] = v
        # Replace any unset values
        if unset:
            for k, v in pdic.items():
                if not v:
                    # Only set 'unset' if None, empty string, empty list;
                    # False is not changed
                    if (v is None) or isinstance(v, str) or isinstance(v, list):
                        pdic[k] = unset
        # Dump dict, with given key format and unwrapping label
        print(
            "#                ------ >>> Parameter settings <<< -----------",
            file=OUTFILE,
        )
        print("#", file=OUTFILE)
        dict_str = utils.report_dict_str(pdic, k_fmt="# {k:<25} ", unwrap="Item")
        print(dict_str, file=OUTFILE)

        # More details?
        run_steps = self.get_run_steps(as_name=True)
        run_str = utils.report_list_str(run_steps, pref="# Run_plan ", unwrap="step")
        geno_str = genes.report_genome_specs_str(self)
        kit_sc_str = kits.report_kit_score_str(self, unset=unset)
        plate_str = plates.report_plate_samples_str(self)

        print("#", file=OUTFILE)
        print(
            "#                ------ >>> Run details <<< ------------------",
            file=OUTFILE,
        )
        print("#", file=OUTFILE)
        print(run_str, file=OUTFILE)
        if geno_str:
            print("#", file=OUTFILE)
            print(geno_str, file=OUTFILE)
        if kit_sc_str:
            print("#", file=OUTFILE)
            print(kit_sc_str, file=OUTFILE)
        if plate_str:
            print("#", file=OUTFILE)
            print(plate_str, file=OUTFILE)

        print("#", file=OUTFILE)
        print("#" * 80, file=OUTFILE)
        print("\n", file=OUTFILE)
        # Flush output and close if new file
        OUTFILE.flush()
        if new_open:
            utils.close_out_fh(OUTFILE)

    def get_bc_info(self):
        """Get barcode information struct"""
        return self._bc_info

    def get_bc_df(self, bc_round=1):
        """Get barcode dataframe for given round

        Return dataframe
        """
        bc_info = self.get_bc_info()
        return bc_info.get_bc_data(bc_round)

    def get_bc_bci_list(self, bc_round=1):
        """Get list of barcode indexes for bc round"""
        df = self.get_bc_df(bc_round)
        return bcutils.get_bc_bci_list(df)

    def get_bc_wind_list(self, bc_round=1, unique=True):
        """Get list of well indexes for bc round"""
        df = self.get_bc_df(bc_round)
        return bcutils.get_bc_wind_list(df, unique=unique)

    def get_bc_well_list(self, bc_round=1, plate=True, unique=True):
        """Get list of (single) well strings (e.g. A4,D10) for bc round
        return str
        """
        df = self.get_bc_df(bc_round)
        return bcutils.get_bc_well_list(df, plate=plate, unique=unique)

    def get_bc_num_wells(self, bc_round=1):
        """Get number of wells for bc round"""
        df = self.get_bc_df(bc_round)
        return bcutils.get_bc_num_wells(df)

    def get_bc_rows_cols(self, from_kit=False, bc_round=1, as_list=False):
        """Get rows and cols for wells

        Return tuple of counts or lists
        """
        if from_kit:
            kit, chem = self.get_kit_chem()
            n_rows, _, _ = kits.kit_num_rows_cols_plates(kit, chem)
            rows, cols = plates.rows_cols_for_96plate(n_rows=n_rows, as_list=as_list)
        else:
            df = self.get_bc_df(bc_round)
            rows, cols = bcutils.get_bc_rows_cols(df, as_list=as_list)
        return rows, cols

    def get_wind_to_sample_dict(self, as_name=True, meta=False):
        """Get dict mapping well indexes to sample

        as_name     = Flag to return sample names rather than objects
        meta        = flag to consider meta samples; Normally no

        Return dict
        """
        # Init to all-well then replace with real samples
        aw_samp = self.get_allwell_sample(as_name=False)
        if as_name:
            # Default for non-covered wells = "all_well"
            aw_samp = "all_well" if aw_samp is None else aw_samp.get_name()

        # Dict for one or two plates; Based on barcode plates, not samples
        winds1 = self.get_bc_wind_list(bc_round=1)
        if self.num_sample_plates() > 1:
            winds2 = self.get_bc_wind_list(bc_round=2)
            new_dict = {
                r1 + "_" + r2: aw_samp for r1, r2 in list(product(winds1, winds2))
            }
        else:
            new_dict = {w: aw_samp for w in winds1}

        for samp in self.get_samples(meta=meta):
            # Ignore all-well except if there are no others
            if samp.is_allwell() and (self.num_samples() > 1):
                continue
            name = samp.get_name()
            for w in samp.get_wind_list(zero_pad=True):
                if as_name:
                    new_dict[w] = name
                else:
                    new_dict[w] = samp
        return new_dict

    def get_bc_num_bci(self, bc_round):
        """Get number of barcodes for bc round"""
        return len(self.get_bc_bci_list(bc_round))

    def get_gene_info(self, key=None):
        """Get genome information structs"""
        gene_info = self._gene_info
        if gene_info and key:
            gene_info = gene_info.get(key)
        return gene_info

    def get_guide_info(self, key=None):
        """Get fb/crispr guide information structs"""
        guide_info = self._guide_info
        if guide_info and key:
            guide_info = guide_info.get(key)
        return guide_info

    def get_targen_info(self, key=None):
        """Get any target enrichment struct"""
        targen = self._targen_info
        if targen and key:
            if key in targen:
                targen = targen[key]
        return targen

    def is_focal_crispr(
        self, focal=False, crispr=False, map_star=False, map_direct=False
    ):
        """Return true if focal barcoding / CRISPR logic

        focal = flag to check focal barcoding specifically
        crispr = flag to check crispr specifically
        map_star = flag to check if using STAR for guide mapping
        map_direct = flag to check if using direct guide mapping

        If no flags, then either use case qualifies

        Return bool
        """
        answer = False

        if (not focal) and (not crispr):
            focal = True
            crispr = True

        # Either current use case or parameters
        if focal:
            if self.get_use_case("focal_bc"):
                answer = True
            if self.get_par_val("focal_barcoding", as_bool=True):
                answer = True
        if crispr:
            if self.get_use_case("crispr"):
                answer = True
            if self.get_par_val("crsp_analysis", as_bool=True):
                answer = True

        # Check guide mapping?
        if answer and map_star:
            if not self.get_par_val("crsp_use_star", as_bool=True):
                answer = False
        if answer and map_direct:
            if self.get_par_val("crsp_use_star", as_bool=True):
                answer = False

        return answer

    def is_tcr(self, parent=False, only=False):
        """Return true if TCR (only or with parent)

        parent = flag to check tcr_parent use case specifically
        only = flag to check tcr_only use case specifically

        If no flags, then either use case qualifies

        Return bool
        """
        answer = False

        if (not parent) and (not only):
            parent = True
            only = True

        is_tcr = self.get_par_val("tcr_analysis", as_bool=True)
        have_parent = self.have_parent()
        if is_tcr:
            if parent and have_parent:
                answer = True
            if only and (not have_parent):
                answer = True
        return answer

    def is_targeted(self, par_only=True):
        """Return true if target enrichment"""
        is_targ = False
        if self.get_targen_info():
            is_targ = True
        if par_only and self.get_par_val("targeted_list", as_str=True):
            is_targ = True
        return is_targ

    def get_genome_list(self):
        """Get list of genomes"""
        genome_list = []
        gene_info = self.get_gene_info()
        guide_info = self.get_guide_info()
        if gene_info:
            genome_list = gene_info["genome_list"]
        elif guide_info:
            genome_list = guide_info["genome_list"]
        return genome_list

    def num_genomes(self):
        """Number of genomes"""
        return len(self.get_genome_list())

    def get_gene_list(self):
        """Get list of genes"""
        gene_list = []
        gene_info = self.get_gene_info()
        guide_info = self.get_guide_info()
        if gene_info:
            all_genes = gene_info["all_genes"]
        elif guide_info:
            all_genes = guide_info["all_genes"]
        gene_list = list(all_genes["gene_name"])
        return gene_list

    def num_genes(self):
        """Number of genes"""
        return len(self.get_gene_list())

    def get_exec(self, key):
        """Get executable path via name key"""
        path = ""
        try:
            xinfo = self._execs[key.lower()]
            if xinfo.status:
                path = xinfo.path
        except Exception:
            pass
        return path

    def check_read_csv(
        self,
        fkey,
        samp=None,
        sep=",",
        top_dir="",
        genome_dir="",
        names=None,
        index_col=0,
        verb=True,
    ):
        """Wrapper for potentially non-existant files"""
        df = None
        fname = self.filepath(fkey, samp, top_dir=top_dir, mkdir=False)
        if utils.check_infile(fname, verb=False):
            df = self.read_csv(
                fkey,
                samp=samp,
                sep=sep,
                top_dir=top_dir,
                genome_dir=genome_dir,
                names=names,
                index_col=index_col,
                verb=verb,
            )
        return df

    def read_csv(
        self,
        fkey,
        samp=None,
        sep=",",
        top_dir="",
        genome_dir="",
        names=None,
        index_col=0,
        verb=True,
    ):
        """Wrapper around df read function"""
        fname = self.filepath(fkey, samp, top_dir=top_dir, genome_dir=genome_dir)
        if names:
            if isinstance(names, str):
                names = names.split(",")
        df = utils.read_csv(fname, sep=sep, index_col=index_col, names=names)
        # df = utils.read_csv(fname)
        if verb:
            if df is None:
                self.report_run_story2(f"Failed to load {fname}")
            else:
                story = utils.report_df_info_str(df)
                self.report_run_story2(f"Read {story} from {fname}")
        return df

    def write_df(
        self,
        df,
        fkey,
        samp=None,
        top_dir="",
        genome_dir="",
        sep=",",
        index=True,
        sort_index=True,
        cols=None,
        mode="w",
        what="",
        quoting=False,
        verb=True,
        dump=True,
        na_rep="NaN",
        n_rows=5,
        row_num=True,
    ):
        """Wrapper around df write function (story goes to log)

        df = dataframe
        fkey = filename key
        samp = samp (in case filename is sample-specific)
        top_dir, genome_dir = also for filename key resolution

        Returns str
        """
        fname = self.filepath(fkey, samp, top_dir=top_dir, genome_dir=genome_dir)

        # If given a 'what' and not append, write header
        if what and (mode != "a"):
            ver = self.get_version()
            utils.file_header(ofile=fname, story=what, ver=ver)
            # append mode for the rest
            mode = "a"
            header = True
        else:
            header = False if mode == "a" else True

        # Make sure dataframe (in case series)
        df = utils.pd_series_to_df(df, cols=cols)
        if index and sort_index:
            df = df.sort_index(key=natsort_keygen())

        # Lower level write wrapper
        story = utils.write_df(
            df,
            fname,
            sep=sep,
            header=header,
            cols=cols,
            index=index,
            mode=mode,
            quoting=quoting,
            na_rep=na_rep,
            verb=verb,
        )
        self.report_run_story(story)

        # Report part of dataframe for output/log
        if dump:
            df_story = utils.report_df_str(
                df, pref=f"#2 {fkey} ", n_rows=n_rows, row_num=row_num
            )
            self.report_run_story(df_story, pref="")

        # Return what we got from writing
        return story

    def write_mtx(self, mtx, path, verb=True):
        """Wrapper around mtx write function (story goes to log)"""
        # Rows = cell, Cols = genes
        story = utils.write_mtx(mtx, path, verb=verb)
        self.report_run_story(story)
        return story

    def get_tscp_cutoff(self, tscp_counts, samp, targeted=False):
        """Get tscp count threshold to call cells

        tscp_counts = dataframe transcript count per barcode

        Return dict with tscp cutoff info
        """
        ok = False
        # If given value, use that, else calculate
        thresh = self.get_par_val("dge_tscp_cut_given", as_int=True)
        if thresh:
            self.report_run_story2(f"Given {thresh} tscp cutoff")
            recipe = {"thresh_given": thresh}
            ok = True
            results = {}
        else:
            # Recipe = parameters for cutoff determination
            # Defaults
            recipe = {
                "min_cnt": 30,
                "use_cnt": 0,
                "min_cells": 10,
                "max_cells": 30000,
                "min_cell_auc": 0.4,
                "use_cells": 0,
                "cv_steps": 1000,
                "cv_smooth_win": 0.02,
                "cv_slope_win": 0.02,
                "edge_pen_x": 0.5,
                "edge_pen_y": 0.25,
                "cnt_scale_fac": 1,
            }
            # Reset as needed
            h_min_c = self.get_par_val("dge_cell_hard_min", as_int=True)
            min_c = self.get_par_val("dge_cell_min", as_int=True)
            max_c = self.get_par_val("dge_cell_max", as_int=True)
            min_cell_auc = self.get_par_val("dge_cell_auc_min", as_float=True)
            use_c = self.get_par_val("dge_cell_given", as_int=True)
            if targeted:
                min_t = self.get_par_val("dge_tscp_cut_targeted_min", as_int=True)
            else:
                min_t = self.get_par_val("dge_tscp_cut_min", as_int=True)
            max_t = self.get_par_val("dge_tscp_cut_max", as_int=True)
            use_t = self.get_par_val("dge_tscp_cut_given", as_int=True)
            scale_fac = self.get_par_val("dge_tscp_cut_scale_fac", as_float=True)

            # Given or calc
            if use_t:
                story = f"Getting tscp cutoff; Given tscp {use_t}"
            elif use_c:
                story = f"Getting tscp cutoff; Given cells {use_c}"
            else:
                story = f"Calculating tscp cutoff; Cell min={min_c} ({h_min_c}), max={max_c}, Tscp min={min_t}"
            self.report_run_story2(story)

            # Scale cell limits by sample fraction?
            if self.get_par_val("dge_cell_scale_by_frac", as_bool=True):
                samp_wells = samp.num_wells()
                tot_wells = self.get_bc_num_wells(1)
                frac = samp_wells / tot_wells
                # Given cells or free calc
                if use_c:
                    use_c = round(use_c * frac)
                    story = f"Given cells scaled by {samp_wells} wells; Target cells={use_c}"
                else:
                    min_c = round(min_c * frac)
                    # Hard lower limit on cell min
                    min_c = max(min_c, h_min_c)
                    max_c = round(max_c * frac)
                    story = f"Cell bounds scaled by {samp_wells} wells; min={min_c}, max={max_c}"

                self.report_run_story2(story)

            recipe["use_cnt"] = use_t
            recipe["min_cnt"] = min_t
            recipe["use_cells"] = use_c
            recipe["min_cells"] = min_c
            recipe["max_cells"] = max_c
            recipe["min_cell_auc"] = min_cell_auc
            recipe["cnt_scale_fac"] = scale_fac

            # Get calculated cutoff
            # First check if possible; i.e. have at least min cells with min tscp?
            tscp_counts = tscp_counts.sort_values("tscp_count", ascending=False)
            tcount = int(tscp_counts.iloc[recipe["min_cells"] - 1].values)
            if tcount < recipe["min_cnt"]:
                mt = recipe["min_cnt"]
                mc = recipe["min_cells"]
                story = f"Cannot get cell cutoff; Min tscp {mt}, cells rank[{mc}] = {tcount} tscp"
                self.report_run_story(story)
                # Copy recipe into results
                results = {k: v for k, v in recipe.items()}
                results["problem"] = story
            else:
                # Calculate with array (series) of counts; Returns dict with results
                ok, results = ccutils.calc_count_cutoff(
                    tscp_counts["tscp_count"], recipe, plot=False, verb=True
                )

            if not ok:
                self.report_run_story("Failed cacl_count_cutoff")
                self.report_run_story2(results["problem"])
                thresh = 0
            else:
                thresh = results["threshold"]
                # Bound if calc and have limits (not for given tscp or cells)
                if use_t or use_c:
                    self.report_run_story2(
                        f"Raw tscp cutoff {thresh} NOT bounded (given tscp / cells)"
                    )
                else:
                    if min_t and (thresh < min_t):
                        self.report_run_story2(
                            f"Raw tscp cutoff {thresh} bounded to min {min_t}"
                        )
                        thresh = min_t
                    if max_t and (thresh > max_t):
                        self.report_run_story2(
                            f"Raw tscp cutoff {thresh} bounded to min {max_t}"
                        )
                        thresh = max_t

        # Counts are int; Also for consistant printing
        thresh = round(thresh)

        # Save input and output to file; Dump so sdtout and log too
        recipe["cell_call"] = "get_tscp_cutoff"
        recipe.update(results)
        recipe["cell_tscp_cutoff"] = thresh

        self.report_run_story2(f"Transcript cutoff status {ok}, thresh {thresh}")
        return ok, recipe

    def get_tscp_assign_df(
        self,
        top_dir="",
        sub_ind=-1,
        fresh=False,
        samp=None,
        targeted=False,
        focal=False,
    ):
        """Get trancript assignment dataframe (from file or memory)

        top_dir = top-level (output) path; If not given, use current 'output_dir'
        sub_ind = sublibrary index; If >= 0, append num to cell_barcode names
        fresh = flag to force fresh file load
        samp = sample to possibly filter by
        targeted = flag to filter by target enrichment genes
        focal = flag to filter by focal min counts

        Return dataframe
        """
        # If want fresh or don't already have, need to load
        if fresh or (not isinstance(self._tscp_assign, pd.DataFrame)):
            # Only need name, not path
            fname = self.filename("PF_TSCP_ASSIGN", None)
            # Get process dir path
            if not top_dir:
                top_dir = self.filepath("DIR_TOP", None)
            proc_dir = self.filepath("DIR_PROC", None, top_dir=top_dir)

            # Get full input file path with/without .gz extension, then load
            read_csv_path = utils.check_file_gz(
                fname, path=proc_dir, verb=True, toxic=False
            )
            self.report_run_story2(
                f"Loading transcript assignments from {read_csv_path}"
            )
            tas_df = utils.load_tscp_assign_df(
                fname, proc_dir, sub_ind=sub_ind, verb=True
            )
            story = utils.report_df_info_str(tas_df, name="tas_df")
            self.report_run_story2(story)
            # If multi-plate, set bc2_wind
            if self.num_sample_plates() > 1:
                tas_df["bc2_wind"] = tas_df["bc_wells"].apply(lambda s: s.split("_")[1])
            # Save for possible next time
            self._tscp_assign = tas_df

        tas_df = self._tscp_assign
        self.report_run_story2(f"Starting tscp records {len(tas_df)}")

        # Filtering by sample wells?
        if samp:
            if samp.is_allwell():
                self.report_run_story2("'all-well' sample so no masking")
            else:
                wspec = samp.get_wspec()
                self.report_run_story2(f"Masking tscp to wells ({wspec})")

                wind1_mask = set(samp.get_wind_list(plate=1, zero_pad=True))
                # Multi plate?
                if samp.is_multi_plate():
                    wind2_mask = set(samp.get_wind_list(plate=2, zero_pad=True))
                    tas_df = tas_df[
                        tas_df["bc1_wind"].isin(wind1_mask)
                        & tas_df["bc2_wind"].isin(wind2_mask)
                    ]
                else:
                    tas_df = tas_df[tas_df["bc1_wind"].isin(wind1_mask)]
                self.report_run_story2(f"After sample-mask tscp records {len(tas_df)}")

        # Masking target enrichment?
        if targeted:
            targen_info = self.get_targen_info()
            if targen_info:
                gene_id_set = set(targen_info["gene_id_list"])
                self.report_run_story2(
                    f"Masking tscp to target genes (n={len(gene_id_set)})"
                )
                tas_df = tas_df[tas_df["gene"].isin(gene_id_set)]
                self.report_run_story2(
                    f"After target-gene-mask tscp records {len(tas_df)}"
                )

        # Transcript min read count?
        if focal:
            min_count = self.get_par_val("dge_tscp_read_fb_umin", as_int=True)
        else:
            min_count = self.get_par_val("dge_tscp_read_min", as_int=True)
        if min_count > 1:
            self.report_run_story2(f"Masking read count per tscp >= {min_count}")
            tas_df = tas_df[tas_df["count"] >= min_count]
            self.report_run_story2(f"After min-read-count tscp records {len(tas_df)}")

        # Gene min tscp count?
        if focal:
            min_count = self.get_par_val("dge_gene_tscp_fb_umin", as_int=True)
        else:
            min_count = self.get_par_val("dge_gene_tscp_min", as_int=True)
        if min_count > 1:
            self.report_run_story2(f"Masking tscp count per gene >= {min_count}")
            # Filter with temp combo col of cell and gene
            tas_df["bc_wells_gene"] = tas_df["bc_wells"] + "___" + tas_df["gene"]
            mask = tas_df.groupby("bc_wells_gene").size()
            mask_set = set(mask[mask >= min_count].index)
            tas_df = tas_df[tas_df["bc_wells_gene"].isin(mask_set)]
            self.report_run_story2(
                f"After min-tscp-gene-count tscp records {len(tas_df)}"
            )
            tas_df = tas_df.drop("bc_wells_gene", axis=1)

        return tas_df

    def free_tscp_assign(self):
        """Clear memory of transcript assignment dataframe"""
        self._tscp_assign = None

    def get_tscp_counts(
        self, top_dir="", sub_ind=-1, samp=None, targeted=False, focal=False
    ):
        """Get sorted high-to-low transcript count per cell (barcode)

        Return sorted dataframe of transcript counts per barcode
        """
        tas_df = self.get_tscp_assign_df(
            top_dir=top_dir, sub_ind=sub_ind, samp=samp, targeted=targeted, focal=focal
        )
        tscp_per_cell = tas_df.groupby("bc_wells").size()
        # Dataframe, sorted by decending counts, ascending index
        df = pd.DataFrame(tscp_per_cell, columns=["tscp_count"])
        return df.sort_values(by=["tscp_count", "bc_wells"], ascending=[False, True])

    def have_in_fastq(self, fq1=True, fq2=True):
        """Check if fastq files are set and good

        Return status; Any missing = False
        """
        ok = True
        fq_list = []
        if fq1:
            fq_list.append(self.get_par_val("fq1", as_path=True))
        if fq2:
            fq_list.append(self.get_par_val("fq2", as_path=True))
        if fq_list:
            if not utils.check_infile(fq_list):
                ok = False
        return ok

    def get_fastq_samp_df(self, which, fresh=False, keep=True, split_bc=True):
        """Get fastq sample data

        which = flag to indicate which fastq; 0 = bc-correcte, 1 = R1, 2 = R2
        fresh = flag to force loading from file; Else may used prior data
        keep = flag to keep dataframe
        split_bc = flag to split barcodes (for R2 only)

        If too few reads in file, return False, None

        Return tuple (status, dataframe)
        """
        targ_df = None
        # Which fastq data? Parse regex True so can have letters (e.g. 'R1' >--> 1)
        # In case barcode corrected without number
        if str(which).lower()[0] == "b":
            which = 0
        ok, which = utils.str_to_num(which, regex=True)
        if ok:
            if which == 0:
                targ_df = self._fq0_samp
                fname = self.filepath("PF_FASTQ_BC", None)
                what = "BC-corrected"
            elif which == 1:
                targ_df = self._fq1_samp
                fname = self.get_par_val("fq1", as_path=True)
                what = "R1"
            elif which == 2:
                targ_df = self._fq2_samp
                fname = self.get_par_val("fq2", as_path=True)
                what = "R2"
            else:
                ok = False
        if not ok:
            story = f"get_fastq_samp_df 'which' not recognized; '{which}'"
            raise ValueError(story)

        # If want fresh or don't already have, need to load
        if fresh or (not isinstance(targ_df, pd.DataFrame)):
            first, last, range_str = self.get_par_slice("fastq_samp_slice", lower=True)
            num = last - first
            self.report_run_story(f"Sampling {num} {what} reads ({range_str})")
            # Attempt to load num reads
            heads, seqs, quals = utils.load_fastq_recs(fname, num, skip=first)
            if len(heads) < num:
                self.set_problem(
                    f"Failed to read {num} reads from {fname}; Got {len(heads)}"
                )
                ok = False
                targ_df = None
            else:
                # Make dataframe and possibly save
                targ_df = pd.DataFrame({"head": heads, "seq": seqs, "qual": quals})
                if keep:
                    if which == 0:
                        self._fq0_samp = targ_df
                    elif which == 1:
                        self._fq1_samp = targ_df
                    else:
                        self._fq2_samp = targ_df

        # For fq2, may want to split out barcode info
        if ok and split_bc and (which == 2):
            amp_seq = self.get_par_val("bc_amp_seq", as_str=True)
            targ_df = bcutils.fq2_df_split_bc(targ_df, amp_seq)

        return (ok, targ_df)

    def free_fastq_samp_df(self, bc=True, r1=True, r2=True):
        """Clear memory for fastq sample dataframes"""
        if bc:
            self._fq0_samp = None
        if r1:
            self._fq1_samp = None
        if r2:
            self._fq2_samp = None

    def get_vbc_count_df(self, fresh=False):
        """Get valid barcode counts as dataframe"""
        # Need to load?
        if fresh or (not isinstance(self._vbc_count, pd.DataFrame)):
            bc_df = self.read_csv("PF_VAL_BC_DATA")
            # Make cols with bc round and well; Index = <bc_round>_<well>
            bc_df["bc_round"] = bc_df.index.str.split("_").str[0].astype(int)
            bc_df["well"] = bc_df.index.str.split("_").str[1]
            self._vbc_count = bc_df
        return self._vbc_count

    def get_2vbc_count_df(self, fresh=False):
        """Get 2 valid barcode counts as dataframe"""
        # Need to load?
        if fresh or (not isinstance(self._2vbc_count, pd.DataFrame)):
            bc_df = self.read_csv("PF_VAL_2BC_DATA")
            # Data loaded as well_well count; Make 'well' col (e.g. 'D12_H8' 58)
            bc_df["well"] = bc_df.index
            self._2vbc_count = bc_df
        return self._2vbc_count

    def get_samp_vbc_count(self, samp, as_frac=False):
        """Get valid barcode count / fraction for given sample

        Return int/float
        """
        # One or 2-plate samples?
        if self.num_sample_plates() > 1:
            df = self.get_2vbc_count_df()
        else:
            df = self.get_vbc_count_df()
            # Sample-specific = round 1
            df = df[df["bc_round"] == 1]

        # Subset of wells or all if no sample given
        if isinstance(samp, spsamp.SpSamp):
            # plate 0 = all plates, combinations for multi-plate (e.g. 'A1_C4')
            well_list = samp.get_well_list(plate=0)
        else:
            well_list = df["well"]

        samp_n = df[df["well"].isin(well_list)]["count"].sum()
        # Fraction of total?
        if as_frac:
            total = df["count"].sum()
            samp_n = samp_n / total
        return samp_n

    def get_seq_pipe_stats_dict(self, top_dir="", verb=True):
        """Get dict with pipeline and sequencing stats data (from file)

        Return dict
        """
        df = self.read_csv("PF_STAT_PIPE", top_dir=top_dir, verb=verb)
        stat_dict = df["value"].to_dict()
        df = self.read_csv("PF_STAT_SEQ", top_dir=top_dir, verb=verb)
        stat_dict.update(df["value"].to_dict())
        # Save raw values as such
        if "number_of_reads" in stat_dict:
            stat_dict["raw_number_of_reads"] = stat_dict["number_of_reads"]
        if "number_of_tscp" in stat_dict:
            stat_dict["raw_number_of_tscp"] = stat_dict["number_of_tscp"]
        return stat_dict

    def get_stats_data(
        self, samp, top_dir="", fresh=False, asum=False, as_dict=False, targeted=False
    ):
        """Get stats data as dict or df, full or analysis_summary subset

        top_dir = top-level (output) path; If not given, use current 'output_dir'
        fresh = flag to force fresh file load
        asum = flag to only return analysis summary subset
        as_df = flag to return dataframe vs dict
        targeted = flag for targeted case

        Return dict or dataframe
        """
        # If given sample is different than remembered name, need to reload
        if samp.get_name() != self._stat_df_sample:
            fresh = True
        # Same for top_dir
        if top_dir != self._stat_df_top_dir:
            fresh = True

        # If want fresh or don't already have, need to load
        if fresh or (not isinstance(self._stat_df, pd.DataFrame)):
            # Only need name, not path
            self._stat_df = self.read_csv("SFR_ALLSTATS", samp=samp, top_dir=top_dir)

        self._stat_df_sample = samp.get_name()
        self._stat_df_top_dir = top_dir
        # Full dataframe or subset
        stat_df = self._stat_df

        clean_index = []
        for k in stat_df.index:
            k = k.replace("_mreads_", "_reads_")
            k = k.replace("_mread_", "_read_")
            clean_index.append(k)
        stat_df.index = clean_index

        stat_df.index.name = "statistic"

        if asum:
            # Returns list of keys to keep and list of discarded keys
            keep, toss = const.get_asum_out_keys(stat_df.index, targeted=targeted)
            stat_df = stat_df.loc[keep, :]

        # Return dict?
        if as_dict:
            # Stats file; For dict, 'final' (combined, whatever) in first col
            stat_df = stat_df.iloc[:, 0].to_dict()

        return stat_df
