#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
import shutil

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import pb_const as const
    import utils
else:
    from splitpipe import pb_const as const
    from splitpipe import utils


# ---------------------------------------------------------------------------
def run_clean(spipe):
    """Handle pipeline cleanup step"""
    spipe.report_proc_step("Clean")
    ok = True

    if spipe.get_par_val("clean_make_summaries_zip", as_bool=True):
        ok = make_summaries_zip(spipe)

    # Collect filenames to delete
    # Per sample
    f_lis = []

    # anndata (not TCR)
    if not spipe.is_tcr():
        for samp in spipe.get_samples():
            key_dict = const.out_filepaths(spipe, samp, as_key=True)
            fkey = key_dict["anndata"]
            # For focal / crispr, anndata treated like temp
            if spipe.is_focal_crispr():
                if not spipe.keep_temp_files():
                    f_lis.append(spipe.filepath(fkey, samp))
            else:
                if not spipe.get_par_val("ana_save_anndata", as_bool=True):
                    f_lis.append(spipe.filepath(fkey, samp))

    # Global
    f_lis.append(spipe.filepath("TMP_TEMP1_FILE", None))
    f_lis.append(spipe.filepath("TMP_TEMP2_FILE", None))
    f_lis.append(spipe.filepath("TMP_TEMP3_FILE", None))
    # Raw bam file and barcode fastq
    if not spipe.keep_temp_files():
        f_lis.append(spipe.filepath("PF_BAM_RAW", None))
        # f_lis.append(spipe.filepath("PF_FASTQ_BC", None))

    # Remove
    spipe.report_run_story2(f"Checking to remove {len(f_lis)} files")
    n = utils.check_and_rm_file(f_lis)
    if n < 0:
        spipe.set_problem(f"rm issue for files: {f_lis}")
        ok = False
    else:
        spipe.report_run_story2(f"  Removed {n}")

    # Compress
    gz_lis = []
    if ok and spipe.compress_temp_files():
        # List of files to possibly compress
        f_lis = []
        f_lis.append(spipe.filepath("PF_TSCP_ASSIGN", None))
        f_lis.append(spipe.filepath("PF_FASTQ_BC_RAW", None))
        # Collect files that exists and aren't already gzipped
        for f in f_lis:
            zf = utils.check_file_gz(f, path=None, verb=False, toxic=False)
            if zf:
                if not zf.endswith(".gz"):
                    gz_lis.append(zf)

    if len(gz_lis) > 0:
        spipe.report_run_story2(f"Compressing {len(f_lis)} files")
        comp_lev = spipe.get_par_val("gzip_compresslevel", as_int=True)
        ok = utils.check_and_gzip_file(gz_lis, comp_lev=comp_lev)
    else:
        spipe.report_run_story2("Not compressing temp files")

    spipe.report_proc_step("Clean", status=ok)
    return ok


def make_summaries_zip(spipe):
    """Make a zip file with all summary files

    Return status
    """
    if spipe.get_mode("mkref"):
        return True

    ok = True

    # Remove any existing zip and temp dir
    zip_file = spipe.filepath("SF_SUMMARIES_ZIP", None)
    utils.check_and_rm_file(zip_file)
    zip_dir = spipe.filepath("TMP_SUM_ZIP_DIR", None)
    utils.check_and_rm_file(zip_dir, isdir=True)

    # Collect all files; sample html and aggregate csv
    f_list = []
    for samp in spipe.get_samples():
        fname = spipe.filepath("SF_ASUM_HTML", samp)
        if utils.check_infile(fname, verb=False):
            f_list.append(fname)
    fname = spipe.filepath("SF_AGG_ASUM_CSV", None)
    if utils.check_infile(fname, verb=False):
        f_list.append(fname)
    # Add FastqQC htmls
    fpath = spipe.filepath("DIR_FASTQC", None, mkdir=False)
    f_list += utils.fnames_for_path(fpath, ext="html")

    # Add log (up to this point); Report current status first
    story = spipe.get_status_story()
    spipe.report_run_story(story)
    f_list.append(spipe.get_log_filename())
    # run_proc json
    fname = spipe.filepath("PF_RUNPROC_DEF", None, mkdir=False)
    f_list.append(fname)
    # Add STAR output log
    fname = spipe.filepath("PF_FASTQ_FIN_OUT", None, mkdir=False)
    if utils.check_infile(fname, verb=False):
        f_list.append(fname)

    # Generate zip (if anything more than log file)
    if len(f_list) > 1:
        # Create new temp dir and copy files
        utils.make_subdir(zip_dir)
        utils.copy_file(f_list, zip_dir)

        # make_archive adds '.zip' to the name, so strip that first
        zip_base = zip_file[:-4]
        shutil.make_archive(zip_base, "zip", zip_dir)

        # Check if we've got zip and report good, else problem
        if utils.check_infile(zip_file, verb=False):
            spipe.report_run_story(f"Created new zip archive: {zip_file}")
        else:
            spipe.set_problem(f"Failed to created zip archive: {zip_file}")
            ok = False

    # Remove temp zip dir
    utils.check_and_rm_file(zip_dir, isdir=True)

    return ok
