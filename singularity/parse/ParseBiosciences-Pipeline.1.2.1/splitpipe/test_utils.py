# Some utils to help testing (e.g. use comargs in notebooks)
#


def cli_to_dict(cli_str, add_speedup=True):
    """Wrapped command line string into dict with opt[vals]

    Given string expected to be like command line:

        --mode comb \
        --rseed 123 \
        --sublib_pref /newvolume/analysis/AllKit-v3/ \
        --sublib_suff _v121a \
        --output_dir /newvolume/analysis/AllKit-v3/div_S6_comb12_v121a \
        --sublibraries div_S6_part1 div_S6_part2

    Return dict[opt] = val
    """
    cli_args = {}
    for raw_line in cli_str.split("--"):
        # Strip everything after hash or amperstand
        raw_line = raw_line.split('#')[0]
        raw_line = raw_line.split('&')[0]
        raw_line = raw_line.replace(">", "")
        # Any extra spaces, new lines
        line = raw_line.strip()
        if not line:
            continue
        # If more than one arg, keep as list
        parts = line.split()
        if len(parts) < 2:
            continue
        elif len(parts) > 2:
            key = parts[0]
            val = parts[1:]
        else:
            key, val = parts
        cli_args[key] = val

    # Add more args?
    if add_speedup:
        speed_dict = speedup_args_dict()
        cli_args.update(speed_dict)

    return cli_args
    

def speedup_args_dict():
    """Return dict with parameters to speed up pipeline init/setup steps

    Return dict
    """
    new_dict = {}
    # Speed up hacks
    new_dict["log_memory_use"] = False
    new_dict["start_timeout"] = 0
    new_dict["kit_score_skip"] = False
    new_dict["fastq_samp_slice"] = [100000, 200000]
    return new_dict


def comb_output_set_pref(output_dir, vertag=""):
    """From given full output path, set comb-mode prefix / suff

    Return dict
    """
    new_dict = {}
    new_dict["sublib_pref"] = "/".join(output_dir.split("/"))[:-1]
    new_dict["sublib_suff"] = vertag
    return new_dict
