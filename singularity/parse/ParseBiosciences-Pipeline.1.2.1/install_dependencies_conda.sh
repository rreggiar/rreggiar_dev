#!/usr/bin/env bash
#   Installs required libraries to get split-pipe running on AWS Ubuntu
#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

#
#   This script does not require sudo privledges
#
#   Some specific conda version and channel settings seem quite brittle;
#       In other words, specific conda calls may hang / fail
#

set -o nounset
set -o errexit


VERSION="Version 1.0; RTK 2024-01-05"

# For version, save program name from first var, first pass
PROGNAME=$0
version() { echo "$PROGNAME $VERSION"; }

LOG_FILE="install_dependencies_conda.log"
TEMP_LOG=$(mktemp)

# Log from other install script
OTHER_LOG_FILE="install_dependencies.log"


# Package handling via conda
# Calls
INSTALL_CALL="conda install --yes"
REMOVE_CALL="conda remove --yes"
UPGRADE_CALL="conda upgrade --yes"


# Channel list; Order may matter; Highest priority first, lowest last
# Bioconda called explicitly, so can be last / lowest here
chanlis="
        defaults
        conda-forge
        bioconda
"

# Package list (Only programs/tools here; Python libraries via pip)
# Comma,sep,values. 
# Processed 'word' at time, so no spaces
#   older versions of samtools don't work
#   numba needs numpy<= 1.24
packlis="
        unzip 
        pigz
        -c,bioconda,star>=2.7
        -c,bioconda,samtools>=1.16      
        "


usage ()
{
    echo " "
    version
    echo " "
    echo "Use: [options]"
    echo "   -i --install   Install packages"
    echo "   -r --remove    Remove installed packages"
    echo "   -u --upgrade   Upgrade installed packages; i.e. latest version"
    echo "   -y --yes       Actually execute commands; Default is just dryrun"
    echo "   --fastQC       Add fastQC to package list; Otherwise skip"
    echo "   -h --help      This help story"
    echo "   -v --version   Version"
    echo " "
    echo "Install command: '${INSTALL_CALL}'"
    echo "Remove command:  '${REMOVE_CALL}'"
    echo "Upgrade command: '${UPGRADE_CALL}'"
    echo " "
    echo "Default behavior is dryrun;"
    echo "Need to specify --yes to actually do anything"
    echo " "
    echo "NOTE: conda operations can take a long time... Please be patient"
    echo " "
    echo " "
}


init_log() {
    # Initialize install log 
    version     > "$TEMP_LOG"
    date        >> "$TEMP_LOG"
    echo " "    >> "$TEMP_LOG"
    # If already have log, copy to backup (which also may already exist)
    if [[ -f $LOG_FILE ]]; then
        cp $LOG_FILE "${LOG_FILE}.old"
    fi
}


check_other_log() {
    # Check if 'other' install script has run and left a log
    if [[ -f $OTHER_LOG_FILE ]]; then
        echo " "
        echo "Problem"
        echo "Found other log file: $OTHER_LOG_FILE"
        echo " "
        echo "Only *one* install script should be run, not both"
        exit 1
    fi
}


main ()
{
    if [[ $# -lt 1 ]]; then
        usage
        exit 1
    fi

    # Init then parse party; Only eat (initial) args starting with '-'
    # Default is dry-run
    dry_run=1
    do_install=0
    do_remove=0
    do_upgrade=0
    do_fastqc=0
    parse_prob=''
    while [[ $# -gt 0 ]] && [[ ${1:0:1} == - ]]; do
        # Simple set / call
        [[ $1 =~ ^-h|--help$ ]] && { usage; shift; exit 0; };
        [[ $1 =~ ^-v|--version$ ]] && { version; shift; exit 0; };
        [[ $1 =~ ^-y|--yes$ ]] && { dry_run=0; shift; continue; };
        [[ $1 =~ ^-i|--install$ ]] && { do_install=1; shift; continue; };
        [[ $1 =~ ^-r|--remove$ ]] && { do_remove=1; shift; continue; };
        [[ $1 =~ ^-u|--upgrade$ ]] && { do_upgrade=1; shift; continue; };
        [[ $1 =~ ^--fast ]] && { do_fastqc=1; shift; continue; };
        # What's this?
        parse_prob="Unknown option: $1"
        break
    done

    # Bail on any parse problem
    [[ -n $parse_prob ]] && { echo "Argument parse problem: $parse_prob"; exit 1; };

    # Check if pkg install command exists
    if ! command -v conda &> /dev/null; then
        echo " "
        echo "Problem:"
        echo "Required 'conda' not found; Giving up"
        exit 1
    fi

    # Need some chosen operation
    if [[ $do_install -eq 0 ]] && [[ $do_remove -eq 0 ]] && [[ $do_upgrade -eq 0 ]]; then
        echo " "
        echo "Problem:"
        echo "Need to specify operation; --install, --remove, or --upgrade"
        exit 1
    fi

    # Check if other version log
    check_other_log

    # FastQC
    if [[ $do_fastqc -gt 0 ]]; then
        packlis+="-c,bioconda,fastqc>=0.12"
    else
        echo " "
        echo "Not considering fastQC"
        echo " "
    fi

    # Init log file ("$TEMP_LOG" at this point)
    init_log

    # Update channel list 
    echo "# Updating conda channels (for install packages)"
    for chan in $chanlis; do
        # call="conda config --add channels $chan"
        # Append puts channel at end of list (lowest priority); add puts at top
        call="conda config --append channels $chan"
        if [[ $dry_run -eq 0 ]]; then
            echo " "                                | tee -a "$TEMP_LOG"
            echo "# $call"                          | tee -a "$TEMP_LOG"
            $call                                   | tee -a "$TEMP_LOG"
        else
            echo "# Dry run; Not calling: $call"
        fi
    done

    # Turn off channel priority to get latest versions:
    #   https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html
    conda config --set channel_priority disabled

    # Collection of current packages
    list_file=$(mktemp)
    echo "# Getting list of currently installed packages"
    conda list > "$list_file"

    # Each listed package
    for pack in $packlis; do
        # Strip any comments
        pack=$(echo "$pack" | cut -d "#" -f 1 | awk '{print $1}')
        if [[ -z $pack ]]; then
            continue
        fi

        # Get package name alone; Last comma,sep,part before any version
        #   e.g. '-c,bioconda,samtools>=1.11'  >-->  'samtools'
        pkg_name=$(echo "$pack" | sed 's/[=<>]/ /g' | awk '{print $1}' | tr ',' ' ' | awk '{print $NF}')

        # Fancier grep so no match doesn't fail out
        have_pkg=$( (grep "$pkg_name" "$list_file" || true) | wc -l)
        # If don't have and remove, skip
        if [[ $have_pkg -lt 1 ]] && [[ $do_remove -gt 0 ]]; then
            echo "# Package $pkg_name is not currently installed via conda; skipping"
            continue
        fi

        # If remove, strip any version from call (i.e. after > or = , etc) 
        if [[ $do_remove -gt 0 ]]; then
            pack=$(echo "$pack" | sed 's/[=<>]/ /g' | awk '{print $1}')
        fi
        # Replace any commas with space
        pack=${pack//,/ }

        if [[ $do_remove -gt 0 ]]; then
            call="$REMOVE_CALL $pack"
        elif [[ $do_upgrade -gt 0 ]]; then
            call="$UPGRADE_CALL $pack"
        elif [[ $do_install -gt 0 ]]; then
            call="$INSTALL_CALL $pack"
        fi

        # Do or just report
        if [[ $dry_run -lt 1 ]]; then
            echo " "                                                        | tee -a "$TEMP_LOG"
            echo "#=======================================================" | tee -a "$TEMP_LOG"
            echo "# $call"                                                  | tee -a "$TEMP_LOG"
            echo "#=======================================================" | tee -a "$TEMP_LOG"
            $call                                                           | tee -a "$TEMP_LOG"
        else
            echo "# $call          (dry run; not calling)"
        fi
    done

    # Save temp log to real thing
    cp "$TEMP_LOG" $LOG_FILE
    echo "# Log file: $LOG_FILE"

}


# Call main with all command line args
main "$@"


