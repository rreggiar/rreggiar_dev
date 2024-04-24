#!/usr/bin/env bash
#   Installs required libraries to get split-pipe running on AWS Ubuntu
#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

# 
#   This script requires sudo privledges
#

set -o nounset
set -o errexit


VERSION="Version 0.7; RTK 2024-01-05"

# For version, save program name from first var, first pass
PROGNAME=$0
version() { echo "$PROGNAME $VERSION"; }

LOG_FILE="install_dependencies.log"
TEMP_LOG=$(mktemp)

# Log from other install script
OTHER_LOG_FILE="install_dependencies_conda.log"


# Package handling commands; Install / purge
# For Ubuntu, apt; Otherwise? Maybe yum?
PKG_INSTALL="apt"

# Calls
UPDATE_CALL="sudo $PKG_INSTALL update"
CLEAN_CALL="sudo $PKG_INSTALL autoremove --yes"

INSTALL_CALL="sudo $PKG_INSTALL install --yes"
REMOVE_CALL="sudo $PKG_INSTALL purge --yes"
UPGRADE_CALL="sudo $PKG_INSTALL upgrade --yes"


# Package list
packlis="build-essential
        samtools
        rna-star
        unzip
        pigz
        zlib1g-dev
        libbz2-dev
        libcurl4-openssl-dev
        libssl-dev
        "


usage ()
{
    echo " "
    version
    echo " "
    echo "Use: [options]"
    echo "   -i --install   Install packages"
    echo "   -r --remove    Remove (purge) installed packages; i.e. remove"
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
    echo " "
}


# Turn off interactive pop up messages to restart deamons  
# i = interactive, a = automatic
APT_MSG_TARG_FILE="/etc/needrestart/needrestart.conf"
APT_MSG_TEMP_FILE="needrestart.conf.orig"

turn_off_apt_deamon_msg() {
    # No target, no action
    if [[ ! -f ${APT_MSG_TARG_FILE} ]]; then
        return
    fi

    echo " "
    echo "# Disabling deamon restart messages (in ${APT_MSG_TARG_FILE})"
    echo " "
    # Copy orig so can restore at end
    sudo cp $APT_MSG_TARG_FILE $APT_MSG_TEMP_FILE
    # Sub in temp files; Switch 'i' for 'a' May or may not start with '#' so two passes
    temp_file1=$(mktemp)
    temp_file2=$(mktemp)

    sed "/#\$nrconf{restart} = 'i';/s/.*/\$nrconf{restart} = 'a';/" "$APT_MSG_TEMP_FILE" > "$temp_file1"
    sed "/\$nrconf{restart} = 'i';/s/.*/\$nrconf{restart} = 'a';/" "$temp_file1" > "$temp_file2"

    sudo cp "$temp_file2" "$APT_MSG_TARG_FILE"
}

turn_on_apt_deamon_msg() {
    # No target, no action
    if [[ ! -f ${APT_MSG_TARG_FILE} ]]; then
        return
    fi

    echo " "
    echo "# Restoring deamon restart messages (in ${APT_MSG_TARG_FILE})"
    echo " "
    sudo cp "$APT_MSG_TEMP_FILE" "$APT_MSG_TARG_FILE"
    sudo rm "$APT_MSG_TEMP_FILE"
}


init_log() {
    version     > "$TEMP_LOG"
    date        >> "$TEMP_LOG"
    echo " "    >> "$TEMP_LOG"
    # If already have log, copy to backup (which also may already exist)
    if [[ -f $LOG_FILE ]]; then
        cp $LOG_FILE "${LOG_FILE}.old"
    fi
}


check_other_log() {
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

    # Init then parse; Only eat (initial) args starting with '-'
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
    if ! command -v $PKG_INSTALL &> /dev/null; then
        echo " "
        echo "Problem:"
        echo "Package install command ($PKG_INSTALL) not found; Giving up"
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
        packlis+="fastqc"
    else
        echo " "
        echo "Not considering fastQC"
        echo " "
    fi

    # Turn off interactive message settings
    if [[ $dry_run -eq 0 ]]; then
        turn_off_apt_deamon_msg
    fi

    # Init log file ($TEMP_LOG at this point)
    init_log

    # Update version list 
    if [[ $dry_run -eq 0 ]]; then
        echo " "                                                    | tee -a "$TEMP_LOG"
        echo "Updating versions (to install packages)"              | tee -a "$TEMP_LOG"
        $UPDATE_CALL                                                | tee -a "$TEMP_LOG"
        echo " "                                                    | tee -a "$TEMP_LOG"
        echo "Calling autoremove (to clean package collection)"     | tee -a "$TEMP_LOG"
        $CLEAN_CALL                                                 | tee -a "$TEMP_LOG"
        echo " "                                                    | tee -a "$TEMP_LOG"
    else
        echo "# $UPDATE_CALL           (dry run; not calling)"
    fi

    # Each listed package
    for pack in $packlis; do
        if [[ $do_remove -gt 0 ]]; then
            call="$REMOVE_CALL $pack"
        elif [[ $do_upgrade -gt 0 ]]; then
            call="$UPGRADE_CALL $pack"
        elif [[ $do_install -gt 0 ]]; then
            call="$INSTALL_CALL $pack"
        fi

        # Do or just report
        if [[ $dry_run -eq 0 ]]; then
            echo " "                                                        | tee -a "$TEMP_LOG"
            echo "#=======================================================" | tee -a "$TEMP_LOG"
            echo "# $call"                                                  | tee -a "$TEMP_LOG"
            echo "#=======================================================" | tee -a "$TEMP_LOG"
            $call                                                           | tee -a "$TEMP_LOG"
        else
            echo "# $call           (dry run; not calling)"
        fi
    done

    # If remove, do final clean up
    if [[ $do_remove -gt 0 ]]; then
        if [[ $dry_run -eq 0 ]]; then
            echo "Calling autoremove (to clean up)"                     | tee -a "$TEMP_LOG"
            $CLEAN_CALL                                                 | tee -a "$TEMP_LOG"
            echo " "                                                    | tee -a "$TEMP_LOG"
        else
            echo "# $CLEAN_CALL         (dry run; not calling)"
        fi
    fi

    # Restore interactive message settings
    if [[ $dry_run -eq 0 ]]; then
        turn_on_apt_deamon_msg
    fi
    
    # Save temp log to real thing
    cp "$TEMP_LOG" $LOG_FILE
    echo "# Log file: $LOG_FILE"

}


# Call main with all command line args
main "$@"

