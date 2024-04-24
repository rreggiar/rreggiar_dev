#!/usr/bin/env bash
#   Installs things required for TCR analysis with split-pipe
#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

set -o nounset
set -o errexit


VERSION="Version 0.4; RTK,AK 2023-11-18"

# For version, save program name from first var, first pass
PROGNAME=$0
version() { echo "$PROGNAME $VERSION"; }

LOG_FILE="install_TCR.log"
TEMP_LOG=$(mktemp)


TCR_TOP="./splitpipe/TCR/"
IGBLAST_DIR="./splitpipe/TCR/IgBlast"
IGBLAST_DBS="./splitpipe/TCR/IgBlast/tr_database/"
IGBLAST_LOG="./splitpipe/TCR/install_IgBlast.log"
IGB_INSTALL_SCRIPT="./splitpipe/scripts/install_IgBlast.sh"

OK_QUERY="Finished successfully"


usage ()
{
    echo " "
    version
    echo " "
    echo "Use: [options]"
    echo "   -i --install   Install TCR packages (IgBlast)"
    echo "   -I --imgt      Install TCR pagagees to use IMGT database"
    echo "   -r --remove    Remove TCR packages"
    echo "   -c --check     Check if TCR is installed"
    echo "   -y --yes       Actually execute commands; Default is just dryrun"
    echo "   -h --help      This help story"
    echo "   -v --version   Version"
    echo " "
    echo "Default behavior is dryrun;"
    echo "Need to specify --yes to actually do anything"
    echo " "
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
    do_check=0
    do_latest=0
    parse_prob=''
    db_type='default'
    while [[ $# -gt 0 ]] && [[ ${1:0:1} == - ]]; do
        # Simple set / call
        [[ $1 =~ ^-h|--help$ ]] && { usage; shift; exit 0; };
        [[ $1 =~ ^-v|--version$ ]] && { version; shift; exit 0; };
        [[ $1 =~ ^-y|--yes$ ]] && { dry_run=0; shift; continue; };
        [[ $1 =~ ^-i|--install$ ]] && { do_install=1; shift; continue; };
        [[ $1 =~ ^-I|--imgt$ ]] && { do_install=1; db_type='imgt'; shift; continue; };
        [[ $1 =~ ^-r|--remove$ ]] && { do_remove=1; shift; continue; };
        [[ $1 =~ ^-c|--check$ ]] && { do_check=1; shift; continue; };
        [[ $1 =~ ^-l|--latest$ ]] && { do_latest=1; shift; continue; };
        # What's this?
        parse_prob="Unknown option: $1"
        break
    done

    # Bail on any parse problem
    [[ -n $parse_prob ]] && { echo "Argument parse problem: $parse_prob"; exit 1; };

    # Need some chosen operation
    if [[ $do_install -eq 0 ]] && [[ $do_remove -eq 0 ]] && [[ $do_check -eq 0 ]]; then
        echo " "
        echo "Problem:"
        echo "Need to specify operation; --install or --remove or --check"
        exit 1
    fi

    # Track any setup modification
    mod_setup=0
    any_error=0

    # Init log file ("$TEMP_LOG" at this point)
    init_log

    echo " "

    # Install
    if [[ $do_install -gt 0 ]]; then
        # Create output if not there
        if [[ ! -d $TCR_TOP ]]; then
            if [[ $dry_run -eq 0 ]]; then
                mkdir -p $TCR_TOP
                echo "# Created TCR dir: $TCR_TOP"              | tee -a "$TEMP_LOG"
            else
                echo "# Dry run; Not creating TCR dir: $TCR_TOP"
            fi
        fi
        # Check if IgBlast dir already exists
        if [[ -d $IGBLAST_DIR ]]; then
            echo " "
            echo "Already have: $IGBLAST_DIR"
            echo "Remove this if you wish a fresh install"
            echo " "
            update_db=1
            ex_args='stable'
            echo "Updating igblast database if needed "
            if [[ $dry_run -eq 0 ]]; then
                all_good=0
                echo "# Calling $IGB_INSTALL_SCRIPT $TCR_TOP $ex_args $db_type $update_db" 
                /usr/bin/env bash $IGB_INSTALL_SCRIPT $TCR_TOP $ex_args $db_type $update_db     | tee -a "$TEMP_LOG"
                if [[ -f $IGBLAST_LOG ]]; then
                    all_good=$( grep -c "$OK_QUERY" $IGBLAST_LOG || true; )
                fi

                if [[ $all_good -gt 0 ]]; then
                    mod_setup=1
                else
                    echo "Problem installing IgBlast"               | tee -a "$TEMP_LOG"
                    any_error=1
                fi
            else
                echo "# Dry run; Not calling $IGB_INSTALL_SCRIPT $TCR_TOP $ex_args $db_type $update_db"
            fi
            exit 1
        fi

        # Extra arg to target latest version
        ex_args='stable'
        if [[ $do_latest -gt 0 ]]; then
            ex_args='latest'
        fi
        # Call the actual install script with <output_top>
        if [[ $dry_run -eq 0 ]]; then
            all_good=0
            echo "# Calling $IGB_INSTALL_SCRIPT $TCR_TOP $ex_args $db_type" 
            /usr/bin/env bash $IGB_INSTALL_SCRIPT $TCR_TOP $ex_args $db_type     | tee -a "$TEMP_LOG"
            if [[ -f $IGBLAST_LOG ]]; then
                all_good=$( grep -c "$OK_QUERY" $IGBLAST_LOG || true; )
            fi

            if [[ $all_good -gt 0 ]]; then
                mod_setup=1
            else
                echo "Problem installing IgBlast"               | tee -a "$TEMP_LOG"
                any_error=1
            fi
        else
            echo "# Dry run; Not calling $IGB_INSTALL_SCRIPT $TCR_TOP $ex_args $db_type"
        fi

        # Save temp log to real thing
        if [[ $dry_run -eq 0 ]]; then
            cp "$TEMP_LOG" $LOG_FILE
            echo "# Log file: $LOG_FILE"
        fi
    
    # Remove
    elif [[ $do_remove -gt 0 ]]; then
        # If have output, remove it
        if [[ -d $TCR_TOP ]]; then
            if [[ $dry_run -eq 0 ]]; then
                rm -r $TCR_TOP
                echo "# Removed TCR dir: $TCR_TOP"
                mod_setup=1
            else
                echo "# Dry run; Not removing TCR dir: $TCR_TOP"
            fi
        else
            echo "# No TCR dir found: $TCR_TOP"
        fi

    # Check what's installed
    elif [[ $do_check -gt 0 ]]; then
        echo "# Checking TCR status"

        # Look for success message in log file
        all_good=0
        if [[ -f $IGBLAST_LOG ]]; then
            all_good=$( grep -c "$OK_QUERY" $IGBLAST_LOG || true; )
        fi

        # Count db files and report first-token names
        if [[ $all_good -gt 0 ]]; then
            echo "# TCR appears to be installed"
            
            temp_file=$(mktemp)
            find $IGBLAST_DBS -name "*.ndb" | tr "/" " " | awk '{print $NF}' > "$temp_file"
            n_files=$(wc "$temp_file" | awk '{print $1}')
            echo "# Found $n_files TCR database indexes; Can run TCR analysis"

            n_imgt=$(grep imgt_ "$temp_file" | wc | awk '{print $1}')
            if [[ $n_imgt -gt 0 ]]; then
                imgt_story="IMGT database files are present; Can use IMGT"
            else
                imgt_story="IMGT database files not found; Cannot use IMGT"
            fi
            echo "# $imgt_story"

            db_names=$(sed 's/\.ndb/ /g' "$temp_file" | awk '{print $1}' | sort -u)
            echo "# DB name(s):"
            echo "$db_names"
        else
            echo "# TCR not installed; Files not found"
        fi
    fi

    if [[ $any_error -gt 0 ]]; then
        echo " "
        echo "Problem encountered running $0"
    # Changed anything requiring install?
    elif [[ $mod_setup -gt 0 ]]; then
        echo " "
        echo "                                 IMPORTANT"
        echo "            Local pipeline package TCR status has been changed"
        echo "You need to   >>>  Reinstall the Pipeline Package  <<<   for TCR to be used:"
        echo " "
        echo "                        pip install . --no-cache-dir"
        echo " "
    fi

    echo " "
}


# Call main with all command line args
main "$@"


