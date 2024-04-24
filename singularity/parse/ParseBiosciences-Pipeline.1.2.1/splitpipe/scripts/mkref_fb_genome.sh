#!/usr/bin/env bash
#   Runs mkref step for focal / CRISPR with --gfasta input (e.g. gRNA seqs)
#
#   Copyright (c) 2023 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

# set -o errexit
set -o nounset


VERSION="Version 0.2; RTK 2023-10-26"

# For version, save program name from first var, first pass
PROGNAME=$0
version() { echo "$PROGNAME $VERSION"; }


main ()
{
    if [[ $# -lt 3 ]]; then
        echo " "
        version
        echo " "
        echo "Usage: <genome_name> <guides> <output_dir> <pars>"
        echo " "
        echo "  <genome_name>   Genome name for splipe-pipe --genome_name"
        echo "  <guides>        Guide file for split-pipe --crsp_guides"
        echo "  <output_dir>    Path for split-pipe mkref --output_dir"
        echo "  <pars>          Path for split-pipe mkref --parfile"
        echo " "
        exit 1
    fi

    genome_name=$1
    guides=$2
    output_dir=$3
    parfile=$4
    
    call="split-pipe \
        --mode mkref \
        --crispr \
        --genome_name $genome_name \
        --crsp_guides $guides \
        --crsp_use_star \
        --output_dir $output_dir \
        --parfile $parfile \
        "
    
    echo " "
    version
    echo "Call: ${call}"

    # Call then get exit code
    eval "$call"

    exit_code=$?

    # Report non-zero exit
    if [[ $exit_code -ne 0 ]]; then
        echo "Problem with call: '${call}'" 
    fi
    echo "Exit code $exit_code" 

    # Explit exit code (Python 'subprocess.run' seems to need this(?))
    exit $exit_code
    
}


# Call main with all command line args
main "$@"

