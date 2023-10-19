#!/usr/bin/bash
#SBATCH --job-name=hisatBuild
#SBATCH --output=out.hisatBuild_%j
#SBATCH --error=err.hisatBuild_%j
#SBATCH --time=24:00:00
#SBATCH -p owners,bigmem,howchang
#SBATCH -c 32
#SBATCH --mem=256GB


#HISAT2 version 2.1.0 by Daehwan Kim (infphilo@gmail.com, http://www.ccb.jhu.edu/people/infphilo)
#Usage: hisat2-build [options]* <reference_in> <ht2_index_base>
#    reference_in            comma-separated list of files with ref sequences
#    hisat2_index_base       write ht2 data to files with this dir/basename
#Options:
#    -c                      reference sequences given on cmd line (as
#                            <reference_in>)
#    --large-index           force generated index to be 'large', even if ref
#                            has fewer than 4 billion nucleotides
#    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
#    -p                      number of threads
#    --bmax <int>            max bucket sz for blockwise suffix-array builder
#    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
#    --dcv <int>             diff-cover period for blockwise (default: 1024)
#    --nodc                  disable diff-cover (algorithm becomes quadratic)
#    -r/--noref              don't build .3/.4.ht2 (packed reference) portion
#    -3/--justref            just build .3/.4.ht2 (packed reference) portion
#    -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)
#    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
#    --localoffrate <int>    SA (local) is sampled every 2^offRate BWT chars (default: 3)
#    --localftabchars <int>  # of chars consumed in initial lookup in a local index (default: 6)
#    --snp <path>            SNP file name
#    --haplotype <path>      haplotype file name
#    --ss <path>             Splice site file name
#    --exon <path>           Exon file name
#    --seed <int>            seed for random number generator
#    -q/--quiet              verbose output (for debugging)
#    -h/--help               print detailed description of tool and its options
#    --usage                 print this usage message
#    --version               print version information and quit

new_REF_DIR="$OAK/users/rreggiar/references"

REF_DIR="/home/groups/howchang/users/rreggiar/references"

hisatPath="$new_REF_DIR/hisat2" 

module load biology hisat2

hisat2-build \
    -p 32 \
    --ss "${hisatPath}"/gencode.v44.annotation.splice_sites.txt \
    --exon "${hisatPath}"/gencode.v44.annotation.exons.txt \
    "${REF_DIR}"/gencode/44/GRCh38.p14.genome.fa \
    "${hisatPath}"/GRCh38.p14.genome
