#!/usr/bin/bash
#SBATCH --job-name=bwaIndex
#SBATCH --output=out.bwaIndex_%j
#SBATCH --error=err.bwaIndex_%j
#SBATCH --time=36:00:00
#SBATCH -p owners,sfgf,howchang
#SBATCH -c 16
#SBATCH --mem=128GB


#Usage:   bwa index [options] <in.fasta>

#Options: -a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]
#         -p STR    prefix of the index [same as fasta name]
#         -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
#         -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* 

#Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
#         `-a div' do not work not for long genomes.

module load biology bwa

REF_DIR="/home/groups/howchang/users/rreggiar/references"

mkdir "${REF_DIR}/bowtie"

cd "${REF_DIR}/bowtie"

bwa index \
 -b 500000000 \
 "${REF_DIR}/gencode/44/GRCh38.p14.genome.fa"
