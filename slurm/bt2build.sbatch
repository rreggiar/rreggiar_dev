#!/usr/bin/bash
#SBATCH --job-name=bt2Build
#SBATCH --output=out.bt2Build_%j
#SBATCH --error=err.bt2Build_%j
#SBATCH --time=36:00:00
#SBATCH -p owners,sfgf,howchang
#SBATCH -c 16
#SBATCH --mem=128GB

module load biology bowtie2

REF_DIR="/home/groups/howchang/users/rreggiar/references"

mkdir "${REF_DIR}/bowtie2"

bowtie2-build \
 --threads 16 \
 "${REF_DIR}/gencode/44/GRCh38.p14.genome.fa" \
 "${REF_DIR}/bowtie2/"
