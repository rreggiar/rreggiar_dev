#!/usr/bin/bash
#SBATCH --job-name=generateGenome
#SBATCH --output=out.generateGenome_%j
#SBATCH --error=err.generateGenome_%j
#SBATCH --time=36:00:00
#SBATCH -p owners,sfgf
#SBATCH -c 16
#SBATCH --mem=128GB

module load biology star/2.7.10b

REF_DIR="/home/groups/howchang/users/rreggiar/references"

mkdir "${REF_DIR}/star"

STAR \
 --runMode genomeGenerate \
 --genomeDir "${REF_DIR}/star/" \
 --genomeFastaFiles "${REF_DIR}/gencode/44/GRCh38.p14.genome.fa" \
 --sjdbGTFfile "${REF_DIR}/gencode/44/gencode.v44.annotation.gtf" \
 --runThreadN 16
