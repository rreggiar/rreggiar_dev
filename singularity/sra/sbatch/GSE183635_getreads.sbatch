#!/usr/bin/bash
#SBATCH --job-name=getreads_GSE183635
#SBATCH --output=getreads_GSE183635.%j.out
#SBATCH --error=getreads_GSE183635.%j.err
#SBATCH --time=48:00:00
#SBATCH -p owners
#SBATCH -c 16
#SBATCH --mem=32GB

SRA_PATH="$HOME/rreggiar_dev/singularity/sra"
OUTPUT_PATH="$GROUP_SCRATCH/users/rreggiar/"

singularity exec ~/rreggiar_dev/singularity/sra/getreads_mamba_AMD.sif \
  nextflow run $SRA_PATH/getreads/main.nf \
  --list $SRA_PATH/data/GSE183635_accession.txt \
  --outdir $OUTPUT_PATH/circRNA_compatible_rnaSeq_20230830/GSE183635 \
  --resume
