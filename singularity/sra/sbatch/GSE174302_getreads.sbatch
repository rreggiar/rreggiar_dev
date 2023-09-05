#!/usr/bin/bash
#SBATCH --job-name=getreads_GSE174302
#SBATCH --output=getreads_GSE174302.%j.out
#SBATCH --error=getreads_GSE174302.%j.err
#SBATCH --time=10:00:00
#SBATCH -p normal
#SBATCH -c 16
#SBATCH --mem=32GB

SRA_PATH="$HOME/rreggiar_dev/singularity/sra"
OUTPUT_PATH="$GROUP_SCRATCH/users/rreggiar/"

singularity exec ~/rreggiar_dev/singularity/sra/getreads_mamba_AMD.sif \
  nextflow run $SRA_PATH/getreads/main.nf \
  --list $SRA_PATH/GSE174302_accession.txt \
  --outdir $OUTPUT_PATH/circRNA_compatible_rnaSeq_20230830/GSE174302 \
  --resume
