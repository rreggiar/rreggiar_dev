#!/usr/bin/bash
#SBATCH --job-name=getreads_PRJNA789110
#SBATCH --output=getreads_PRJNA789110.%j.out
#SBATCH --error=getreads_PRJNA789110.%j.err
#SBATCH --time=24:00:00
#SBATCH -p owners,howchang
#SBATCH -c 1
#SBATCH --mem=5GB

SRA_PATH="$HOME/rreggiar_dev/singularity/sra"
OUTPUT_PATH="$GROUP_SCRATCH/users/rreggiar/"

module load java/18.0.2 biology nextflow

#singularity exec ~/rreggiar_dev/singularity/sra/getreads_mamba_AMD.sif \
nextflow run $SRA_PATH/getreads/main.nf \
 -profile slurm \
 -with-singularity "/scratch/groups/howchang/users/rreggiar/nf_singularity_cache/getreads.img" \
 --list $SRA_PATH/data/PRJNA789110_accession.txt \
 --outdir $OUTPUT_PATH/circRNA_compatible_rnaSeq_20230830/PRJNA789110
