#!/usr/bin/bash
#SBATCH --job-name=getreads_GSE174302
#SBATCH --output=getreads_GSE174302.%j.out
#SBATCH --error=getreads_GSE174302.%j.err
#SBATCH --time=24:00:00
#SBATCH -p owners,howchang,sfgf
#SBATCH -c 1
#SBATCH --mem=5GB

SRA_PATH="$HOME/rreggiar_dev/singularity/sra"
FETCH_NGS="$GROUP_HOME/users/rreggiar/fetch-ngs/nf-core-fetchngs_1.10.0/1_10_0"
OUTPUT_PATH="$GROUP_SCRATCH/users/rreggiar/"

module load java/18.0.2 biology nextflow

#singularity exec ~/rreggiar_dev/singularity/sra/getreads_mamba_AMD.sif \
#nextflow run $SRA_PATH/getreads/main.nf \
# -resume \
# -profile singularity,slurm \
# -with-singularity "/scratch/groups/howchang/users/rreggiar/nf_singularity_cache/getreads.img" \
# --list $SRA_PATH/"data"/GSE174302_accession.txt \
# --outdir $OUTPUT_PATH/circRNA_compatible_rnaSeq_20230830/GSE174302

nextflow run $FETCH_NGS \
 -profile singularity,slurm \
 --input $SRA_PATH/"data"/GSE174302_accession.txt \
 --outdir $OUTPUT_PATH/circRNA_compatible_rnaSeq_20230830/GSE174302
