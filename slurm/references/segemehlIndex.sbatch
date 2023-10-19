#!/usr/bin/bash
#SBATCH --job-name=segemehlIndex
#SBATCH --output=out.segemehlIndex_%j
#SBATCH --error=err.segemehlIndex_%j
#SBATCH --time=24:00:00
#SBATCH -p owners
#SBATCH -c 32
#SBATCH --mem=128GB

new_REF_DIR="$OAK/users/rreggiar/references"

REF_DIR="/home/groups/howchang/users/rreggiar/references"

segemehlPath="$new_REF_DIR/segemehl" 

module load contribs poldrack anaconda/2009.03

conda activate rreggiar

mkdir -p ${segemehlPath}

segemehl.x \
    -t 32 \
    -d "${REF_DIR}"/gencode/44/GRCh38.p14.genome.fa \
    -x "${segemehlPath}"/GRCh38.p14.genome.idx
