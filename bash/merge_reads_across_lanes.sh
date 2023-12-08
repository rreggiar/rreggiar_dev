#!/bin/bash
# rreggiar
# 2023-12-07
# see: https://www.biostars.org/p/317385/#317471
# purpose: combines lane-split output from illumina in read-pairs

case $1 in
"")
	echo "ERROR: missing expected output arg in position 1"
	;;
*)
	outdir=$1 && mkdir -p $1
	;;
esac

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort -u); do

	echo "merging R1"

	cat "$i"_L00*_R1_001.fastq.gz >${outdir}/"$i"_L001_R1_001.fastq.gz

	echo "merging R2"

	cat "$i"_L00*_R2_001.fastq.gz >${outdir}/"$i"_L001_R2_001.fastq.gz

done
