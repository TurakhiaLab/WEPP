#!/bin/bash -ex
# Script for aligning amplicon reads to SARS-CoV2 reference
aligner="$1"
reads="$2"
cores=`grep -c ^processor /proc/cpuinfo`

# Use minimap2 as aligner (default)
if [[ "${aligner}" == "minimap2" ]]; then
  minimap2 -ax sr "wuhCor1.fa" "${reads}" > "data/amplicons.sam"
fi

# Use bwamem2 as alternate aligner (experimental)
if [[ "${aligner}" == "bwamem2" ]]; then
  ./bwamem2/bwa-mem2 index "wuhCor1.fa" 
	# Mapping (defaults to number of cores)
  ./bwamem2/bwa-mem2 mem -t "${cores}" "${reads}" > "data/amplicons.sam"
fi



