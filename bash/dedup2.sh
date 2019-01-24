#!/bin/bash

#Remove duplicates
#code from Helen Lindsay

picard='java -jar /home/Shared_penticton/data/seq/mirco_mets_mbdseq/src2/picard/build/libs/picard.jar'

for f in ../bam/*_s.bam
do
  base=$(basename $f)
  out="../bam_dedup/"${base}
  metrics="../metrics_new/"${base}"_dedup.txt"
  $picard MarkDuplicates INPUT=$f OUTPUT=${out} METRICS_FILE=${metrics} REMOVE_DUPLICATES=true
done
