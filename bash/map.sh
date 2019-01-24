#!/bin/bash

#Align debarcoded FASTQ files with bwa mem
#Code from Helen Lindsay

annot_dir=/home/Shared_sherborne/data/annotation/Human/hg19/
idx=${annot_dir}bwa_mem/hg19.fa

while read filename sample condition patient
do 
  base=../bam/${sample}
  r1="../FASTQ/"${filename}".fastq.gz"
  bwa mem -t 16 $idx $r1 | samtools view -Sb - > ${base}.bam
  samtools sort ${base}.bam -@ 10 -o ${base}_s.bam 
  samtools index ${base}_s.bam
done < <(awk '(NR > 1)' ../metadata.txt)
