#!/bin/bash

#Debarcode FASTQ files
#code from Helen Lindsay

b1=TGACCAAT
b2=CTTGTAAT
b3=GCCAATAT

snm=${1%%.*}
zcat $1 | awk -v b1=$b1 -v b2=$b2 -v b3=$b3 -v snm=$snm '(NR % 4 == 1){if ($1 ~ b1)
  print substr($1,2,length($1)) >> snm"_"b1".txt";
else if ($1 ~ b2) print substr($1,2,length($1)) >> snm"_"b2".txt";
else if ($1 ~ b3) print substr($1,2,length($1)) >> snm"_"b3".txt";
}'

snm=${1%%.*}
for nms in $b1 $b2 $b3
do
  ( seqtk subseq $1 ${snm}_${nms}.txt > ${snm}_${nms}.fastq && rm ${snm}_${nms}.txt
  gzip ${snm}_${nms}.fastq )
done
