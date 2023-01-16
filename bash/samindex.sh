#!/bin/bash
# a shell script to index sam files
# input: /output/hisat2ss/bams
# output bai files
# 11/20/2020
# Xunliang Liu

jobdt=$(date "+%Y%m%d")
outpath="../output/hisat2ss/"
bampath=${outpath}bams/
sum_file="../output/summary.txt"
logpath="../logs/"

log=$logpath${jobdt}_samindex.log

for f in "$bampath"*.bam;
do
    echo indexing $f | tee -a $log
    samtools index -b $f
done
