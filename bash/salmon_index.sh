#!/bin/bash
# to build index for salmon

cdna="../genomes/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
index="../genomes/cdna_index2_salmon_k10"
mkdir $index

salmon index -t $cdna -i $index -k 11

