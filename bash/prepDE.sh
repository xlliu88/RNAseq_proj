#!bin/bash

jobdt=$(date "+%Y%d%m")
path="../output/hisat2ss/"
stpath=${path}stringtie/
ctpath=${path}stringcount/

# sample_list
EWR_sample_list=${path}prepde_sample_list_EWR-rf.txt
LCM_sample_list=${path}prepde_sample_list_LCM-rf.txt

# output files
EWR_gene_count=${ctpath}/EWR_gene_count.csv
EWR_trans_count=${ctpath}/EWR_transcript_count.csv

LCM_gene_count=${ctpath}/LCM_gene_count.csv
LCM_trans_count=${ctpath}/LCM_transcript_count.csv

# test if output directory exists
test -e ${ctpath} && echo "directory $depath exists" || mkdir $depath

# prep matrix for DEseq2
python prepDE.py -i $EWR_sample_list -g $EWR_gene_count -t $EWR_trans_count
python prepDE.py -i $LCM_sample_list -g $LCM_gene_count -t $LCM_trans_count
