#!/bin/bash
# qualimap quality control

inpath="../output/hisat2_drRNA/"

conf_bamqc="${inpath}conf_multi_bamqc.txt"
conf_bamqc_EWR="${inpath}conf_multi_bamqc_EWR.txt"
conf_bamqc_LCM="${inpath}conf_multi_bamqc_LCM.txt"

out_bamqc="${inpath}bamqc_ALL/"
out_bamqc_EWR="${inpath}bamqc_EWR/"
out_bamqc_LCM="${inpath}bamqc_LCM/"

echo "multiple bamqc...  ALL SAMPLES"
echo "out put folder: $out_bamqc"
qualimap multi-bamqc -d $conf_bamqc -outdir $out_bamqc

echo "out put folder: $out_bamqc_EWR"
echo "multiple bamqc...  EWR SAMPLES"
qualimap multi-bamqc -d $conf_bamqc_EWR -outdir $out_bamqc_EWR

echo "out put folder: $out_bamqc_LCM"
echo "multiple bamqc...  LCM SAMPLES"
qualimap multi-bamqc -d $conf_bamqc_LCM -outdir $out_bamqc_LCM

conf_counts="${inpath}conf_counts.txt"
conf_counts_EWR="${inpath}conf_counts_EWR.txt"
conf_counts_LCM="${inpath}conf_counts_LCM.txt"

out_counts="${inpath}countsQC_ALL/"
out_counts_EWR="${inpath}countsQC_EWR/"
out_counts_LCM="${inpath}countsQC_LCM/"


qualimap counts -c -d $conf_counts -outdir $out_counts
qualimap counts -c -d $conf_counts_EWR -outdir $out_counts_EWR
qualimap counts -c -d $conf_counts_LCM -outdir $out_counts_LCM
