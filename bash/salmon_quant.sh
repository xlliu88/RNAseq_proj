#!/bin/bash
# map RNAseq data to transcript_index

idx="../genomes/cdna_index2_salmon"
reads="../giantcell/reads/"
test_reads="../reads_drRNA/EWR_7-1_S1_R1_001_FC1_drRNA.fastq"
out="../giantcell/quant/"
out_exist=$(test -d $out && echo TRUE || echo FALSE)
if [[ $out_exist == FALSE ]];
then
    mkdir $out
fi

for fn in ${reads}DRR1870{33..40}.fastq;
do
  echo "processing...$fn"
  samp=`basename ${fn}`
  samp=$(echo $samp | cut -d"." -f 1)
  res=${out}${samp}_quant
  echo "result file: $res"
  salmon quant -i $idx -l U -r $fn -p 4 --validateMappings -o $res
done

# echo mapping testing reads
# test_out=${out}_test_quant
#salmon quant -i $idx -l U -r $test_reads --validateMappings -o $test_out
