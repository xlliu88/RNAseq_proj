#!/bin/bash
# a shell script to map 18 samples to arabidopsis genome
#map with hisat2_default, index built by hisat2-built

idx="/mnt/e/genomes/dna_index_hisat2/hisat2_index"
data_path="/mnt/e/reads/"
out_path="/mnt/e/output/"
sum_file="/mnt/e/output/hisat2/summary.txt"

m1_mark="R1"
m2_mark="R2"

for rd in "$data_path"*R1*;
do
    base_name=$(echo $rd|cut -d"/" -f 5)
	if [[ $base_name == *LCM_* ]];
	then
		echo "skipping LCM samples"
		continue
	fi
	
	sample=$(echo $base_name|cut -d"_" -f 1)
	flow_cell=$(echo $base_name|cut -d"_" -f 2)
	out_file="/mnt/e/output/hisat2/${sample}_${flow_cell}_default"
	mate="${rd/R1/$m2_mark}"
	
	echo "processing... $rd" | tee -a $sum_file
	echo "	m1 file: $rd" | tee -a $sum_file
	echo "	m2 file: $mate" | tee -a $sum_file
	echo "	output file: $out_file" | tee -a $sum_file
	
	hisat2 -x $idx -1 $rd -2 $mate -S $out_file \
		--min-intronlen 20 \
		--time	| tee -a $sum_file
		
	echo "-------------------------" | tee -a $sum_file
  #fi
done
#hisat2 -x /mnt/e/genomes/dna_index_hisat2/hisat2_index -1 LCM_1-1_S19_R1_001.fastq -2 LCM_1-2_S20_R2_001.fastq -S ../output/hisat2_default --min-intronlen 5 -I 60 --time



