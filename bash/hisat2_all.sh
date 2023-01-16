#!/bin/bash
# a shell script to map 18 samples to arabidopsis genome
#map with hisat2_default, index built by hisat2-built

idx="/mnt/r/genomes/dna_index_hisat2/hisat2_index"
data_path="/mnt/r/reads/"
out_path="/mnt/r/output/"
sum_file="/mnt/r/output/hisat2/summary.txt"

m1_mark="R1"
m2_mark="R2"
m1=""
m2=""
nm1=0
nm2=0

for rd in "$data_path"*R1*;
do
    base_name=$(echo $rd|cut -d"/" -f 5)
	# if [[ $base_name == *LCM_* ]];
	# then
		# echo "skipping LCM samples"
		# continue
	# fi
	
	sample=$(echo $base_name|cut -d"_" -f 1)
	flow_cell=$(echo $base_name|cut -d"_" -f 2)
	mate="${rd/R1/$m2_mark}"
	
	echo "adding... $rd to m1"
	m1="${m1},${rd}"
	nm1=nm1+1
	echo "adding $mate to m2"
	m2="${m2},${mate}"
	nm2=nm2+1
	echo "-------------------------"
done


out_file="/mnt/e/output/hisat2/hisat2_all"
m1=${m1:1}
m2=${m2:1}
echo "m1 files: ${nm1};"
echo "m2 files: ${nm2};"

echo "running hist2 mapping..."
start_time=(echo date "+%Y-%m-%d %H:%M:%S"）
echo "  Started at: ${start_time}"

#hisat2 -x $idx -1 $m1 -2 $m2 -S $out_file \
#		--min-intronlen 20 \
#		--time

finish_time=(echo date "+%Y-%m-%d %H:%M:%S"）
echo "  Finished at: ${finish_time}"




