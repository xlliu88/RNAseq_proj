#!/bin/bash
#samtools sort
#sort sam from hisat2 for stringtie

log_dt=$(date "+%Y%m%d")
dt=$(date "+%Y-%m-%d %H:%M:%S")

sfold="../output/hisat2fc/sams/"
bfold="../output/hisat2fc/bams/"
logfold="../logs/"
count=0
log=${logfold}${log_dt}_samsort_log.txt

test -e $bfold && echo "output fold exists: $$bfold" || mkdir $bfold

echo "=============JOB STARTED AT [$dt] ================" | tee -a $log
count=0

for sam in ${sfold}*_sx;
do
    started=$(date "+%Y-%m-%d %H:%M:%S")
	echo "processing $sam ...[ $started ]"  | tee -a $log

	base_name=$(echo $sam | cut -d"/" -f 5) #LCM_1-1_ss_ex
	n1=$(echo $base_name | cut -d"_" -f 1-3)
	sorted=$bfold${n1}_sorted.bam

	fileexist=$(test -e $sorted && echo TRUE || echo FALSE)

	if [[ $fileexist == TRUE ]];
	then
		echo "~~~ $sorted exists" | tee -a $log
		echo "~~~ skiping $sam" | tee -a $log
		continue
	fi

	echo "||| sorting: 	    $sam" | tee -a $log
	echo "    base name:    $base_name" | tee -a $log
	echo "    out_put:	    $sorted" | tee -a $log

	if [[ $1 == wet ]]
	then
	  samtools view -Su $sam | samtools sort -o $sorted
#	  samtools view -Su $sam | samtools sort $sorted
	fi

    count=$((count+1))
    finished=$(date "+%Y-%m-%d %H:%M:%S")
	echo "-------[ $finished] -----" | tee -a $log
	echo  | tee -a $log
#	break
done

echo "<<< -- JOB FINISHED - total files processed: $count -- >>>" | tee -a $log
echo "<<< -- sorted bam files saved at: $bfold -- >>>" | tee -a $log
echo
