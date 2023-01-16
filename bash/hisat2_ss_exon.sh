#!/bin/bash
# a shell script to map 30 samples to arabidopsis genome
# map with hisat2_ss, index built by hisat2-built

jobdt=$(date "+%Y%m%d")
idx_ss_exon="/mnt/r/genomes/dna_index/index_hisat2_ss-exon"
data_path="/mnt/r/reads_drRNA/"
out_path="/mnt/r/output/hisat2_ss/sams/"
log_path="/mnt/r/logs/"
opt="_ss_ex"

log=${log_path}${jobdt}_hisat2mapping.txt

test -e $out_path && echo "output path $out_path exist." || mkdir $out_path
presample=""
count=0
echo "============== NEW JOB STARTED AT $(date "+%Y-%m-%d %H:%M:%S") ===========================" | tee -a $log

for rd in "$data_path"*R1*.fastq;
do
	base_name=$(echo $rd|cut -d"/" -f 5)
	sample=$(echo $base_name|cut -d"_" -f 1-2)    #LCM_4-3_S11_R1_001_FC4.fastq
	flow_cell=$(echo $base_name | cut -d"_" -f 6 | cut -d"." -f 1)

        if [[ $sample == $presample ]];
        then
             echo " ~~~ $rd has been processed. skipped in this round" | tee -a $log
             echo " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee -a $log
             echo | tee -a $log
             continue
        else
            presample=$sample
        fi

        m1reads=$(ls ${data_path}${sample}*R1*.fastq)
        m1=""
        m2=""
        for r in $m1reads;
        do
           m1=${m1},$r
           m2=${m2},${r/R1/R2}
        done
        m1=${m1:1}
        m2=${m2:1}

	out_file="$out_path${sample}$opt"
	new_ss="$out_path${sample}_summary.txt"
	log="${log_path}${dt}_log.txt"

	started=$(date "+%Y-%m-%d %H:%m%s")
	echo "processing... $rd [$started]" | tee -a $log
	echo "	+ m1 file: $m1"  | tee -a $log
	echo "	+ m2 file: $m2"  |  tee -a $log
	echo "	+ output file: $out_file" | tee -a $log

    if [[ $1 == wet ]];
    then
	   hisat2 -p 4 -x $idx_ss_exon -1 $m1 -2 $m2 -S $out_file --time --dta -I 50 #--new-summary $new_ss
          #--novel-splicesite-outfile $new_ss 
          #  --min-intron 20
    fi
    finished=$(date "+%Y-%m-%d %H:%M:%S")
    ((count=$count+1))
    echo "--------------- [ $finished ]-------------------"  | tee -a $log
    echo | tee -a $log
done

echo "<<<<<<<<< JOB FINISHED, TOTAL SAMPLES: $count >>>>>>>>>>>>>>>>>" | tee -a $log
echo | tee -a $log
#hisat2 -x /mnt/e/genomes/dna_index_hisat2/hisat2_index -1 LCM_1-1_S19_R1_001.fastq -2 LCM_1-2_S20_R2_001.fastq -S ../output/hisat2_default --min-intronlen 5 -I 60 --time



