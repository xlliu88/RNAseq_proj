#!/bin/bash

jobdt=$(date "+%Y%m%d")

ref_gtf="../genomes/gtf/Arabidopsis_thaliana.TAIR10.59.gtf"
path="../output/hisat2_ss/"
inpath="${path}transAssembly/"
outpath="${path}transAssembly/"
log="../logs/${jobdt}_gtfmerge_log.txt"

mergelist="${path}mergelist_ALL.txt"
mergelist_EWR="${path}mergelist_EWR.txt"
mergelist_LCM="${path}mergelist_LCM.txt"

#commands used in this script
#comm1=echo $(stringtie --merge -p 4 -G $ref_gtf -o ${outpath}stringtie_merged.gtf $mergelist)
#comm2=echo $(stringtie --merge -p 4 -G $ref_gtf -o ${outpath}stringtie_merged.gtf $mergelist_EWR)
#comm3=echo $(stringtie --merge -p 4 -G $ref_gtf -o ${outpath}stringtie_merged.gtf $mergelist_LCM)

dt=$(date "+%Y-%m-%d %H:%M:%S")
echo "========== JOB STARTED AT $dt =============" | tee -a $log

# check if outpath exist. if not, make the directory
outdirexist=$(test -e $outpath && echo TRUE || FALSE)
if [[ $outdirexist == FALSE ]];
then
    echo "making output directory. $outpath"  | tee -a $log
    mkdir $outpath
    echo   | tee -a $log
fi

dt=$(date "+%Y-%m-%d %H:%M:%S")
echo "merging gtfs -- ALL SAMPLES [$dt]"  | tee -a $log
echo   | tee -a $log
stringtie --merge -p 4 -G $ref_gtf -o ${outpath}stringtie_merged.gtf $mergelist

dt=$(date "+%Y-%m-%d %H:%M:%S")
echo "merging gtfs -- EWR SAMPLES [$dt]"   | tee -a $log
echo  | tee -a $log
stringtie --merge -p 4 -G $ref_gtf -o ${outpath}stringtie_merged_EWR.gtf $mergelist_EWR

dt=$(date "+%Y-%m-%d %H:%M:%S")
echo "merging gtfs -- LCM SAMPLES [$dt]"   | tee -a $log
echo
stringtie --merge -p 4 -G $ref_gtf -o ${outpath}stringtie_merged_LCM.gtf $mergelist_LCM

dt=$(date "+%Y-%m-%d %H:%M:%S")
echo "<<<<< JOB FINISHED [$dt] >>>>>>"  | tee -a $log
echo   | tee -a $log
