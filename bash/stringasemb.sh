#!/bin/bash
#stringtie quantification test

# setup
jobdate=$(date "+%Y%m%d")
ref_gtf="../genomes/gtf/Arabidopsis_thaliana.TAIR10.59.gtf"
bfile_path="../output/hisat2_ss/bams/"
out_path="../output/hisat2_ss/transAssembly/"
log_path="../logs/"
log=${log_path}${jobdate}_stringtie_assembly.log

echo  ============ JOB STARTED $(date "+%Y-%m-%d %H:%M:%S") ========= | tee -a $log
count=0

# make a fold if the output fold is not exist
assemfoldexist=$(test -e ${out_path} && echo TRUE || echo FALSE)
if [[ $assemfoldexist == FALSE ]];
then
    echo making directory: ${out_path} | tee -a $log
    mkdir ${out_path}
fi

# start to loop through bam files
for bam in ${bfile_path}*.bam;
do
   #get basic file info
   base_name=$(echo $bam | cut -d"/" -f 5)
   label=$(echo $base_name | cut -d"_" -f 1-2)
   #FC=$(echo $base_name | cut -d"_" -f 4)
   #gene_adb=${out_path}${label}_abund.ctab
   out_gtf=${out_path}${label}.gtf
   full_cov=${out_path}${label}_fcr.gtf
   echo current file  $bam
   echo sample:       $label
   echo out put:      $out_gtf
   #skip the file if the corresponding gtf file already exist
   fileexist=$(test -e $out_gtf && echo TRUE || echo FALSE)
   if [[ $fileexist == TRUE ]];
   then
      echo "~~~ $bam has been processed. " | tee -a $log
      echo "~~~ skipping $bam" | tee -a $log
      echo  | tee -a $log
      continue
   else
      start=$(date "+%Y-%m-%d %H:%M:%S")
      echo "||| processing: $bam  [$start]" | tee -a $log
      echo "    output_gtf: $out_gtf" | tee -a $log
   fi

   if [[ $1 == wet ]];
   then
  #  stringtie $bam -eB -m 100 -G $ref_gff -A $gene_adb -o $gtf -p 4
     stringtie $bam -p 4 -m 75 -G $ref_gtf -C $full_cov -l $label -o $out_gtf #-A $gene_adb
   else
     echo " ^_^ I did nothing...... "
   fi
   ((count+=1))
   finished=$(date "+%H:%M:%S")

   echo ----------[ $start - $finished ]---------- | tee -a $log
   echo | tee -a $log
done

echo "<<<< JOB FINISHED, Total Files: $count >>>>" | tee -a $log
echo | tee -a $log
