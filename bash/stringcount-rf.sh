#!/bin/bash
#stringtie quantification test
jobdt=$(date "+%Y%m%d")

ref_gff="../genomes/gtf/Arabidopsis_thaliana.TAIR10.59.gtf"
merged_gtf="../output/hisat2ss/transAssembly/stringtie_merged.gtf"
#merged_gtf_EWR="../output/hisat2ss/transAssembly/stringtie_merged_EWR.gtf"
#merged_gtf_LCM="../output/hisat2ss/transAssembly/stringtie_merged_LCM.gtf"
bfile_path="../output/hisat2ss/bams/"
out_path="../output/hisat2ss/stringtie/"
log="../logs/${jobdt}_stringcount.log"

#outdirexist=$(test -e $out_path && echo TRUE || echo FALSE)
#if [[ $outdirexist == FALSE ]];
#then
#    echo "making output directory: " | tee -a $log
#    echo "   $out_path" | tee -a $log
#    mkdir $out_path
#fi
test -e $out_path && echo "$out_path Exists" || mkdir $out_path

echo "=========== JOB STARTED AT: $(date "+%Y-%m-%d %H:%M:%S") =========="  | tee -a $log

count=0
for bam in ${bfile_path}*7-1*sorted.bam;
do 
   start=$(date "+%Y-%m-%d %H:%M:%S")
   echo "||| processing: $bam" [$start] | tee -a $log

   base_name=$(echo $bam | cut -d"/" -f 5)    #LCM_1-1_S1_FC3_sorted.bam
   sample=$(echo $base_name | cut -d"_" -f 1-2)
   #FC=$(echo $base_name| cut -d"_" -f 4)
   label=${sample}
   sample_out=${out_path}${label}/
   gene_adb=${sample_out}${label}_abund.ctab
   out_gtf=${sample_out}${label}.gtf

   echo "  + checking if $base_name has been processed..."
   fileexist=$(test -e $gene_adb && echo TURE || echo FALSE)
   if [[ $fileexist == TURE ]]
   then
        echo "  file exist: $gene_adb" | tee -a $log
        echo "  ~~~ skipping: $base_name" | tee -a $log
        echo "  ~~~~~~~~~~~~~~~~" | tee -a $log
        fileexist=FALSE
        continue
   fi

   # if [[ $label == LCM* ]];
   # then
      # gtf=$merged_gtf_LCM
   # elif [[ $label == EWR* ]];
   # then
      # gtf=$merged_gtf_EWR
   # else
      # gtf=$ref_gff
   # fi

   mkdir $sample_out
   echo "||| stringtieing $bam" | tee -a $log
   echo "  + out directory:    $sample_out"| tee -a $log
   echo "  + out gene_abund:   $gene_adb"| tee -a $log
   echo "  + out gene_gtf:     $gtf"| tee -a $log

   if [[ $1 == wet ]];
   then
     #gtf=$ref_gff
     gtf=$merged_gtf
     stringtie $bam -e -B --rf -G $gtf -A $gene_adb -o $out_gtf -p 4
   fi

   ((count+=1))
   finish=$(date "+%Y-%m-%d %H:%M:%S")
   echo "--------------[ $finish ] -------------"   | tee -a $log
   echo  | tee -a $log
   #break
done

echo "<<<<< JOB DONE. Total Samples: $count >>>>>" | tee -a $log

