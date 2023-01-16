#!bash/bin

ref_gtf="../genomes/gtf/Arabidopsis_thaliana.TAIR10.59.gtf"
outpath="../output/hisat2_ss/gffcompare/"
merged_gtf="../output/hisat2_ss/transAssembly/stringtie_merged.gtf"
merged_EWR_gtf="../output/hisat2_ss/transAssembly/stringtie_merged_EWR.gtf"
merged_LCM_gtf="../output/hisat2_ss/transAssembly/stringtie_merged_LCM.gtf"

test -e $outpath && echo "outpath exists" || mkdir $outpath

#gffcompare -r $ref_gtf -G -o ${outpath}merged $merged_gtf
#gffcompare -r $ref_gtf -G -o ${outpath}merged_EWR $merged_EWR_gtf
#gffcompare -r $ref_gtf -G -o ${outpath}merged_LCM $merged_LCM_gtf
gffcompare -r $merged_gtf -G -o ${outpath}EWRvsMerged $merged_EWR_gtf
gffcompare -r $merged_gtf -G -o ${outpath}LCMvsMerged $merged_LCM_gtf
gffcompare -r $merged_EWR_gtf -G -o ${outpath}LCMvsEWR $merged_LCM_gtf
