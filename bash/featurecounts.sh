#!bin/bash
# use feature to count reads.
# use Arabidopsis_thaliana.TAIR10.59.gtf as reference gtf
# all samples counted together.

jobdt=$(date "+%Y%m%d")
bampath="../output/hisat2fc/bams/"
gtf="../genomes/gtf/Arabidopsis_thaliana.TAIR10.59.gtf"
saf="../genomes/gtf/Arabidopsis_thaliana.TAIR10.59.saf"

outpath="../output/hisat2fc/featureCounts/"
counts_EWR=${outpath}counts_EWR.txt
counts_LCM=${outpath}counts_LCM.txt
counts_ALL=${outpath}counts_all.txt

log="../logs/${jobdt}_featureCounts_log.txt"

EWR_in=""
LCM_in=""
ALL_in=""
Ecount=0
Lcount=0
total=0

for bam in ${bampath}*.bam;
do
    if [[ $bam == *EWR*.bam ]];
    then
        echo "adding $bam to EWR samples"
        EWR_in="${EWR_in} ${bam}"
        #echo $infiles
        ((Ecount++))
    elif [[ $bam == *LCM*.bam ]];
    then
        echo "adding $bam to LCM samples"
        LCM_in="${LCM_in} ${bam}"
        ((Lcount++))
    fi
    ALL_in="${ALL_in} ${bam}"
    ((total++))
done
EWR_in=${EWR_in:1}
LCM_in=${LCM_in:1}
ALL_in=${ALL_in:1}

echo "$Ecount bam files for EWR samples"
echo "$Lcount bam files for LCM samples"
echo "$total bam files for ALL samples"
echo
echo "   " $EWR_in
echo "   " $LCM_in

test -e $outpath && echo output directory exists: $outpath || mkdir $outpath

# it will count single-end mapped reads
# will not count chimeric fragments
# multiple mapped reads will be assigned to primary mapped site
#featureCounts -T 4 -p -C -M --primary -s 2 -g gene_id -t exon -a $gtf -o $counts_EWR $EWR_in
#featureCounts -T 4 -p -C -M --primary -s 2 -g gene_id -t exon -a $gtf -o $counts_LCM $LCM_in
featureCounts -T 4 -p -C -M --primary -s 2 -g gene_id -t exon -a $gtf -o $counts_ALL $ALL_in
  # -T threads
  # -p count fragments instead of reads
  # -P check validaty of pair-end distance default: -d 50 -D 600
  # -C don't count pairs map to different chr, or different strands of same chr
  # -B count pairs both mapped
  # -M count multiple mapped pairs 
  # --primary
  # -s strandness. 0 unstranded, 1 forward specific, 2 reverse specific

## ref: gtf
## feature: gene -t gene#
#featureCounts -T 4 -p -C -M --primary -s 2 -g gene_id -t gene -a $gtf -o $counts_ALL_gene $ALL_in

## ref: saf
## 
#featureCounts -T 4 -p -C -M --primary -s 2 -g gene_id -F SAF -a $saf -o $counts_ALL_saf $ALL_in
