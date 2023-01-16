#!/bin/bash
# qualimap quality control

inpath="../output/hisat2_drRNA/"

ingff="../genomes/gffa3/Arabidopsis_thaliana.TAIR10.59.gff3"
ingtf="../genomes/gtf/Arabidopsis_thaliana.TAIR10.59.gtf"

conf="${inpath}exp_config.txt"
outpath_multi_bamqc="${inpath}multi-bamqc/"

c=0

for bf in ${inpath}*.bam;
do
    echo "processing $bf"

	base_name=$(echo $bf | cut -d"/" -f 4 | cut -d"_" -f 1-4)
    outpath_bamqc=${inpath}${base_name}_bamqc/
    geno_cove=${outpath_bamqc}genome_coverage.txt
	outpath_Rsqc=${inpath}${base_name}_seqqc/
    counts=${outpath_Rsqc}counts.txt

	fileexist=$(test -e $counts && echo TRUE || echo FALSE)
        echo "testing if file exists..."

	if [[ $fileexist == TRUE ]];
    then 
        echo "~~~~ file $bf has been processed: (file status: $fileexist)"
	    echo "~~~~ skipping .... ~~~"
	    continue
        fi

    if [[ $1 == wet ]];
    then
        echo "    bam-QC for:           $base_name"
        echo "    bamqc results in:     $outpath_bamqc"
        echo "    genome coverage file: $geno_cove"
            qualimap bamqc -bam $bf -c -gff $ingtf -ip -nr 1000 -nt 4 -nw 400 -oc $geno_cove -outdir $outpath_bamqc --java-mem-size=6G

        echo "    RNA-seq-QC for:       $base_name"
        echo "    bamqc results in:     $outpath_bamqc"
            echo "    counting file:        $counts"
            qualimap rnaseq -bam $bf -gtf $ingtf -oc $counts -outdir $outpath_Rsqc -p non-strand-specific --paired --sorted
            ((c+=1))
    else
        echo "dry run, nothing done"
	fi

    echo "---------------------------------------"
    #break
done

echo "<<<<<<<<<<< DONE  TOTAL FILES PROCESSED: $c >>>>>>>>>>>>"




