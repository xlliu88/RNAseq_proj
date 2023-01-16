#!/bin/env bash
#SBATCH -J build_index_hisat2
#SBATCH -o job_files/hisat_index.sbatch.o_%j
#SBATCH -e job_files/hisat_index.sbatch.e_%j
#SBATCH --partition BioCompute,Lewis,General
#SBATCH --mem 50G
#SBATCH --ntasks 2
#SBATCH --nodes 2
#SBATCH --time 03:59:00
#SBATCH --cpus-per-task 10

#hisat2-build
ref_path="/storage/hpc/data/liuxu/genomes/genomes/"
ref_genome="dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
exons="gtf/Arabidopsis_thaliana.TAIR10.59.exons.txt"
splice_sites="gtf/Arabidopsis_thaliana.TAIR10.59.ss.txt"

out="dna_index/index_hisat2"
out_ss_exon="dna_index/index_hisat2_ss-exon"

#echo ${ref_path}$ref_genome

#echo ${ref_path}$splice_sites
#echo ${ref_path}$exons 
#echo ${ref_path}$out_ss_exon

#build un-annotated index
srun hisat2-build -f ${ref_path}$ref_genome ${ref_path}$out
#build annotated index
srun hisat2-build --ss ${ref_path}$splice_sites --exon ${ref_path}$exons -f ${ref_path}$ref_genome ${ref_path}$out_ss_exon
