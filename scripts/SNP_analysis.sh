#!/bin/bash
#$ -ckpt restart
#$ -q sam,sam128,abio,abio128,pub64
#$ -pe openmp 1

#module load sorenar/miniconda/3
#source activate

module load samtools
module load bcftools

ref='/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/ref/grch38/genome.fa'

#samtools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --output-tags DP,AD -f ${ref} --BCF "/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/GM12878/20190930/G30/G30_sorted.bam" | bcftools call --multiallelic-caller --variants-only -Ob > GM.bcf

samtools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --output-tags DP,AD -f ${ref} --BCF "/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/GM12878/20190930/G10/G10_sorted.bam" | bcftools call --multiallelic-caller --variants-only -Ob > GM2.bcf


bcftools view -Ov GM2.bcf | /data/apps/bcftools/1.9/bin/vcfutils.pl varFilter -d 18 -w 1 -W 3 -a 1 -1 0.05 -2 0.05 -3 0.05 -4 0.05 -e 0.05 -p > GM2.filt.vcf

awk '{if($4=="T" && $5=="C") {print $1"\t"$2"\t"$4"\t"$5}}' GM2.filt.vcf > TCSNP_G10.txt

