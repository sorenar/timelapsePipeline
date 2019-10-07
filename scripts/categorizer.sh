#!/bin/bash
#$ -ckpt restart
#$ -q sam128,sam,abio,abio128,pub64
#$ -pe openmp 1
#$ -t 1-25
#$ -tc 5




#module load sorenar/miniconda/3
#source activate

fastq_path=$1
exp_name=$2
result_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}"


module load R/3.6.0

input="${result_path}_chr.txt"
#chrom=$(sed '${SGE_TASK_ID}q;d' ${input})

chrom=$(head -n ${SGE_TASK_ID} ${input}| tail -n 1)

chr="$(echo ${chrom} | cut -d " " -f1)"
len="$(echo ${chrom} | cut -d " " -f2)"

echo ${chr}_${len}
Rscript mismatch_analysis.R ${fastq_path} ${exp_name} ${chr} ${len}
