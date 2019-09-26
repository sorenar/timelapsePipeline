#!/bin/bash
#$ -ckpt restart
#$ -q sam128,sam,abio,abio128,pub64
#$ -pe openmp 4


#module load sorenar/miniconda/3
#source activate

result_path=$1
exp_name=$2
chr=$3
len=$4

module load R/3.6.0

Rscript mismatch_analysis.R ${result_path} ${exp_name} ${chr} ${len} 



