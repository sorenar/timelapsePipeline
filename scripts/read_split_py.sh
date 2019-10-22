#!/bin/bash
#$ -ckpt restart
#$ -q sam,sam128,abio,abio128,pub64
#$ -pe openmp 1



module load sorenar/miniconda/3
source activate

fastq_path=$1
exp_name=$2
TCN=$3
#chr=$4

#result_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${chr}/${exp_name}"

result_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}"

#mkdir -p /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${chr}

#for mouse
# annot_file="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/ref/Mus_musculus.GRCm38.97.gtf"

#for human
annot_file="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/ref/gencode.v30.primary_assembly.annotation.gtf"

#chr-specific
#awk -v x=${TCN} '$2>=x {print $1}' /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}_read_profile.txt > ${result_path}_read_list_${TCN}.txt

#awk -v x=${TCN} '$2>=x {print $1}' ${result_path}_read_profile.txt > ${result_path}_read_list_${TCN}.txt

#samtools view /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}_sorted.bam | grep -f ${result_path}_read_list_${TCN}.txt > ${result_path}_${TCN}.sam

#awk 'NR==FNR{c[$1]++;next};c[$1]' ${result_path}_read_list_${TCN}.txt /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}_filtered.sam > ${result_path}_${TCN}.sam

#cat /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}_header.txt ${result_path}_${TCN}.sam > ${result_path}_temp_${TCN}.sam


samtools view -bh ${result_path}_py_${TCN}.sam > ${result_path}_${TCN}.bam

samtools sort ${result_path}_${TCN}.bam > ${result_path}_${TCN}_sorted.bam

samtools index ${result_path}_${TCN}_sorted.bam

bamCoverage -b ${result_path}_${TCN}_sorted.bam -o ${result_path}_${TCN}.bw

htseq-count -f bam -r pos -s no -i transcript_id ${result_path}_${TCN}_sorted.bam ${annot_file} > ${result_path}_${TCN}_counts.txt



