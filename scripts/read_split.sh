#!/bin/bash
#$ -ckpt restart
#$ -q sam,sam128,abio,abio128,pub64
#$ -pe openmp 1



module load sorenar/miniconda/3
source activate

result_path=$1
TCN=$2

annot_file="/dfs3/samlab/sorenar/OsO-seq/timelapse_pipeline/ref/Mus_musculus.GRCm38.97.gtf"

awk -v x=${TCN} '$2>=x {print $1}' ${result_path}_read_profile.txt > ${result_path}_read_list_${TCN}.txt

#samtools view ${exp_name}/${exp_name}_sorted.bam | grep -f ${exp_name}/${exp_name}_read_list_${TCN}.txt > ${exp_name}/${exp_name}_${TCN}.sam

awk 'NR==FNR{c[$1]++;next};c[$1]' ${result_path}_read_list_${TCN}.txt ${result_path}_filtered.sam > ${result_path}_${TCN}.sam

cat ${result_path}_header.txt ${result_path}_${TCN}.sam > ${result_path}_temp_${TCN}.sam

samtools view -bh ${result_path}_temp_${TCN}.sam > ${result_path}_${TCN}.bam

samtools sort ${result_path}_${TCN}.bam > ${result_path}_${TCN}_sorted.bam

samtools index ${result_path}_${TCN}_sorted.bam

bamCoverage -b ${result_path}_${TCN}_sorted.bam -o ${result_path}_${TCN}.bw

htseq-count -f bam -r pos -s no -i transcript_id ${result_path}_${TCN}_sorted.bam ${annot_file} > ${result_path}_${TCN}_counts.txt



