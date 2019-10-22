#!/bin/bash
#$ -ckpt restart
#$ -q sam,sam128,abio,abio128,pub64
#$ -pe openmp 1

fastq_path=$1
exp=$2
result_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp}"
fastq1=$(sed '1q;d' ../data/${fastq_path}/fastq_list.txt)
fastq2=$(sed '2q;d' ../data/${fastq_path}/fastq_list.txt)

#printf "${fastq1}\n${fastq2}" > ${exp}/exp_stat.tsv

printf "pre_uniq_fastq1\t$(($(wc -l ${fastq1} | cut -d " " -f 1) / 4))\n" > ${result_path}_stat.tsv
printf "pre_uniq_fastq2\t$(($(wc -l ${fastq2} | cut -d " " -f 1) / 4))\n" >> ${result_path}_stat.tsv


printf "uniq_fastq1\t$(($(wc -l ${result_path}_R1_uniq.fastq | cut -d " " -f 1) / 4))\n" >> ${result_path}_stat.tsv
printf "uniq_fastq2\t$(($(wc -l ${result_path}_R2_uniq.fastq | cut -d " " -f 1) / 4))\n" >> ${result_path}_stat.tsv


printf "trim_fastq1\t$(($(wc -l ${result_path}_R1_trim.fastq | cut -d " " -f 1) / 4))\n" >> ${result_path}_stat.tsv
printf "trim_fastq2\t$(($(wc -l ${result_path}_R2_trim.fastq | cut -d " " -f 1) / 4))\n" >> ${result_path}_stat.tsv

module load samtools

printf "total_sam\t$(samtools view ${result_path}_Aligned.sam | cut -f1 | sort | uniq | wc -l | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv

printf "aligned_sam\t$(samtools view -F 4 ${result_path}_Aligned.sam | cut -f1 | sort | uniq | wc -l | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv

printf "fixmate_sam\t$(samtools view ${result_path}_fixmate.sam | cut -f1 | sort | uniq | wc -l | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv

printf "Aligned_fixmate_sam\t$(samtools view -F 4 ${result_path}_fixmate.sam | cut -f1 | sort | uniq | wc -l | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv

printf "filtered_bam\t$(samtools view ${result_path}_sorted.bam | cut -f1 | sort | uniq | wc -l | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv


printf "sub0_read\t$(wc -l ${result_path}_read_list_0.txt | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub1_read\t$(wc -l ${result_path}_read_list_1.txt | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub2_read\t$(wc -l ${result_path}_read_list_2.txt | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub3_read\t$(wc -l ${result_path}_read_list_3.txt | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub4_read\t$(wc -l ${result_path}_read_list_4.txt | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
#printf "sub5_read\t$(wc -l ${result_path}_read_list_5.txt | cut -d " " -f 1)\n" >> ${exp}/exp_stat.tsv


printf "sub0_sam\t$(wc -l ${result_path}_0.sam | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub1_sam\t$(wc -l ${result_path}_1.sam | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub2_sam\t$(wc -l ${result_path}_2.sam | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub3_sam\t$(wc -l ${result_path}_3.sam | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
printf "sub4_sam\t$(wc -l ${result_path}_4.sam | cut -d " " -f 1)\n" >> ${result_path}_stat.tsv
#printf "sub5_sam\t$(wc -l ${result_path}_5.sam | cut -d " " -f 1)\n" >> ${exp}/exp_stat.tsv


printf "total_count0\t$(cut -f 2 ${result_path}_0_counts.txt | head -n -5 | paste -sd+ | bc)\n" >> ${result_path}_stat.tsv
printf "total_count1\t$(cut -f 2 ${result_path}_1_counts.txt | head -n -5 | paste -sd+ | bc)\n" >> ${result_path}_stat.tsv
printf "total_count2\t$(cut -f 2 ${result_path}_2_counts.txt | head -n -5 | paste -sd+ | bc)\n" >> ${result_path}_stat.tsv
printf "total_count3\t$(cut -f 2 ${result_path}_3_counts.txt | head -n -5 | paste -sd+ | bc)\n" >> ${result_path}_stat.tsv
printf "total_count4\t$(cut -f 2 ${result_path}_4_counts.txt | head -n -5 | paste -sd+ | bc)\n" >> ${result_path}_stat.tsv
#printf "total_count5\t$(cut -f 2 ${result_path}_5_counts.txt | head -n -5 | paste -sd+ | bc)\n" >> ${exp}/exp_stat.tsv


printf "nofeature_count0\t$(grep no_feature ${result_path}_0_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "nofeature_count1\t$(grep no_feature ${result_path}_1_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "nofeature_count2\t$(grep no_feature ${result_path}_2_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "nofeature_count3\t$(grep no_feature ${result_path}_3_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "nofeature_count4\t$(grep no_feature ${result_path}_4_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
#printf "nofeature_count5\t$(grep no_feature ${result_path}_5_counts.txt | cut -f 2)\n" >> ${exp}/exp_stat.tsv


printf "ambiguous_count0\t$(grep ambiguous ${result_path}_0_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "ambiguous_count1\t$(grep ambiguous ${result_path}_1_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "ambiguous_count2\t$(grep ambiguous ${result_path}_2_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "ambiguous_count3\t$(grep ambiguous ${result_path}_3_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
printf "ambiguous_count4\t$(grep ambiguous ${result_path}_4_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv
#printf "ambiguous_count5\t$(grep ambiguous ${result_path}_5_counts.txt | cut -f 2)\n" >> ${result_path}_stat.tsv




