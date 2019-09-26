#!/bin/bash
#$ -ckpt restart
#$ -q sam128,sam,abio,abio128,pub64
#$ -pe openmp 4


module load sorenar/miniconda/3
conda activate timelapse_env

#module load hisat2/2.1.0

fastq_path=$1
exp_name=$2


result_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}"
input_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${fastq_path}"
HISAT_ref_genome="/dfs3/samlab/sorenar/OsO-seq/ref/grcm38/genome"
HISAT_ref_trans="/dfs3/samlab/sorenar/OsO-seq/ref/GRCH38/grch38_tran/genome_tran"
STAR_ref=""
STAR_annot=""
aligned_sam=""

#mkdir -p ../results/${fastq_path}

#sorted_bam_file='/dfs3/samlab/sorenar/OsO-seq/timelapse/filtering/'${exp_name}'/'${exp_name}'_sorted.bam'
#mismatch_file='/dfs3/samlab/sorenar/OsO-seq/timelapse/filtering/'${exp_name}'/'${exp_name}'_actb_mismatch.txt'
#ref='/dfs3/samlab/sorenar/OsO-seq/ref/HISAT/grcm38/genome.fa'
#regions='/dfs3/samlab/sorenar/OsO-seq/ref/HISAT/grcm38/mouse_Actb_ENSMUST00000100497_exons.bed'


#ls -ld ${input_path}*.fastq | cut -d ' ' -f 9 > ${input_path}fastq_list.txt

#fastuniq -i ${input_path}fastq_list.txt -o ${result_path}_R1_uniq.fastq -p ${result_path}_R2_uniq.fastq

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 -o ${result_path}_R1_trim.fastq -p ${result_path}_R2_trim.fastq ${result_path}_R1_uniq.fastq ${result_path}_R2_uniq.fastq

# using hisat2 with genome index
#hisat2 --mp 4,2 -p 8 -x ${HISAT_ref_genome} -1 ${result_path}_R1_trim.fastq -2 ${result_path}_R2_trim.fastq -S ${result_path}_Aligned.sam

# using hisat with transcriptome index
#hisat2 --mp 4,2 -p 8 -x ${HISAT_ref_trans} -1 ${result_path}_R1_trim.fastq -2 ${result_path}_R2_trim.fastq -S ${result_path}_Aligned.sam

# using star with genome index and transcriptome annotations



#java -jar /data/apps/picard/2.18.4/picard.jar FixMateInformation I=${result_path}_Aligned.sam O=${result_path}_fixmate.sam

#samtools view -H ${result_path}_fixmate.sam > ${result_path}_header.txt
#samtools view -q 2 -f 83 ${result_path}_fixmate.sam > ${result_path}_filtered_1.sam
#samtools view -q 2 -f 163 ${result_path}_fixmate.sam > ${result_path}_filtered_2.sam
#samtools view -q 2 -f 99 ${result_path}_fixmate.sam > ${result_path}_filtered_3.sam
#samtools view -q 2 -f 147 ${result_path}_fixmate.sam > ${result_path}_filtered_4.sam
#cat ${result_path}_header.txt ${result_path}_filtered_*.sam > ${result_path}_filtered.sam
#samtools sort ${result_path}_filtered.sam > ${result_path}_sorted.bam

#samtools index ${result_path}_sorted.bam

# using bamreadcount to look at the mismatch profiles 
#./mismatch_bam_readcount.sh ${result_path}_sorted.bam ${ref} ${result_path}_mismatch.txt ${regions}

# generating SNP dataset from he control runs
# qsub -N SNPanalysis

# getting the chromosome lengths (use the first 25 for human)
#grep "^@SQ" ${result_path}_header.txt | awk -F '[\t:]' '{print $3"\t"$5}' | head -n 22 > ${result_path}_chr.txt

#module load R/3.6.0

input="${result_path}_chr.txt" 

while IFS= read -r line
do

#chr="15"
#len="104043685"
  chr="$(echo $line | cut -d " " -f1)"
  len="$(echo $line | cut -d " " -f2)"
  qsub -N ${exp_name}_chr${chr} categorizer.sh ${fastq_path} ${exp_name} ${chr} ${len}  
#Rscript mismatch_analysis.R ${fastq_path} ${exp_name} ${chr} ${len} 

done < "$input"


#cat ${result_path}_read_profile_*.txt > ${result_path}_read_profile.txt

#for TC in 0 1 2 3 4
#do
#   qsub -N ${exp_name}_split_${TC} read_split.sh ${result_path} ${TC}
#done

#STAR --runMode inputAlignmentsFromBAM --inputBAMfile ${exp_name}/${exp_name}_${TCN}.bam --outWigType bedGraph


