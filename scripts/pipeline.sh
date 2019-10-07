#!/bin/bash
#$ -ckpt restart
#$ -q sam128,sam,abio,abio128,pub64
#$ -pe openmp 1


module load sorenar/miniconda/3
source activate

#module load hisat2/2.1.0

fastq_path=$1
exp_name=$2


result_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${fastq_path}${exp_name}"
input_path="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${fastq_path}"
# mouse references:
#HISAT_ref_genome="/dfs3/samlab/sorenar/OsO-seq/ref/grcm38/genome"
#HISAT_ref_trans="/dfs3/samlab/sorenar/OsO-seq/ref/GRCm8/grch38_tran/genome_tran"

#human refs:
HISAT_ref_genome="/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/ref/grch38/genome"
#HISAT_ref_trans="/dfs3/samlab/sorenar/OsO-seq/ref/GRCH38/grch38_tran/genome_tran"

STAR_ref="/dfs3/samlab/sorenar/OsO-seq/ref/hg38"
STAR_annot="/dfs3/samlab/sorenar/OsO-seq/ref/gencode.v30.primary_assembly.annotation.gtf"
aligned_sam=""

#mkdir -p ../results/${fastq_path}

#sorted_bam_file='/dfs3/samlab/sorenar/OsO-seq/timelapse/filtering/'${exp_name}'/'${exp_name}'_sorted.bam'
#mismatch_file='/dfs3/samlab/sorenar/OsO-seq/timelapse/filtering/'${exp_name}'/'${exp_name}'_actb_mismatch.txt'
#ref='/dfs3/samlab/sorenar/OsO-seq/ref/HISAT/grcm38/genome.fa'
#regions='/dfs3/samlab/sorenar/OsO-seq/ref/HISAT/grcm38/mouse_Actb_ENSMUST00000100497_exons.bed'


#ls -ld ${input_path}*.fastq | cut -d ' ' -f 10 > ${input_path}fastq_list.txt

#fastuniq -i ${input_path}fastq_list.txt -o ${result_path}_R1_uniq.fastq -p ${result_path}_R2_uniq.fastq

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 -o ${result_path}_R1_trim.fastq -p ${result_path}_R2_trim.fastq ${result_path}_R1_uniq.fastq ${result_path}_R2_uniq.fastq

# using hisat2 with genome index
#hisat2 --mp 4,2 -p 8 -x ${HISAT_ref_genome} -1 ${result_path}_R1_trim.fastq -2 ${result_path}_R2_trim.fastq -S ${result_path}_Aligned.sam

# using hisat with transcriptome index
#hisat2 --mp 4,2 -p 8 -x ${HISAT_ref_trans} -1 ${result_path}_R1_trim.fastq -2 ${result_path}_R2_trim.fastq -S ${result_path}_Aligned.sam

# using star with genome index and transcriptome annotations

#module load STAR/2.6.0c

#STAR --runThreadN 8 --genomeDir ${ref} --readFilesIn ${result_path}_R1_trim.fastq ${result_path}_R2_trim.fastq --sjdbGTFfile ${STAR_annot} --outFileNamePrefix ${result_path}_ --outFilterMismatchNmax 15 --outFilterMismatchNoverReadLmax 0.07 --outFilterMultimapNmax 10 --outSAMunmapped None --outSAMattributes MD NM --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM

#mv ${result_path}_Aligned.sortedByCoord.out.bam ${result_path}_Aligned.sam

#java -jar /data/apps/picard/2.18.4/picard.jar FixMateInformation I=${result_path}_Aligned.sam O=${result_path}_fixmate.sam

#samtools view -H ${result_path}_fixmate.sam > ${result_path}_header.txt
#samtools view -q 2 -f 83 -F 256 ${result_path}_fixmate.sam > ${result_path}_filtered_1.sam
#samtools view -q 2 -f 163 -F 256 ${result_path}_fixmate.sam > ${result_path}_filtered_2.sam
#samtools view -q 2 -f 99 -F 256 ${result_path}_fixmate.sam > ${result_path}_filtered_3.sam
#samtools view -q 2 -f 147 -F 256 ${result_path}_fixmate.sam > ${result_path}_filtered_4.sam
#cat ${result_path}_header.txt ${result_path}_filtered_*.sam > ${result_path}_filtered.sam
#samtools sort ${result_path}_filtered.sam > ${result_path}_sorted.bam

#samtools index ${result_path}_sorted.bam

# using bamreadcount to look at the mismatch profiles 
#./mismatch_bam_readcount.sh ${result_path}_sorted.bam ${ref} ${result_path}_mismatch.txt ${regions}

# generating SNP dataset from he control runs
# qsub -N SNPanalysis

# getting the chromosome lengths (use the first 25 for human)
grep "^@SQ" ${result_path}_header.txt | awk -F '[\t:]' '{print $3"\t"$5}' | head -n 25 > ${result_path}_chr.txt

#module load R/3.6.0

#input="${result_path}_chr.txt" 

#while IFS= read -r line
#do

#chr="8"
#len="145138636"
#begin="127736594"
#end="127740958"
#  chr="$(echo $line | cut -d " " -f1)"
#  len="$(echo $line | cut -d " " -f2)"
qsub -l h_vmem=2G -N ${exp_name}_R categorizer.sh ${fastq_path} ${exp_name} ${chr}  
# Rscript mismatch_analysis.R ${fastq_path} ${exp_name} ${chr} ${len} 

# using specific region
#Rscript mismatch_analysis.R ${fastq_path} ${exp_name} ${chr} ${begin} ${end}

#done < "$input"


#cat ${result_path}_read_profile_*.txt > ${result_path}_read_profile.txt

for TC in 0 1 2 3 4
do
#  qsub -N ${exp_name}_split_${TC} read_split.sh ${fastq_path} ${exp_name} ${TC} ${chr}
done

#STAR --runMode inputAlignmentsFromBAM --inputBAMfile ${exp_name}/${exp_name}_${TCN}.bam --outWigType bedGraph



