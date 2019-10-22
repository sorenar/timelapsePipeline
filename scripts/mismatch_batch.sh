#!/bin/bash

input=$1

while IFS= read -r line
do
  #old_dir="$(echo $line | cut -d " " -f1)"
  exp_name="$(echo $line | cut -d " " -f1)"
  exp_dir="$(echo $line | cut -d " " -f2)"

  #mkdir -p /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${exp_dir}/
  #mv /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${exp_name}* /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${exp_dir}/
  #gunzip /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${exp_dir}/*
 # rm /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${new_dir}*
 # cp /dfs3/samlab/sorenar/OsO-seq/timelapse/${old_dir}/${old_dir}_R*.fastq /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${new_dir} 
  #gunzip /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${exp_dir}*.gz
  qsub -N $exp_name pipeline.sh ${exp_dir} ${exp_name}
  #qsub -N $exp_name read_stat.sh ${exp_dir} ${exp_name}
done < "$input"

