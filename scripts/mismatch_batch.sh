#!/bin/bash

input=$1

while IFS= read -r line
do
  old_dir="$(echo $line | cut -d " " -f1)"
  exp_name="$(echo $line | cut -d " " -f2)"
  new_dir="$(echo $line | cut -d " " -f3)"

 # mkdir -p /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${new_dir}
 # rm /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${new_dir}*
 # cp /dfs3/samlab/sorenar/OsO-seq/timelapse/${old_dir}/${old_dir}_R*.fastq /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/data/${new_dir} 
  qsub -N $exp_name pipeline.sh ${new_dir} ${exp_name}

done < "$input"

