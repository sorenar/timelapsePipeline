#!/bin/bash

input=$1

touch tracks.txt

while IFS= read -r line
do

  exp_name="$(echo $line | cut -d " " -f1)"
  exp_dir="$(echo $line | cut -d " " -f2)"
  exp_desc="$(echo $line | cut -d " " -f3)"
  mkdir -p /pub/public-www/sorenar/bws/${exp_dir}
  cp /dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/${exp_dir}${exp_name}*.bw /pub/public-www/sorenar/bws/${exp_dir}
  printf "track type=bigWig name=\"${exp_desc}_0\" visibility=full color=220,220,220 bigDataUrl=https://hpc.oit.uci.edu/~sorenar/bws/${exp_dir}${exp_name}_0.bw\n" >> tracks.txt
  printf "track type=bigWig name=\"${exp_desc}_1\" visibility=full color=230,200,230 bigDataUrl=https://hpc.oit.uci.edu/~sorenar/bws/${exp_dir}${exp_name}_1.bw\n" >> tracks.txt
  printf "track type=bigWig name=\"${exp_desc}_2\" visibility=full color=220,160,220 bigDataUrl=https://hpc.oit.uci.edu/~sorenar/bws/${exp_dir}${exp_name}_2.bw\n" >> tracks.txt
  #printf "track type=bigWig name=\"${exp_desc}_3\" visibility=full color=240,110,240 bigDataUrl=https://hpc.oit.uci.edu/~sorenar/bws/${exp_dir}${exp_name}_3.bw\n" >> tracks.txt
  #printf "track type=bigWig name=\"${exp_desc}_4\" visibility=full color=240,10,240 bigDataUrl=https://hpc.oit.uci.edu/~sorenar/bws/${exp_dir}${exp_name}_4.bw\n" >> tracks.txt

done < "$input"


