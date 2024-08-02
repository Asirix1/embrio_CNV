#!/bin/bash

fastq_dir=$1
reference_path=$2
work_dir=$3
threads=$4



time bwa mem -t $threads -o $work_dir/out.sam $fastq_dir/*/_R1.fastq $fastq_dir/*/_R2.fastq 

time samtools sort -@ $threads $work_dir/out.sam -o $work_dir/out_sorted.sam

time samtools index -@ $threads $work_dir/out_sorted.sam.bai

time bamCoverage -p $threads $work_dir/out_sorted.sam -o $work_dir/out.bw

