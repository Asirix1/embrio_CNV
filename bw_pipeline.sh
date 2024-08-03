#!/bin/bash

fastq_dir=$1
reference_path=$2
work_dir=$3
threads=$4

time bwa mem -t $threads $reference_path $fastq_dir/*/_R1.fastq $fastq_dir/*/_R2.fastq |samtools view  -b -@ $threads -h -F 4 > $work_dir/out.bam

time samtools sort -@ $threads $work_dir/out.bam -o $work_dir/out_sorted.bam

time samtools index -@ $threads $work_dir/out_sorted.bam.bai

time bamCoverage -p $threads $work_dir/out_sorted.bam.bai -o $work_dir/out.bw
