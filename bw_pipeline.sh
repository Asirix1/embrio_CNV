#!/bin/bash

fastq_dir=$1
reference_path=$2
work_dir=$3
threads=$4
MAPQ=$5
sample=$6

time bwa mem -t $threads $reference_path $fastq_dir/*R1*.fastq.gz $fastq_dir/*R2*.fastq.gz |samtools view  -b -@ $threads -h -F 4 |samtools sort -@ $threads -o $work_dir/${sample}_sorted.bam

time samtools index -@ $threads $work_dir/${sample}_sorted.bam

time bamCoverage -b $work_dir/${sample}_sorted.bam -o $work_dir/${sample}.bw -p $threads --minMappingQuality $MAPQ
