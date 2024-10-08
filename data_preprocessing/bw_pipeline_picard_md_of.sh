#!/bin/bash

fastq_dir=$1
reference_path=$2
work_dir=$3
threads=$4
MAPQ=$5
sample=$6
md_py=$7

time bwa mem -t $threads $reference_path $fastq_dir/*R1*.fastq.gz $fastq_dir/*R2*.fastq.gz |samtools view  -b -@ $threads -h -F 4 |samtools sort -@ $threads -o $work_dir/${sample}_sorted_pre.bam

cd ./picard
time java -jar build/libs/picard.jar AddOrReplaceReadGroups -I ${work_dir}/${sample}_sorted_pre.bam -O ${work_dir}/${sample}_sorted_with_rg.bam -RGID group1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM sample1

time java -jar build/libs/picard.jar MarkDuplicates -I ${work_dir}/${sample}_sorted_with_rg.bam -O ${work_dir}/${sample}_marked_duplicates.bam -M ${work_dir}/${sample}_marked_dup_metrics.txt --REMOVE_DUPLICATES true

time samtools index -@ $threads ${work_dir}/${sample}_marked_duplicates.bam

time bamCoverage -b ${work_dir}/${sample}_marked_duplicates.bam -o $work_dir/${sample}_md_offset.bw -p $threads --minMappingQuality $MAPQ --Offset 1 1 --binSize 1
