#!/bin/bash

data_dir=$1

find "$data_dir" -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) -delete

