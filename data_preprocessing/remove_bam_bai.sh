#!/bin/bash

data_dir=$1

find "$data_dir" -type f \( -name "*.bam" -o -name "*.bam.bai" \) -delete