#!/usr/bin/env bash
set -e

cd ../..

TBS=64
BS=32
NP=1

export CUDA_HOME="$HOME/.local/cuda/"
export PATH="$HOME/.local/cuda/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/.local/cuda/lib64/:$LD_LIBRARY_PATH"

CUDA_VISIBLE_DEVICES=1 python /mnt/nfs_dna/popov/CNV_prediction/Models_scripts/Coverage_Prediction_Finetune.py \ 
		--experiment_config "CoveragePrediction.yaml" \  
        --batch_size $BS --gradient_accumulation_steps $(($TBS/($BS*$NP)))

echo "done"