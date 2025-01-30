import argparse
import torch
from Coverage_Prediction_Model import BertForCoveragePrediction
import json
from transformers import AutoConfig
from pathlib import Path
from hydra import compose, initialize_config_dir
import importlib
import inspect
import pysam
import random
from hydra.utils import instantiate
from lm_experiments_tools.utils import get_cls_by_name
from tqdm import tqdm
import csv
import h5py
import re
import numpy as np
import os
from torch import nn
import torch

def main(keys_path, hdf5_path, experiment_config_path, config_path, checkpoint_path, output_file, batch_size, labels):
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    torch.cuda.empty_cache()
    torch.cuda.synchronize()

    with open(keys_path, 'r') as f:
        keys = [line.strip() for line in f]

    hdf5_file = h5py.File(hdf5_path, 'r')

    with open(config_path, 'r') as f:
        config = json.load(f)

    with initialize_config_dir(str(experiment_config_path.parent)):
        experiment_config = compose(config_name=experiment_config_path.name)

    model_cfg = AutoConfig.from_pretrained(config["model_cfg"])
    print(f"num_labels after AutoConfig.from_pretrained: {model_cfg.num_labels}")

    if "model_kwargs" in experiment_config:
        model_kwargs = instantiate(experiment_config["model_kwargs"])
    else:
        model_kwargs = {}

    if "config" in model_kwargs:
        for k, v in model_kwargs["config"].items():
            model_cfg.__setattr__(k, v)

    print(f"num_labels after applying model_kwargs: {model_cfg.num_labels}")

    model_kwargs["config"] = model_cfg
    print(f"num_labels before model instantiation: {model_cfg.num_labels}")

    model_cfg.num_labels = labels

    model_cls = get_cls_by_name(config["model_cls"])
    importlib.reload(inspect.getmodule(model_cls))
    model = model_cls(**model_kwargs)

    labels_shape = (batch_size, labels)

    model.add_dynamic_layers(labels_shape)
    output_weights_file = "model_weights_and_params.txt"
    with open(output_weights_file, 'w') as f:
        f.write("=== Model Weights and Parameters ===\n")
        for name, param in model.state_dict().items():
            f.write(f"Layer: {name}, Shape: {param.shape}\n")

        checkpoint = torch.load(checkpoint_path, map_location='cpu')
        f.write("\n=== Weights from Checkpoint ===\n")
        for name, param in checkpoint["model_state_dict"].items():
            f.write(f"Layer: {name}, Shape: {param.shape}\n")

        total_params = sum(p.numel() for p in model.parameters())
        trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        non_trainable_params = total_params - trainable_params

        f.write("\n=== Parameter Counts ===\n")
        f.write(f"Total parameters: {total_params}\n")
        f.write(f"Trainable parameters: {trainable_params}\n")
        f.write(f"Non-trainable parameters: {non_trainable_params}\n")

        f.write("\n=== Model Structure ===\n")
        f.write(str(model) + "\n")

    print(f"All weights and parameters saved to: {output_weights_file}")
    print(f"num_labels after model instantiation: {model.num_labels}")

    checkpoint = torch.load(checkpoint_path, map_location='cpu')
    model_state_dict = model.state_dict()
    missing_keys, unexpected_keys = model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    model.eval()

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        header = ['Chromosome', 'Start', 'End'] + [f'Prediction_{i}' for i in range(1, 12)]
        writer.writerow(header)

        batch_keys = []
        batch_inputs = []

        for key in tqdm(keys, desc="Working on keys"):
            chromosome, start, end = key.split(',')
            start, end = int(start), int(end)

            token_key = f'{chromosome},{start},{end}'
            if token_key not in hdf5_file:
                continue

            token_ids = hdf5_file[token_key][()]
            batch_keys.append((chromosome, start, end))
            batch_inputs.append(token_ids)

            if len(batch_inputs) == batch_size:
                inputs = torch.tensor(np.array(batch_inputs)).to(device)
                with torch.no_grad():
                    outputs = model(input_ids=inputs)
                    batch_predictions = outputs.logits.cpu().numpy()

                for i, (chromosome, start, end) in enumerate(batch_keys):
                    row = [chromosome, start, end] + batch_predictions[i].tolist()
                    writer.writerow(row)

                batch_keys = []
                batch_inputs = []

        if batch_inputs:
            inputs = torch.tensor(batch_inputs).to(device)
            with torch.no_grad():
                outputs = model(input_ids=inputs)
                batch_predictions = outputs.logits.cpu().numpy()

            for i, (chromosome, start, end) in enumerate(batch_keys):
                row = [chromosome, start, end] + batch_predictions[i].tolist()
                writer.writerow(row)

    print("Saved to:", output_file)

    hdf5_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run coverage prediction model.")
    parser.add_argument("--keys_path", type=str, required=True, help="Path to the keys CSV file.")
    parser.add_argument("--hdf5_path", type=str, required=True, help="Path to the HDF5 file.")
    parser.add_argument("--experiment_config_path", type=str, required=True, help="Path to the experiment config YAML file.")
    parser.add_argument("--config_path", type=str, required=True, help="Path to the model config JSON file.")
    parser.add_argument("--checkpoint_path", type=str, required=True, help="Path to the model checkpoint file.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output TSV file.")
    parser.add_argument("--batch_size", type=int, required=True, help="Batch size for inference.")
    parser.add_argument("--labels", type=int, required=True, help="Number of labels for the model.")

    args = parser.parse_args()

    main(args.keys_path, args.hdf5_path, Path(args.experiment_config_path), Path(args.config_path), args.checkpoint_path, args.output_file, args.batch_size, args.labels)