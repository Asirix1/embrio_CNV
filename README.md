# Embryo CNV detection

# FASTQ to BigWig pipeline (data preprocessing)

## Requirements
1) Python >= 3.10.14 with Pandas and OS
2) Java = 17   (should be installed befor picard istallation) (it can be installad with ``` mamba install openjdk=17.0.3 ```)
3) 
  	- [bwa](https://github.com/open2c/cooler) (for aligning genome) (it can be installad with ``` mamba install bwa -c bioconda ```)
    - [samtools](https://github.com/samtools/samtools) (for sorting genome) (it can be installad with ``` mamba install samtools -c bioconda ```)
  	- [picard](https://github.com/broadinstitute/picard) (for marking duplicates)

*picard* can be installed as follow:
```
git clone https://github.com/broadinstitute/picard.git
```
```
cd picard/
```
```
./gradlew shadowJar
```

Check picard is installed:
```
cd picard/
```
```
java -jar build/libs/picard.jar
```
For more information see https://github.com/broadinstitute/picard .

# Quick Start

1) Create your csv-file as in example CNV_filt.csv with download links to R1 and R2 in *.fastq.gz* or *fq.gz* format (characters *R1* or *R2* in filenames are not obligatory).

(!) There should be the same header as in example CNV_example.csv.

You can add several acts of sequencing or sequencing data from several slots with the same sample name. In this case file pairs will be merged by pipeline.

For test run use CNV_example.csv.

3) Create a work directory (work_dir).
4) Run
```
cd embrio_CNV/data_preprocessing
python3 python_script_writer.py -csv CNV_example.csv -d work_dir -bwp bw_pipeline_picard_md_of.sh -t threads -g reference_genome.fa -q required_MAPQ
```
5) Look output shell-scripts in the directory work_dir/scripts. There should be one script to each sample.
6) You can run scripts separately:
```
bash sample.sh
```
or you can run scripts for all samples at once:
```
bash scr_pr.sh work_dir/scripts
```
7) Output *.bam* and *.bw* files will be located in work_dir/data/bw .
8) You can remove all fastq-files from your directory with
```
bash remove_fastq.sh dir
```
8) You can remove all *.bam* and *.bam.bai* files from your directory with
```
bash remove_bam_bai.sh dir
```
10) You can remove all bw-files from your directory with
```
bash remove_bw.sh dir
```
11) You can remove all *metrics.txt* files from your directory with
```
bash remove_metrics.sh dir
``` 

# Coverage Prediction Model

This repository contains code for fine-tuning a coverage prediction model using PyTorch and BigWig data. The model and dataset configurations are customizable, with training and validation processes logged to TensorBoard.

## Setup

To set up the environment and install all required packages:

1. Create a Conda environment from the configuration file:
   ```bash
   conda env create -f cnv_prediction.yml
   ```
## Configuration
Before training, update the CoveragePrediction.yaml file:

Set HOME_PATH to your working directory path where all training files will be stored.
And unzip file in */data/bigwigs/* and */data/*
## Training the Model
To start model fine-tuning, use the following command in the terminal:
```
python Coverage_Prediction_Finetune.py --experiment_config CoveragePrediction.yaml --batch_size <your_batch_size>
```
## Fine-tuning Data Configuration
To use custom data for fine-tuning, modify the *CoveragePrediction.yaml* file:
* train_dataset and valid_dataset sections:
	* targets_path: Path to a text file listing BigWig files for sample coverage. This file should be formatted with:
		* Column 1: Sample name
		* Column 2: Path to BigWig file for that sample
	* key_file_path: Path to your key file for train and validation keys, formatted as **chr, start, end** with sequences of 510 tokens. You can randomly select keys from the file /data/test_for_git_tokens.hdf5 for convenience.

In *CoveragePrediction.yaml*, you may also customize:
* save_interval: Model save interval
* valid_interval: Validation interval
* or other parameters
## Model Output
* Models are saved in the Models directory, with a uniquely named folder created for each fine-tuning session.
## Monitoring Training Progress
```
tensorboard --logdir <model_save_directory> --port <port>
```
# Coverage Prediction Model by predictions_by_the_model.py

This script runs a coverage prediction model using a pre-trained transformer model. It loads input data from an HDF5 file and processes it in batches to generate predictions.


## Usage

Run the script predictions_by_the_model.py with required parameters:

```bash
python predictions_by_the_model.py \
  --keys_path /path/to/sorted_keys.csv \
  --hdf5_path /path/to/combined_510_tokens.hdf5 \
  --experiment_config_path /path/to/CoveragePrediction.yaml \
  --config_path /path/to/config.json \
  --checkpoint_path /path/to/model_best.pth \
  --output_file results.tsv \
  --batch_size 512 \
  --labels 10
```

## Arguments

- `--keys_path`: Path to the CSV file containing the keys.
- `--hdf5_path`: Path to the HDF5 file with tokenized genome data.
- `--experiment_config_path`: Path to the YAML configuration file for the experiment.
- `--config_path`: Path to the JSON configuration file for the model.
- `--checkpoint_path`: Path to the model checkpoint file.
- `--output_file`: Output file name for storing predictions.
- `--batch_size`: Number of samples to process at once (default: 512).
- `--labels`: Number of labels (samples) for model prediction (default: 10).

## Output

The script produces:

- A TSV file (`output_file`) with predictions in the following format:
  ```
  Chromosome  Start  End  Sample_1  Sample_2  ...
  chr1        1000   2000  0.12          0.45         ...
  ```

## Implementation Details

- Loads the model using `transformers.AutoConfig`.
- Uses `torch.load()` to restore weights from the checkpoint.
- Processes data in batches for efficiency.
- Uses `torch.no_grad()` for inference to reduce memory usage.
- Handles GPU acceleration automatically.

# HMM Segmentator - hmm_segmentator.py

This code implements a Hidden Markov Model (HMM)-based segmentation tool for genomic coverage data. The tool processes both predicted and real coverage data, computes differences, filters regions based on a provided BED file, and trains an HMM with fixed transitions across different dispersion thresholds. The resulting segmentation predictions for each threshold are aggregated and saved into a single output file.

## Features

- **Data Preprocessing:**  
  Merges predicted and real coverage data, calculates a difference metric, and applies filtering based on excluded regions provided in a BED file (e.g., `T2T.excluderanges.bed`).

- **Row Combination:**  
  Combines rows into groups (e.g., every 50 rows) and computes mean values for subsequent analysis.

- **Visualization:**  
  Generates scatter plots to visualize normalized values alongside the predicted states. The plots are created and then closed immediately to prevent display overload.

- **Aggregation of Predictions:**  
  Aggregates segmentation predictions from multiple dispersion thresholds into one final output file.
## Usage
Run the script from the command line with the following parameters:

- --prediction_coverage: Path to the prediction coverage file (TSV format) - produced by predictions_by_the_model.py.
- --real_coverage: Path to the real coverage file (CSV format) - produced by real_coverage_count.py.
- --output_file: Path to the output file where all aggregated predictions will be saved (TSV format).
Example command:
```bash
python hmm_segmentator.py --prediction_coverage [path/to/coverage_predictions.tsv] \
--real_coverage [path/to/real_coverage.csv] \
--output_file [path/to/output_predictions.tsv]
```
## Output
The final output is a TSV file containing aggregated segmentation predictions. Each row includes:

- Column: Base name of the processed prediction column.
- Chromosome: Chromosome identifier.
- Start/End: Start and end positions of the aggregated segment.
- Class: Predicted class (1 - depletion, 3 - duplication).
- Threshold: The dispersion threshold (disp value) used during HMM prediction.