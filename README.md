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

# Genomic Coverage Prediction and Segmentation Pipeline

This repository contains a bash script (`cnv_prediction_pipeline.sh`) that automates a pipeline for genomic coverage prediction, real coverage computation, and segmentation using a Hidden Markov Model (HMM). The pipeline consists of three main steps:

1. **Coverage Prediction**: Predicts genomic coverage using a pre-trained model.
2. **Real Coverage Computation**: Computes real coverage values from bigWig files for specified genomic regions.
3. **HMM Segmentation**: Segments the genome using a Hidden Markov Model based on the predicted and real coverage values.



# Usage
## Running the Pipeline
To run the pipeline, use the following command:
```
bash
./cnv_prediction_pipeline.sh \
    --keys_path path/to/keys.csv \
    --hdf5_path path/to/data.hdf5 \
    --experiment_config_path path/to/experiment_config.yaml \
    --labels num_samples \
    --config_path path/to/model_config.json \
    --checkpoint_path path/to/checkpoint.pt \
    --sample_file path/to/samples.txt \
    --regions_file path/to/regions.csv
```
## Command-Line Arguments

| Argument                     | Description                                                                 | Default Value           |
|------------------------------|-----------------------------------------------------------------------------|-------------------------|
| `--keys_path`                | Path to the genomic regions in CSV format used for the model.                                    | **Required**            |
| `--hdf5_path`                | Path to the HDF5 file for the model.                                        | **Required**            |
| `--experiment_config_path`   | Path to the experiment config YAML file for the model.                      | **Required**            |
| `--config_path`              | Path to the model config JSON file.                                         | **Required**            |
| `--checkpoint_path`          | Path to the model checkpoint file.                                          | **Required**            |
| `--model_output_file`        | Output file for the model predictions.                                      | `model_predictions.tsv` |
| `--batch_size`               | Batch size for the model inference.                                         | `32`                    |
| `--labels`                   | Number of samples used to finetune the model.                               | **Required**            |
| `--sample_file`              | Path to the file with sample names and bigWig file paths.                   | **Required**            |
| `--real_coverage_output`     | Output file for real coverage.                                              | `real_coverage.csv`     |
| `--hmm_output_file`          | Output file for HMM predictions.                                            | `hmm_predictions.tsv`   |
| `--help`                     | Show the help message and exit.                                             | N/A                     |

## Output Files

1. **Model Predictions** (`model_predictions.tsv` or custom name):
   - Contains predicted coverage values for the specified genomic regions.

2. **Real Coverage** (`real_coverage.csv` or custom name):
   - Contains computed real coverage values from bigWig files.

3. **HMM Predictions** (`hmm_predictions.tsv` or custom name):
   - Contains the final genome segmentation results from the HMM.

## Script Workflow

1. **Coverage Prediction**:
   - Runs `predictions_by_the_model.py` to generate predicted coverage values.
   - Output: `model_predictions.tsv`.

2. **Real Coverage Computation**:
   - Runs `real_coverage_count.py` to compute real coverage values from bigWig files.
   - Output: `real_coverage.csv`.

3. **HMM Segmentation**:
   - Runs `hmm_segmentator.py` to segment the genome using the predicted and real coverage values.
   - Output: `hmm_predictions.tsv`.

----------------------------------------------------------------------------------------------------------------------------------------------------

# Getting Coverage Prediction by predictions_by_the_model.py

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


# BigWig Real Coverage Calculator - real_coverage_count.py

This Python tool computes the total coverage for each sample over specified genomic regions using bigWig files. It reads sample information from a text file, regions from a CSV file, calculates the coverage for each region using numpy operations, and outputs the results as a CSV file. A progress bar (via `tqdm`) is used to track the computation progress.

## Features

- **Sample Loading:**  
  Reads sample names and corresponding bigWig file paths from a tab-separated file.

- **Region Loading:**  
  Loads genomic regions (chromosome, start, end) from a CSV file.

- **Coverage Calculation:**  
  For each region, computes the sum of coverage for every sample using bigWig files.  
  If an error occurs while reading a sample for a region, it records a coverage value of 0.


- **Output:**  
  Generates a CSV file with columns: `chrom`, `start`, `end`, followed by one column per sample containing the total coverage in that region.
## Usage
Prepare the following input files:

1. Sample File:
A tab-separated text file where each line contains a sample name and the path to its bigWig file.
Example (`samples.txt`):
```bash
sample1	/path/to/sample1.bw
sample2	/path/to/sample2.bw
```
2. Regions File:
A CSV file without header, where each line specifies a genomic region with three columns: chrom, start, end.
Example (`regions.csv`):
```bash
chr1,10000,10500
chr2,20000,20500
```
Run the script from the command line:
```bash
python real_coverage_count.py --samples samples.txt \
--regions regions.csv \
--output coverage_output.csv
```

- `--samples`: Path to the sample file containing sample names and bigWig file paths.
- `--regions`: Path to the CSV file with genomic regions.
- `--output`: Desired output CSV file name.


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

- `--prediction_coverage`: Path to the prediction coverage file (TSV format) - produced by predictions_by_the_model.py.
- `--real_coverage`: Path to the real coverage file (CSV format) - produced by real_coverage_count.py.
- `--output_file`: Path to the output file where all aggregated predictions will be saved (TSV format).
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