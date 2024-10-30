# embrio_CNV

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
