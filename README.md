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

(!) There should be the same header as in example CNV_filt.csv.

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
