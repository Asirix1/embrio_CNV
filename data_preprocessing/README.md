# FASTQ to BigWig pipeline

## Requirements
1) Python >= 3.10.14 with Pandas and OS
2) Java = 17   (should be installed befor picard istallation)
3) 
  	- [bwa](https://github.com/open2c/cooler) (for aligning genome) (it can be istallad with ``` conda install bwa -c bioconda ```)
  	- [picard](https://github.com/broadinstitute/picard)) (for marking duplicates)
    - [samtools](https://github.com/samtools/samtools) (for sorting genome) (it can be istallad with ``` conda install bwa -c samtools ```)

# Quick Start

1) Create your csv-file with downloads links
Example:


```
Sample (е-эмбрион. к-биоптат),R1,R2,,,,
[test,trn,chr1,0,,1000000,!>,der1,1](https://genedev.bionet.nsc.ru/ftp/_RawReads/2022-11-07_BGI_Moscow/Demultiplexed/embryo6/embryo6_R1.fastq.gz,https://genedev.bionet.nsc.ru/ftp/_RawReads/2022-11-07_BGI_Moscow/Demultiplexed/embryo6/embryo6_R2.fastq.gz),,,,
test,trn,chr1,2000000,7000000,>>,der1,1
```
3) Duplicate [EXAMPLE.ini file](testdataset/EXAMPLE.ini) and modify:
