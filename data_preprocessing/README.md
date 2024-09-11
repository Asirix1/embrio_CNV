# FASTQ to BigWig pipeline

## Requirements
1) Python >= 3.7 with NumPy
2)
  	- [Juicer Tools](https://github.com/aidenlab/juicer) and Java (for dumping the contacts from the existing .hic-files and/or creating the new .hic-files)
  	- [Cooler](https://github.com/open2c/cooler) (for dumping the contacts from the existing .mcool-files and/or creating the new .mcool-files)
3) the generation of the database (this is a main file required for SV simulation) based on Hi-C map with 600 million contacts requires around 48-60Gb RAM and ~9 hours
4) the simulation of rearrangement for *one chromosome pair* using a pre-computed database requires around 6-36Gb RAM and 1-2 hours
