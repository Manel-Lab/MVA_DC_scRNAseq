# Alignment and feature count

## Introduction

The bash script ``count.sh`` aligns the sequencing data to the reference transcritpome and count the number of UMI for each of our feature, human or viral, for one sample. The arguments are :

1. Path to transcritpome of reference
2. Path to the fastq files for a particular sample. Fastq files for this project are available through SRA (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153835). Each sample will have three files (R1, R2 and I1)
3. Number of cells expected. 1000 for MVA samples, 2000 for Uninfected samples.
4. Torque template file (to launch CellRanger count in cluster mode)
5. Output directory
6. Name of the sample

## Usage

You can execute it through the bash script inside or outside the container, for each sample :

``bash count.sh /path/to/GRCh38_MVA_mcherry_merged /path/to/Sample1.fa 1000 /path/to/template_file /path/to/output Sample1``

Or using the bash script provided to run it from outside the container, for all sample :

``bash count_singularity.sh /abolute/path/to/project``
