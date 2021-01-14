# Indexing the transcritpome

## Introduction

The bash script ``mkref.sh`` indexes the transcriptome of reference for the analysis with CellRanger's mkref. The arguments are :

1. Name you want to give to the transcritpome of reference
2. Path to the concatenated fasta
3. Path to the concatenated gtf
4. Output directory

## Usage

You can execute it through the bash script inside or outside the container :

``bash mkref.sh GRCh38_MVA_mcherry_merged /path/to/GRCh38_MVA_mcherry_merged.fa /path/to/GRCh38_MVA_mcherry_merged.gtf /path/to/output``

Or using the bash script provided to run it from outside the container. For this to work the results of the genome concatenation must be in ``/abolute/path/to/project/results`` :

``bash mkref_singularity.sh /abolute/path/to/project``
