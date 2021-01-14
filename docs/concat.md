# Creation of the concatenated genome

## Introduction

The bash script ``create_merged_genome.sh`` download the human and MVA genome aswell as the annotation, and concatenate them with the mcherry sequence to create our reference genome used in the analysis. The mcherry sequence and gtf have to be provided through the arguments. The arguments are :

1. Fasta file with the mcherry sequence
2. Gtf of the mcherry sequence
3. Output directory of the created genome

The fasta and gtf of the mcherry sequence are in the ``data`` directory.

## Usage

``bash create_merged_genome.sh /path/to/mcherry.fa /path/to/mcherry.gtf /path/to/output``
