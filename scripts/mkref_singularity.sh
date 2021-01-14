#!/bin/bash

#################################
## Indexing the reference transcriptome
##################################

#PBS -N singularity_MVADC
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=5
#PBS -l mem=40gb
#PBS -j oe
#PBS -o Singularity_mkref.log
#PBS -V

##setting up path
DATADIR=$1
cd $DATADIR

##run singularity
singularity exec --bind /data/tmp/kdeazeve/project/:/project MVA_DC_scRNAseq.sif bash -c "/project/mkref.sh GRCh38_MVA_mcherry_merged /project/results/mkref/GRCh38_MVA_mcherry_merged.fa /project/results/mkref/GRCh38_MVA_mcherry_merged.gtf /project/results/mkref"
