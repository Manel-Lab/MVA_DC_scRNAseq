#!/bin/bash

#################################
## Script to run the downstream R analysis
##################################

#PBS -N singularity_Rscript_MVA
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -j oe
#PBS -o Singularity_Rscript.log
#PBS -V

## Setting up path
DATADIR=$1
cd $DATADIR

##run singularity
singularity exec --bind /data/tmp/kdeazeve/project/:/project --no-home MVA_DC_scRNAseq.sif bash -c "Rscript --vanilla /project/scripts/mvadcs_downstream.R /project/results/counts /project/data/MVA_gene_map.csv /project/data/mva-class.tsv TRUE /project/results/Rscript args.json /project/scripts/downstream_analysis"
