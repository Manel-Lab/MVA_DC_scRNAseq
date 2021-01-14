#!/bin/bash

#################################
## Genome alignment ans feature count
##################################

#PBS -N singularity_MVADC
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=10
#PBS -l mem=50gb
#PBS -j oe
#PBS -o Singularity_count.log
#PBS -V

##setting up path
DATADIR=$1
cd $DATADIR

##run singularity
singularity exec --bind /data/tmp/kdeazeve/project/:/project MVA_DC_scRNAseq.sif bash -c "/project/count.sh /project/results/mkref/GRCh38_MVA_mcherry_merged/ /project/data/fastqs/D1689_MVA_05/ 1000 /project/torque.template /project/results/101_v4_reprod MVA-rep1"

singularity exec --bind /data/tmp/kdeazeve/project/:/project MVA_DC_scRNAseq.sif bash -c "/project/count.sh /project/results/mkref/GRCh38_MVA_mcherry_merged/ /project/data/fastqs/D1690_MVA_05/ 1000 /project/torque.template /project/results/101_v4_reprod MVA-rep2"

singularity exec --bind /data/tmp/kdeazeve/project/:/project MVA_DC_scRNAseq.sif bash -c "/project/count.sh /project/results/mkref/GRCh38_MVA_mcherry_merged/ /project/data/fastqs/D1689_uninf/ 2000 /project/torque.template /project/results/101_v4_reprod Uninfected-rep1"

singularity exec --bind /data/tmp/kdeazeve/project/:/project MVA_DC_scRNAseq.sif bash -c "/project/count.sh /project/results/mkref/GRCh38_MVA_mcherry_merged/ /project/data/fastqs/D1690_uninf/ 2000 /project/torque.template /project/results/101_v4_reprod Uninfected-rep2"
