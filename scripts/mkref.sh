#!/bin/bash
# author: Francesca Nadalin
# usage: qsub run_cellranger_mkref_merged.sh -F "$1 $2 $3 $4 $5 $6 $7" / bash run_cellranger_mkref_merged.sh $1 $2 $3 $4 $5 $6 $7
# goal : Create a custom reference for CellRanger by merging the human and the MVA genome

#### REQUIREMENTS ####
# CellRanger 4.0.0
# gffread
######################

################################################
## PBS ARGUMENTS
################################################

#PBS -V
#PBS -l nodes=1:ppn=5
#PBS -l mem=40gb
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -o /data/tmp/kdeazeve/mkref.log
#PBS -N mkref

################################################
## ARGUMENTS AND EXAMPLE
################################################

name=$1
genome=$2
gtf=$3
OUT_DIR=$4
CELLRANGER_PATH="/usr/local/cellranger-4.0.0/cellranger"

################################################
## SCRIPT
################################################

[ -d $OUT_DIR ] || mkdir -p $OUT_DIR
cd $OUT_DIR

echo Running on `hostname`
echo Start time is `date`
echo Directory is `pwd`
$CELLRANGER_PATH mkref --nthreads=5 --memgb=30 --genome=$name --fasta=$genome --genes=$gtf
cp $gtf "${OUT_DIR}/${name}/genes/genes.gtf"
echo End time is `date`
#qstat -f $PBS_JOBID > "cluster_mkref.log"
