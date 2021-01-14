#!/bin/bash
# author: Francesca Nadalin
# usage: qsub cellranger_count.sh -F "$1 $2 $3 $4 $5" / bash cellranger_count.sh $1 $2 $3 $4 $5
# goal : Alignment and cell counting with cellranger 4.0.0

#### REQUIREMENTS ####
# CellRanger 4.0.0
######################

################################################
## PBS ARGUMENTS
################################################

#PBS -V
#PBS -l nodes=1:ppn=10
#PBS -l mem=50gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o /data/tmp/kdeazeve/count.log
#PBS -N count

################################################
## ARGUMENTS AND EXAMPLE
################################################

# Arguments
TRANSCRIPTOME=$1
CELLRANGER_PATH="/usr/local/cellranger-4.0.0/cellranger"
FASTQS=$2
EXPECT_CELLS=$3
TEMPLATE=$4
OUT_DIR=$5
SAMPLE=$6

################################################
## SCRIPT
################################################

[ -d $OUT_DIR ] || mkdir -p $OUT_DIR
cd $OUT_DIR
echo Running on `hostname`
echo Start time is `date`
echo Directory is `pwd`
$CELLRANGER_PATH count --id $SAMPLE --project DCINNATERESPMVA --fastqs $FASTQS --transcriptome "$TRANSCRIPTOME" \
--expect-cells $EXPECT_CELLS --localmem=50 --localcores=10 --nosecondary --chemistry SC3Pv2
echo End time is `date`
#qstat -f $PBS_JOBID > "cluster_cellranger_count.log"
