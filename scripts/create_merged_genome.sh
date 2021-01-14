#!/usr/bin/bash
# author: Francesca Nadalin
# usage: qsub run_cellranger_mkref_merged.sh -F "$1 $2 $3 $4 $5 $6 $7" / bash run_cellranger_mkref_merged.sh $1 $2 $3 $4 $5 $6 $7
# goal : Create a custom reference for CellRanger by merging the human and the MVA genome

#### REQUIREMENTS ####
# CellRanger 4.0.0
######################

################################################
## PBS ARGUMENTS
################################################

#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -o /data/tmp/kdeazeve/create_merged_genome.log
#PBS -N create_merged_genome

################################################
## ARGUMENTS AND EXAMPLE
################################################

genome1="GRCh38.fa"
annotation1="GRCh38.gtf"
genome2="MVA.fa"
annotation2="MVA.gtf"
genome3=$1
annotation3=$2
OUT_DIR=$3
genome1_name="GRCh38"
genome2_name="MVA"
genome3_name="mcherry"
CELLRANGER_PATH="/usr/local/cellranger-4.0.0/cellranger"

################################################
## SCRIPT
################################################

[ -d $OUT_DIR ] || mkdir -p $OUT_DIR
cd $OUT_DIR

if [ ! -f "GRCh38.fa"] ; then
  # Human genome download and unzip
  wget http://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -c > GRCh38.fa
  gzip -d GRCh38.fa
fi
if [ ! -f "GRCh38.gtf"] ; then
  # Human annotation download and unzip
  wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
  gzip -d Homo_sapiens.GRCh38.101.gtf.gz -c > GRCh38.gtf
  gzip -d GRCh38.gtf
fi
if [ ! -f "MVA.fa"] ; then
  # Viral genome download and unzip
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Vaccinia_virus/latest_assembly_versions/GCF_000860085.1_ViralProj15241/GCF_000860085.1_ViralProj15241_genomic.fna.gz
  gzip -d GCF_000860085.1_ViralProj15241_genomic.fna.gz -c > MVA.fa
fi
if [ ! -f "MVA.gtf"] ;then
  # Viral annotation download andunzip
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Vaccinia_virus/latest_assembly_versions/GCF_000860085.1_ViralProj15241/GCF_000860085.1_ViralProj15241_genomic.gtf.gz
  gzip -d GCF_000860085.1_ViralProj15241_genomic.gtf.gz -c > MVA.gtf
fi

# Filter out unwanted information

echo "${genome1_name}"
filt_ann1="${genome1_name}-filtered.gtf"
ann2_exons="${genome2_name}-exons.gtf"
filt_ann2_exons="${genome2_name}-filtered_exons.gtf"


# Keep only protein coding genes from the human genome
$CELLRANGER_PATH mkgtf $annotation1 $filt_ann1 --attribute=gene_biotype:protein_coding

# Rename entries for viral genome (exons are required for CellRanger!!)
# And replace transcrip_id 'unknown_transcript' with gene_id
sed 's/CDS/exon/g' $annotation2 > $ann2_exons
perl -wpe 's/gene_id "([^"]+)"; transcript_id "[^"]+"/gene_id "$1"; transcript_id "$1"/;' $ann2_exons > $filt_ann2_exons

#### Merge #####

prefix="${genome1_name}_${genome2_name}_${genome3_name}"

genome="${prefix}_merged.fa"
filt_ann="${prefix}_merged-filtered_exons.gtf"

cat $genome1 $genome2 $genome3 > $genome

cat $filt_ann1 $filt_ann2_exons $annotation3 > $filt_ann
