# MVA-DCs single cell RNA-seq

## Introduction

This repository contains bash and R scripts analyzing single cell aligned data from dentritic cells infected by modified vaccinia virus and healthy matched controls. It is comprised of four steps.

1. [Downloading the genome, annotations, and concatenating them](docs/concat.md)
2. [Creating reference transcriptome with cellRanger](docs/mkref.md)
3. [Alignment and feature count with cellRanger](docs/alignment.md)
4. [Downstream analysis with R](docs/mva_dc_downstream.md)
  * [Cell calling and annotation of Seurat object](docs/object_creation.md)
  * [Donor effect removal](docs/removing_donor_effect.md)
  * [Clustering](docs/clustering.md)
  * [Differential expression](docs/diff_exp.md)
  * [Paper figures](docs/paper_figures.md)

The four steps are executed separately thanks to a main bash script for each of them. These bash scripts can be executed locally but are also configured to work in a Torque cluster configuration. For any other configuration, the user will have to modify these main scripts.

## Singularity usage

A Singularity container is provided, with every tool pre-installed to run the pipeline, and it is recommended for reproducibility. You can get the container from Sylab :

``singularity pull library://kdeazevedo/default/mva_dc_scrnaseq``

The Singularity version used was 3.5.3. If you wish, you can also build the Singularity container yourself with the Singularity file provided, in ``conf/Singularity_file.txt`` (superuser rights necessary) :

``sudo https_proxy=$https_proxy singularity build --sandbox My_container Singularity.KD01.txt``

You can then enter the singularity container to launch the scripts the same way you would on any machine :

``singularity exec --no-home My_container bash``

You can mount a local directory to the container as well, for instance to mount the input data. The container comes with two mount points called /mnt and /project :

``singularity exec --bind [local/dir]:/project --no-home My_container bash``

For more information on installing and using Singularity, see [the official website](https://sylabs.io/guides/3.5/user-guide/quick_start.html).
