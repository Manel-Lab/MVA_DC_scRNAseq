# Downstream analysis

## Introduction

This downstream analysis is comprised of several steps that are usually better executed separately, as they provide evaluation metrics alongside their run which make the user able to know wether to continue on with the analysis or rethink their design plan. However, for reproducibility purposes, a R script called ``mvadcs_downstream.R`` is available to re-run the full pipeline from start to finish with parameters used for the published results. It takes as arguments :

1. The input path where the count tables of CellRanger count are
2. A csv mapping the official gene name in the MVA genome to a "common" gene name
3. A tsv indicating the class of MVA gene (early, intetmediary, late).
4. A boolean indicating wether to normalize features separately depending of their nature
5. The output path
6. The path to the json containing most other arguments.
7. Optional, the working directory, otherwise it is the current directory. Needs to be the directory where the other scripts are, since they are imported in the ``mvadcs_downstream.R`` file.

## Args json file

Each part of the pipeline fetch most of its arguments from a json file. The same json file is shared through all parts, the structure of the file separating arguments from one part from the others. For more informations on what the arguments mean, refer to the proper README.

## Usage

You can execute it as a R script, from inside or oustidethe container :

``Rscript --vanilla mvadcs_downstream.R ../data/MVA_gene_map.csv ../data/mva-class.tsv TRUE results/Rscript ../args.json``

Or execute it from outside the container with the ``Rscript_singularity.sh`` file :

``bash Rscript_singularity.sh /abolute/path/to/project``
