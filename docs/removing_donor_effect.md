# Donor effect removal

## Rationale
The goal of this R pipeline is to compare several integration methods to try and determine the method that remove the variability between our two donors without removing too much of the biological variability.

## Steps

### 1) Donor effect removal

The first script is ``main_integration.R`` runs MNN and/or CCA on the dataset with varying parameters.  Some paramaters for the script are passed directly through command line :
1. relative or absolute path to json file
2. The path to the infected Seurat object
3. The path to the uninfected Seurat object
4. The path to the whole dataset Seurat object
5. Writing directory for infected results
6. Writing directory for uninfected results
7. Writing directory for whole dataset results

Other parameters are defined as in a json file called ``args.json``. In this json, the _"red"_ section correspond to arguments relating to the reduction methods themsleves, while the _"integration"_ section correspond to arguments of the pipeline itself. The arguments are :
- __batch__ : Seurat metadata of the batch effect (ex : metadata specifying the donor)
- __k__ : A list with the different neighbors to test for MNN
- __d__ : A list with the different dimensions to test for MNN
- __cc__ : A list with the different dimensions to test for CCA
- __red__ : A list of reductions to compute, i.e. ["mnn"] or ["mnn", "cca"] for instance
- __method__ : How to compute variable features (vst, disp  or mean.var.plot)
- __nvf__ : Number of variable features to compute (for vst and disp methods).
- __regressvars__ : List of variable to regress while scaling the data. If no regression is needed, it is "NULL"
- __analysis__ : What object(s) should the script produce results for (valid values are _infected_, _uninfected_ and _both_)
- __intersect__ : A boolean indicating wether the variable features computed should be the intersection of those of the infected and uninfected datasets ("TRUE") or just those computed on the dataset ("FALSE")

The outputs, written in the relevant output paths, will be :
1. A new Seurat object in RDS format with as many new reduction slots as there were combinations of parameters for the integration method. The keys of these reduction slots will be unique (ex : _mnn1_). To know which key correspond to which parameters, see point 3.
2. A PDF with each page corresponding to a set of parameters, and containing the UMAP of the MNN with the samples and the Silhouette plot.
3. A diagnostic table summarizing for each set of parameter the Silhouette score, and (for MNN) the percentage of variance lost at each merge for each sample.
4. A tsv summarizing the different arguments used for the run

Those informations should be sufficient to determine what are the best parameters for CCA and MNN, but to compare the two, it might be necessary to see how they do against clustering.

### 2) Testing solution against clustering

The second script ``main_cluster_red_solution.R`` will do a clustering of the dataset based on the reduction, with default parameters and varying levels of resolutions. Then, the Silhouette score of the solutions varying in cluster numbers will be computed, put in another list, as well as the corresponding Silhouette plot. Then, several other diagnostic plots will be computed. Parameters for the script that are passed through command line are :

1. Path to the RDS of a Seurat object
2. Output path
3. A json argument file

Arguments in the json are :

- __red__ : Key of the reduction of interest (e.g. _mnn1__)
- __list_genes__ : list of genes of interests to plot
- __res__ : list of resolutions to test for the clustering
- __d__ : number of dimensions from the dimension reduction to use
- __is_norm__ : boolean, wether to use normalized data
- __is_log__ : boolean, wether to log transform the number of UMIs before plotting
- __threshold__ : threshold value of detection for a gene : genes below this value will be considered "Absent". This only has effect in a particular plot showing the viral UMi aganst the cellular UMI in each cell.
- __clustering__ : boolean, wether to do the clustering part
- __geneplot__ : boolean, wether to output the plots

The different outputs will be in the path defined in command line.

The results of the clustering as well as the plots should be enough to decide which methods to use.
