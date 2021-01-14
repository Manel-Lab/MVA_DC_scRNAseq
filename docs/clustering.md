# Clustering

## Rationale
These scripts will cluster a given reduction of a Seurat dataset with a varying number of parameters, and evaluate the validity of the clustering.

## Steps

### 1) Clustering

The script is ``main_clustering.R`` will create a number of clustering solution from a Seurat object with a dimensionality reduction (like MNN for instance). If parallelization is enabled, the script will spawns as many parallel processes as there are _k_ parameters, then, for each _k_, the different combination of parameters will be evaluated sequentially. The evaluation consist of computing the Silhouette score of the cluster, aswell as the Jaccard index of a bootstrapping for a clustering solution. Some arguments are passed through command line :

1. path to the Seurat Robj file with dimensionality reduction
2. path to the json arguments file
3. output path

Other parameters are defined as in ``args.json``, in the _clust_ section of the file :
- __clustcomp__ : boolean, wether to use parallization or not
- __k__ : list of number of neighbors for clustering
- __red__ : key of the reduction of interest (ex: _mnn1__)
- __res__ : list of resolutions for clustering
- __d__ : list of number of dimensions from the reduction to use : watch out to not exceed the number of dimensions actually computed
- __alg__ : list of clustering algorithm to use
- __prune__ : list of pruning thresholds
- __rand__ : list of random seeds
- __start__ : list of random starts
- __iter__ : list of number of iterations
- __n__ : number of iterations for bootstrapping
- __size__ : size of the sub samples during bootstrapping
- __template__ : template for the _future_ plan

The output of this script will be as many RDS file as there are k parameters, each containing several clustering solutions (combinations of this particular k and all other arguments). Each will correspond to a diagnostic tibble (summarizing the arguments used, Jaccard index and Silhouette score) and a PDF with the UMAP and the Silhouette profile. There will also be a a diagnostic tibble for all solutions, and a tsv with the list of arguments used. Each clustering solution will have a unique key (e.g. : _sol1017_)

### 2) Further evaluations

It is possible for a select solution to run a few more evaluations, namely a barplot depicting for each cluster the number of cells and the repartition of the two replicates, aswell as violin plots with features defined by the users (number of UMIs, genes etc.). The script for that is ``main_clust_eval.R``. Some arguments are passed through command line :

1. path to R RDS Seurat object with clustering solution of interest
2. arguments json file path
3. output path
