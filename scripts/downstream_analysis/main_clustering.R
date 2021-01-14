## ---------------------------
##
## Purpose of script: Clustering a Seurat dataset
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2019-12-12
##
## ---------------------------



## Arguments

if (sys.nframe() == 0){
  .libPaths(c("~/R_packages",.libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    print("Not enough arguments provided")
    print("Usage : Rscript robj jfile <path> <WORKDIR>")
    q()
  } else if (length(args) > 4) {
    print("Too many arguments provided")
    print("Usage : Rscript robj jfile <path> <WORKDIR>")  
    q()
  }
  
  robj <- args[1]     # Seurat R Object Path
  jfile <- args[2]    # Path to json file
  path_to_cluster <- args[3]     # Output path

  
## Global variables

# Read from json files
result <- rjson::fromJSON(file = jfile)
Arguments <- c(unlist(result$clust), unlist(result$integration))
}

## ---------------------------

## Packages

library(Seurat)
library(ggplot2)
library(gridExtra)
library(clusteval)
library(foreach)
library(doFuture)
library("future.batchtools")
source("utils/clustering.R")

## ------ FUNCTIONS ------- ##

# This function takes a list with Seurat objects 
# (all with different clustering solutions)
# and create one Seurat object with all solutions
# path_to_cluster = output path
# l = list containing Seurat object and arguments
customClusters <- function(path_to_cluster, l){
    args <- l[[1]]
    # Name of the new reduction : red + arguments
    allargs <- paste(paste(names(args), args, sep = "_"), collapse = "_")
    # Create the new reduction
    print(path_to_cluster)
    print(args$k)
    new <- readRDS(file.path(path_to_cluster, paste0(args$k, "k_Clusters_All-samples.rds")))
    newkey <- paste0("sol", args$k, dim(new@meta.data)[2])
    new <- AddMetaData(new, l[[8]]$seurat_clusters, col.name = newkey)
    # Save new Seurat object with all reductions
    saveRDS(new, file.path(path_to_cluster, paste0(args$k, "k_Clusters_All-samples.rds")))
    return(newkey)
}

# Run clustering, Silhouette and save in list
# robj = Seurat object
# path_to_cluster = output path
# n = number of iterations for bootstrapping
# size = size of the sub-samples
# other args  : clustering arguments
RunClustering <- function(robj, path_to_cluster, k = 30, red = "umap", d = 20, prune = 1/15,
                          res = 1, alg = 1, rand = 0, start = 10, iter = 10,
                          n = 10, size = 500){
  
  msg = paste0("Clustering, d = ", d, ", k = ", k, ", prune = ", prune, ", res = ", res, ", alg = ", alg,
               ", rand = ", rand, ", start = ", start, ", iter =  ", iter)
  print(msg)
  
  # Clustering
  robj <- Neighborparams(robj, k = k, red = red, d = d, prune = prune)
  robj <- Clusterparams(robj, res = res, alg = alg, rand = rand, start = start, iter = iter)
  dimplot <- UMAPplot(robj, reduction = stringr::str_replace(red, "mnn", "umap"), 
                      group = "orig.ident", title = "UMAP plot")
  clusterplot <- UMAPplot(robj, reduction = stringr::str_replace(red, "mnn", "umap"), 
                          group = "seurat_clusters", title = "UMAP plot")

  # Jaccard
  jaccard <- BootstrapClustering(robj, n, size,
                                 k = k, red = red, d = d, prune = prune,
                                 res = res, alg = alg, rand = rand, start = start, iter = iter)
  
  # Silhouette : don't compute if number of clusters is 1
  if(length(levels(robj$seurat_clusters)) != 1){
    scoreSil <- computeSilhouette(robj, red)
    plotSil <- plotSilhouette(scoreSil)
  }
  else{
    scoreSil <- NA
    plotSil <- NA
  }
  
  argslist <- list(k, red, d, prune, res, alg, rand, start, iter)
  names(argslist) <- c("k", "red", "d", "prune", "res", "alg", "rand", "start", "iter")
  l <- list(argslist, dimplot, clusterplot, scoreSil, plotSil, jaccard, msg, robj)
  newkey <- customClusters(path_to_cluster, l)
  print("Removing robj")
  l[8] <- newkey
  return(l)
}

# This functions creates several sub-samples of our Seurat object 
# Run the same clustering, and compute the Jaccard index between them 
# To test the robustness of the solution
# robj = Seurat object
# n = number of iterations
# size = size of subsamples
# other args = clustering arguments
BootstrapClustering <- function(robj, n, size,
                                k, red, d, prune,
                                res, alg, rand, start, iter){
  print("Jaccard index")
  # Create a list of Seurat objects with n cells
  v <- list()
  for(i in 1:n) {
    # Choose cells by draw without replacement
    cells <- sample(names(robj@active.ident), size, replace = FALSE)
    # Create a subset Seurat object based on these cells
    s <- SubsetData(robj, cells = cells)
    s <- Neighborparams(s, k = k, red = red, d = d, prune = prune)
    s <- Clusterparams(s, res = res, alg = alg, rand = rand, start = start, iter = iter)
    v[[i]] <- s
  }
  
  scores <- c()
  # Calculate the Jaccard index between each of the solutions
  for(i in seq_along(v)){
    for(j in i:length(v)){
      s <- cluster_similarity(v[[i]]$seurat_clusters, v[[j]]$seurat_clusters, 
                              similarity = c("jaccard"), method = "independence")
      scores <- append(scores, s, after = length(scores))
    }
  }
  rm(v)
  return(c(mean(scores), var(scores), min(scores), max(scores)))
}

# Apply the clustering for several parameter combinations
# robj = Seurat object
# path_to_cluster = output path
# my_vector = contain all clustering solutions
# other args = list of values for clustering parameters
ClusteringApply2 <- function(k = c(30), red = "umap", d = c(20), prune = c(1/15),
                            res = c(1), alg = c(1), rand = c(0), start = c(10), iter = c(10),
                            n = 10, size = 500, robj, path_to_cluster, my_vector){
  print("Starting ClusteringApply")
  
  # Create all combinations of argyments
  args <- expand.grid(k = k, d = d, prune = prune, res = res, alg = alg,
                      rand = rand, start = start, iter = iter)
  
  # Will send a job on cluster for each interation
  for(i in 1:length(args$k)){ 
    # Create the filename
    filename <- file.path(path_to_cluster,
                          paste("Clustering", args[i,1], red, args[i,2], args[i,3],
                                args[i,4], args[i,5], args[i,6], args[i,7], args[i,8],
                                ".rds", sep="_"))
    my_list <- RunClustering(robj, path_to_cluster, args[i,1], red, args[i,2], args[i,3], args[i,4], 
                             args[i,5], args[i,6], args[i,7], args[i,8], n, size)
    print("Saving to vector")
    my_vector[[i]] <- my_list
    print("Saved")
  }
  print("Returning vector")
  return(my_vector)
}

# Diagnostic tibble outputing in table format 
# the Jaccard index, the Silhouette score, and the args used
# clustvect = list of clustering solutions
# outpath = output path of the tibble
clustering_tibble <- function(clustvect, outpath){
  print("Creating diagnostic tibble")
  
  silhouette <- c()
  mean_jaccard <- c()
  var_jaccard <- c()
  min_jaccard <- c()
  max_jaccard <- c()
  k <- c()
  d <- c()
  prune <- c()
  res <- c()
  alg <- c()
  rand <- c()
  start <- c()
  iter <- c()
  key <- c()
  
  # Iterate through the list of Seurat object
  for(i in seq_along(clustvect)){
    obj <- clustvect[[i]]
    # If the number of cluster is 1, there is no Silhouette score
    if(!(is.na(obj[[4]]))){
      silhouette[i] <- mean(obj[[4]][,3])
    }
    else{
      silhouette[i] <- c(NA)
    }
    key[i] <- obj[[8]]
    k[i] <- obj[[1]]$k
    d[i] <- obj[[1]]$d
    prune[i] <- obj[[1]]$prune
    res[i] <- obj[[1]]$res
    alg[i] <- obj[[1]]$alg
    rand[i] <- obj[[1]]$rand
    start[i] <- obj[[1]]$start
    iter[i] <- obj[[1]]$iter
    mean_jaccard[i] <- obj[[6]][1]
    var_jaccard[i] <- obj[[6]][2]
    min_jaccard[i] <- obj[[6]][3]
    max_jaccard[i] <- obj[[6]][4]
  }
  
  print("Writing tibble")
  my_tibble <- tibble(key=key, k=k, d=d, prune=prune, res=res, alg=alg,
                      rand=rand, start=start, iter=iter,
                      mean_jaccard=mean_jaccard, var_jaccard=var_jaccard,
                      min_jaccard=min_jaccard, max_jaccard=max_jaccard,
                      silhouette=silhouette)
  # Write a tibble for this particular k
  write_tsv(my_tibble, file.path(outpath, paste0(obj[[1]]$k, "k_clustering_diag_table.tsv")))
  return(my_tibble)
}

# Create the recap PDF for a particular k
# clustvect = vector of clustering solutions
recapPDF <- function(clustvect){
  print("Doing PDF")
  
  # Iterate over list of Seurat objects for this k
  for(i in seq_along(clustvect)){
    obj <- clustvect[[i]]
    # Create grid
    # The grid is different if we have a Silhouette plot or not
    if(!(is.na(obj[[4]]))){
      g <- grid.arrange(obj[[2]], obj[[3]], obj[[5]],
                 ncol=1, nrow=3, top=obj[[7]])
      print(g)
      rm(g)
    }
    else{
      g <- grid.arrange(obj[[2]], obj[[3]],
                   ncol=1, nrow=2, top=obj[[7]])
      print(g)
      rm(g)
    }
  }
}

# Main function for doing clustering with parallel future
# robj = Seurat object
# path_to_cluster = output path
# other args = list of values for clustering parameters
do_clustering_parallel <- function(robj, k, red, d, prune,
                            res, alg, rand, start, iter,
                            n, size, path_to_cluster){
  
  # Cluter template : either a path to a template
  # or an predefined template from future
  template = result$clust$template
  
  # Create output directory if it doesn't already exists
  dir.create(file.path(path_to_cluster))
  
  registerDoFuture()
  options(future.globals.maxSize = 1000*1024^2)
  #TODO: TEMPLATE AS ARGUMENT !
  plan(batchtools_torque, template=template, resources=list(ppn=1, memory=5000, walltime=1800,
                                                                   queue="q_small"))

  list_tibble <- list()
  # Create a PDF, a tibble, and an object different for every k parameter
  list_tibble <- foreach(i = 1:length(k), .export=c("list_tibble")) %dopar% {
    
    my_vector <- list()
    saveRDS(robj, file.path(path_to_cluster, paste0(k[i], "k_Clusters_All-samples.rds")))
    newRobj <- readRDS(file.path(path_to_cluster, paste0(k[i], "k_Clusters_All-samples.rds")))
    my_vector <- ClusteringApply2(k = k[i], red = red, d = d, prune = prune,
                                 res = res, alg = alg, rand = rand, start = start, iter = iter,
                                 n = n, size = size, robj, path_to_cluster, my_vector)
    pdf(file.path(path_to_cluster, paste0(k[i], "k_Clusters_results.pdf")))
    recapPDF(my_vector)
    dev.off()
    list_tibble <- clustering_tibble(my_vector, path_to_cluster)
    return(list_tibble)
  }
  # Write all diagnostic tibbles in 1
  final_tibble <- dplyr::bind_rows(list_tibble)
  write_tsv(final_tibble, file.path(path_to_cluster, "Clustering_diag_table.tsv"))
  # Write the table with arguments used for this run at the end
  #write.table(as.data.frame(Arguments), file.path(path_to_cluster, "list_arguments_clustering.tsv"), sep="\t")
}

# Main function for doing clustering sequentially
# robj = Seurat object
# path_to_cluster = output path
# other args = list of values for clustering parameters
do_clustering_seq <- function(robj, k, red, d, prune,
                              res, alg, rand, start, iter,
                              n, size, path_to_cluster){
  
  print("Entering do clustering local")
  # Create output directory if it doesn't already exists
  dir.create(file.path(path_to_cluster))
  list_tibble <- list()
  # Create a PDF, a tibble, and an object different for every k parameter
  for(i in 1:length(k)){
    
    print("Entering loop")
    my_vector <- list()
    print(path_to_cluster)
    saveRDS(robj, file.path(path_to_cluster, paste0(k[i], "k_Clusters_All-samples.rds")))
    newRobj <- readRDS(file.path(path_to_cluster, paste0(k[i], "k_Clusters_All-samples.rds")))
    my_vector <- ClusteringApply2(k = k[i], red = red, d = d, prune = prune,
                                  res = res, alg = alg, rand = rand, start = start, iter = iter,
                                  n = n, size = size, robj, path_to_cluster, my_vector)
    pdf(file.path(path_to_cluster, paste0(k[i], "k_Clusters_results.pdf")))
    recapPDF(my_vector)
    dev.off()
    list_tibble <- append(list_tibble, clustering_tibble(my_vector, path_to_cluster))
  }
  # Write all diagnostic tibbles in 1
  final_tibble <- dplyr::bind_rows(list_tibble)
  write_tsv(final_tibble, file.path(path_to_cluster, "Clustering_diag_table.tsv"))
  # Write the table with arguments used for this run at the end
  #write.table(as.data.frame(Arguments), file.path(path_to_cluster, "list_arguments_clustering.tsv"), sep="\t")
}

# Main clustering function
# robj = Seurat object
# path_to_cluster = output path
# result = from json args file
main_clustering <- function(robj, path_to_cluster, result){
  
  # Vector of number of neighbors
  k = result$clust$k
  # Name of the reduction
  red = result$clust$red
  # Vector of dimensions
  d = result$clust$d
  # Vector of pruning threshold 
  prune = result$clust$prune
  # Vector of number of resolutions
  res = result$clust$res
  # Vector of algorithms
  alg = result$clust$alg
  # Vector of seeds
  rand = result$clust$rand
  # Vector of random starts
  start = result$clust$start
  # Vector of number of iterations
  iter = result$clust$iter
  # Number of bootstrap to do
  n = result$clust$n
  # Size of bootstrap sub-sample
  size = result$clust$size
  # Parallel or no
  clustcomp = result$clust$clustcomp
  
  robj <- readRDS(robj)
  if(clustcomp){
    do_clustering_parallel(robj, k, red, d, prune,
                  res, alg, rand, start, iter,
                  n, size, path_to_cluster)
  }
  else{
    do_clustering_seq(robj, k, red, d, prune,
                  res, alg, rand, start, iter,
                  n, size, path_to_cluster)
  }
}


## --------- RUN ---------- ##

if (sys.nframe() == 0){
  main_clustering(robj, path, result)
}