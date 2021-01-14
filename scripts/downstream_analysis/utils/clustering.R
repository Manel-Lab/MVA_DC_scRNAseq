#!/usr/bin/Rscript
## ---------------------------
##
## Seurat functions regarding clustering
## Including UMAP, and DimPlots for clustering
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Based on Tarek Gharsalli's scripts
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-11-14

## ---------------------------

## Packages

library(Seurat)
library(ggplot2)
source("utils/integration_eval.R")

## ------ FUNCTIONS ------- ##

# Run UMAP reduction
UMAPparams <- function(robj, red = "pca", d = 20, redname = "umap"){
  robj <- RunUMAP(robj, reduction = red, dims = 1:d, reduction.name	= redname)
  return(robj)
}

# Find number of neighbors for each cell
Neighborparams <- function(robj, k = 30, red = "umap", d = 20, prune = 1/15){
  robj <- FindNeighbors(robj, k.param = k, force.recalc = TRUE, 
                        reduction = red, dims = 1:d, prune.SNN = 0.07)
  return(robj)
}

# Find clusters for our data
Clusterparams <- function(robj, res = 1, alg = 1, rand = 0, start = 10, iter = 10){
  robj <- FindClusters(robj, resolution = res, algorithm = alg,
                       random.seed = rand, n.start = start, n.iter = iter)
  return(robj)
}

# Plot UMAP reduction with DimPlot
UMAPplot <- function(robj, reduction = "umap", group = NULL, shape = NULL, title = "UMAP plot"){
  d = DimPlot(robj, reduction = reduction, group.by = group, shape.by = shape) +
    labs(title = title)
  return(d)
}

# Compute Silhouette score
computeSilhouette <- function(object, red){
  Emb <- Embeddings(object, reduction = red)
  mat <- as.matrix(Emb)
  clust <- as.vector(as.numeric(object@meta.data$seurat_clusters))
  silscore <- sil_estimate(mat, clust)
  return(silscore)
}

# Plot Silhouette score
plotSilhouette <- function(sil_score){
  silplot <- sil_plotting(sil_score, "Silhouette")
  return(silplot)
}
