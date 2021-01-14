#!/usr/bin/Rscript
## ---------------------------
##
## Functions for MNN and CCA reduction
##
## Author: Kévin De Azevedo (kdeazevedo@github, kdeazeve@gitlab.curie.fr)
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-11-14

## ---------------------------

## Packages

library(Seurat)
library(SeuratWrappers)
library(foreach)
library(doFuture)
library("future.batchtools")
source("utils/integration_eval.R")
source("utils/clustering.R")

## ------ FUNCTIONS ------- ##

# Compute the MNN reduction with a select choice of parameters
# Output a Seurat object with MNN reduction and the
# UMAP of this reduction
# k = the number of neighbors
# d = number of dimensions to compute for the PCA
# ndist = kernel threshold
# robj = Seurat object
# batch = name of the Seurat metadata describing the batches
MNNparams <- function(k=20, d=50, ndist=3, robj, batch){

  # Wrapper function for fastMNN
  x <- RunFastMNN(SplitObject(robj, split.by = batch), 
                  k=k, d=d, ndist=ndist, 
                  # VariableFeatures must have already been computed
                  features = VariableFeatures(robj),
                  correct.all=TRUE
                  )
  # Run UMAP on this reduction
  print("run umap")
  x <- RunUMAP(x, reduction = "mnn", dims = 1:d)
  print("run umap done")
  return(x)
}

# Alternate way of computing MNN without wrapper
# (not used)
batchelorMNN <- function(k=20, d=50, ndist=3, robj, sub1, sub2, batch, path){
  
  mnn <- fastMNN(sub1[["RNA"]]@scale.data, sub2[["RNA"]]@scale.data, batch=batch, k = k, d = d, ndist = ndist,
                 subset.row	= VariableFeatures(sub1), cos.norm = FALSE)
  
  robj[["mnn"]] <- CreateDimReducObject(embeddings = SingleCellExperiment::reducedDims(mnn)$corrected, key = "mnn_")
  print("Finished MNN")
  robj <- RunUMAP(robj, "mnn", dims = 1:d)
  print("Finished UMAP")
  return(robj)
}

# Apply the MNNparams function with several successive values of 
# neighbor numbers, dimensions or distance cutoff
# k = list of int
# d = list of int
# ndist = list of int
# robj = Seurat R object
# batch = Seurat metadata with the batch
# path = output path
MNNapply <- function(k = c(20), d = c(20), ndist = c(3), 
                     robj, batch, path){
  print("Starting MNNapply")
  
  # Create all combinations
  args <- expand.grid(k = k, d = d, ndist = ndist)
  
  diag_tibble <- tibble(mnnkey = character(), k = numeric(),
                        d = numeric(), ndist = numeric(),
                        silhouette = numeric(), 
                        variance_lost_1 = numeric(),
                        variance_lost_2 = numeric())
  
  cnt = 1
  list_seurat <- list()
  for(i in seq_along(args$k)){
    # Save each Seurat object with new MNN reduction in a list
    list_seurat[[i]] <- MNNparams(args[i,1], args[i,2], args[i,3], robj, batch)
  }
  #for(i in seq_along(args$k)){
  #  sub1 <- SplitObject(robj, split.by = batch)$`1`
  #  sub1 <- ScaleData(sub1, features = rownames(sub1), model.use = "negbinom")
  #  VariableFeatures(sub1) <- VariableFeatures(robj)
  #  sub2 <- SplitObject(robj, split.by = batch)$`2`
  #  sub2 <- ScaleData(sub2, features = rownames(sub2), model.use = "negbinom")
  #  VariableFeatures(sub2) <- VariableFeatures(robj)
  #  batchvect <- robj@meta.data[[batch]]
  #  list_seurat[[i]] <- batchelorMNN(args[i,1], args[i,2], args[i,3], robj, sub1, sub2, batchvect, path)
  #}
  
  pdf(file.path(path, "MNN_results.pdf"))
  for(i in seq_along(list_seurat)){
    # For each object, add to the original object its
    # mnn and umap reduction with unique keys
    mnnkey <- paste0("mnn", i, "_")
    umapkey <- paste0("umap", i, "_")
    x <- list_seurat[[i]]
    robj[[mnnkey]] <- CreateDimReducObject(embeddings = Embeddings(x, reduction="mnn"), 
                                          key = mnnkey)
    robj[[umapkey]] <- CreateDimReducObject(embeddings = Embeddings(x, reduction="umap"), 
                                           key = umapkey)
    
    # For each MNN reduction
    for(j in seq_along(x@reductions)[c(TRUE, FALSE)]){
        Emb <- x@reductions[[j]]@cell.embeddings
        # Calculate the silhouette score
        batches <- as.factor(x@meta.data[[batch]])
        sil <- sil_estimate(Emb, batches)
        p1 <- UMAPplot(x, group = "orig.ident", title = umapkey, 
                       reduction = names(x@reductions[j+1]))
        p2 <- sil_plotting(sil, "Silhouette plot")
        g <- grid.arrange(p1, p2,
                          ncol=1, nrow=2)
        # Print the UMAP + Silhouette on the PDF
        print(g)
        # Create the diagnostic tibble
        diag_tibble[cnt, "mnnkey"] <- mnnkey
        diag_tibble[cnt, "silhouette"] <- round(mean(sil[,3]), digits = 2)
        diag_tibble[cnt, "variance_lost_1"] <- round(x@tools$RunFastMNN@metadata$merge.info$lost.var[1,1], digits = 2)
        diag_tibble[cnt, "variance_lost_2"] <- round(x@tools$RunFastMNN@metadata$merge.info$lost.var[1,2], digits = 2)
        cnt = cnt + 1
    }
  }
  dev.off()
  
  diag_tibble$k <- args$k
  diag_tibble$d <- args$d
  diag_tibble$ndist <- args$ndist
  # Save the final object with all MNN reductions, and the diagnostic table
  saveRDS(robj, file.path(path, "New_MNN.rds"))
  write_tsv(diag_tibble, file.path(path, "MNN_table.tsv"))
}


# Use RunCCA
CCAparams <- function(robj1, robj2, batch, cc=20, path){
  
  # RunCCA
  x <- RunCCA(robj1, robj2,
              num.cc = cc)
  x <- RunUMAP(x, reduction = "cca", dims = 1:cc)
  return(x)
}

# Apply the CCAparams function with several successive values of 
# neighbor numbers, dimensions or distance cutoff
CCAapply <- function(cc = c(20), sub1, sub2, batch, path){
  
  print("Starting CCAapply")
  
  diag_tibble <- tibble(ccakey = character(), cc = numeric(),
                        silhouette = numeric())
  
  cnt = 1
  # Will send a job on cluster for each interation
  list_seurat <- list()
  for(i in seq_along(cc)){
    list_seurat[[i]] <- CCAparams(sub1, sub2, batch, cc[i], path)
  }
  
  robj <- list_seurat[[1]]
  pdf(file.path(path, "CCA_results.pdf"))
  for(i in seq_along(list_seurat)){
    ccakey <- paste0("cca", i, "_")
    umapkey <- paste0("umap", i, "_")
    x <- list_seurat[[i]]
    robj[[ccakey]] <- CreateDimReducObject(embeddings = Embeddings(x, reduction="cca"), 
                                           loadings = Loadings(x, reduction="cca"),
                                           projected = Loadings(x, reduction="cca", projected=TRUE),
                                           key = ccakey)
    robj[[umapkey]] <- CreateDimReducObject(embeddings = Embeddings(x, reduction="umap"), 
                                            key = umapkey)
    
    for(j in seq_along(x@reductions)[c(TRUE, FALSE)]){
      Emb <- x@reductions[[j]]@cell.embeddings
      # Calculate the silhouette score
      batches <- as.factor(x@meta.data[[batch]])
      sil <- sil_estimate(Emb, batches)
      p1 <- UMAPplot(x, group = "orig.ident", title = umapkey, 
                     reduction = names(x@reductions[j+1]))
      p2 <- sil_plotting(sil, "Silhouette plot")
      g <- grid.arrange(p1, p2,
                        ncol=1, nrow=2, top = names(x@reductions[j+1]))
      print(g)
      
      diag_tibble[cnt, "ccakey"] <- ccakey
      diag_tibble[cnt, "silhouette"] <- round(mean(sil[,3]), digits = 2)
      cnt = cnt + 1
    }
  }
  dev.off()
  
  diag_tibble$cc <- cc
  
  saveRDS(robj, file.path(path, "New_CCA.rds"))
  write_tsv(diag_tibble, file.path(path, "CCA_table.tsv"))
}
