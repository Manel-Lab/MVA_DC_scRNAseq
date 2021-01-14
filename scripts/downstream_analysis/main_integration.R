#!/usr/bin/Rscript
## ---------------------------
##
## This script will run CCA and MNN on conditions
## separately, or on the whole dataset
##
## Author: Kévin De Azevedo (kdeazevedo@github, kdeazeve@gitlab.curie.fr)
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-11-22

## ---------------------------

## Arguments
if (sys.nframe() == 0){
  
  .libPaths(c("~/R_packages",.libPaths()))
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 7) {
    print("Not enough arguments provided")
    print("Usage : Rscript main_integration.R batch jfile [WORKDIR]")
    q()
  } else if (length(args) > 7) {
    print("Too many arguments provided")
    print("Usage : Rscript main_integration.R batch jfile [WORKDIR]")  
    q()
  }
  
  jfile <- args[1]    # Path to json file
  path_i <- args[2]   # Path to the infected R object
  path_u <- args[3]   # Path to the uninfected R oject
  path_a <- args[4]   # Path to the merged R object
  dir_i <- args[5]    # Output directory for infected results
  dir_u <- args[6]    # Output directory for uninfected results
  dir_a <- args[7]    # Output directory for merged results
  
  ## Global variables
  
  # Read from json files
  result <- rjson::fromJSON(file = jfile)
}

## ---------------------------

## Package

library(Seurat)
source("utils/integration_methods.R")
source("utils/sc_utils.R")
source("utils/integration_eval.R")

## ------ FUNCTIONS ------- ##

# Take two Seurat objects 
# (two uninfected samples and two infected in our case)
# Compute common variable features
# Return the new objects in a list
# result = from the json arguments file
common_variableFeatures <- function(result, path_u, path_i){
  # Method for variable features
  method = result$integration$method
  # Number of variable features to compte
  nvf = result$integration$nvf
  # boolean, wether to intersect vf or not
  intersection = result$integration$intersect
  
  uninf <- get(load(path_u))
  inf <- get(load(path_i))
  uninf <- FindVariableFeatures(uninf, selection.method = method, nfeatures = nvf)
  inf <- humanVarFeatures(inf, method = method, nfeatures = nvf)

  # The variable features of the object
  # will be the intersection of those of the infected samples 
  # and those of the uninfected samples
  if(intersection){
    varfeatures <- intersect(VariableFeatures(inf), VariableFeatures(uninf))
    VariableFeatures(uninf) <- varfeatures
    VariableFeatures(inf) <- varfeatures
  }
  
  return(list(uninf, inf))
}

# Main function that computes CCA and MNN for an dataset
# robj = Seurat object
# samples = name of the samples
# outpath = output path
# result = from the json argument file
main_integration <- function(robj, samples, outpath, result){
  
  # Name of the batch effect, must be a Seurat metadata
  batch = result$integration$batch
  # Vector of different k neighbors to test
  k = result$integration$k
  # Vector of different components for MNN
  dmnn = result$integration$d
  # Number of components for CCA
  ndist = result$integration$ndist
  # Number of components for CCA
  cc = result$integration$cc
  # Arguments for integration 
  k.anchor = result$integration$k.anchor
  k.filter = result$integration$k.filter
  k.score = result$integration$k.score
  k.weight = result$integration$k.weight
  d1 = result$integration$d1
  d2 = result$integration$d2
  # Vector of reductions
  red = result$integration$red
  # Regression variables
  regressvars = result$integration$regressvars
  if(regressvars=="NULL"){regressvars<-NULL}
  
  if("mnn" %in% red){
    # Compute MNN
    dir.create(file.path(outpath))
    saveRDS(robj, file.path(outpath, "New_MNN.rds"))
    newRobj <- readRDS(file.path(outpath, "New_MNN.rds"))
    MNNapply(k = k, d = dmnn, ndist = ndist, robj = newRobj, batch = batch, path = outpath)
  }
  
  if("cca" %in% red){
    # Subset for RUN CCA
    dir.create(file.path(outpath))
    saveRDS(robj, file.path(outpath, "New_CCA.rds"))
    newRobj <- readRDS(file.path(outpath, "New_CCA.rds"))
    # VariableFeatures of each subset = VariableFeatures of the original dataset
    sub1 <- SplitObject(newRobj, split.by = batch)$`1`
    sub1 <- ScaleData(sub1, features = rownames(sub1), use.umi = TRUE, vars.to.regress = regressvars)
    VariableFeatures(sub1) <- VariableFeatures(newRobj)
    sub2 <- SplitObject(newRobj, split.by = batch)$`2`
    sub2 <- ScaleData(sub2, features = rownames(sub2), use.umi = TRUE, vars.to.regress = regressvars)
    VariableFeatures(sub2) <- VariableFeatures(newRobj)
    CCAapply(cc, sub1, sub2, batch, outpath)
  }
}

# Run integration successivelty on infected, uninfected, an all
# result : from the json argument file
integrate_all <- function(result, path_i, path_u, path_a,
                          dir_i, dir_u, dir_a){

  # Directory creation
  dir.create(dir_i)
  dir.create(dir_u)
  dir.create(dir_a)
  # Which analysis (one set of sample or all samples)
  analysis = result$integration$analysis
  # Take the intersection for varf (TRUE) 
  intersection = result$integration$intersect
  # Variable to regress for ScaleData (can be NULL)
  regressvars = result$integration$regressvars
  
  Arguments <- c(unlist(result$red), unlist(result$integration))
  
  
  # New objects with common variable features
  l <- common_variableFeatures(result, path_u, path_i)
  # Compute CCA and MNN with new objects
  if("uninfected" %in% analysis){
    obju <- l[[1]]
    print("Integration Uninfected samples")
    main_integration(obju, c("Uninf-rep1", "Uninf-rep2"), dir_u, result)
    print(paste0("Saved uninfected samples integration at ", dir_u))
    write.table(as.data.frame(Arguments), file.path(dir_u, "list_arguments_integration.tsv"), sep="\t")
  }
  if("infected" %in% analysis){
    obji <- l[[2]]
    print("Integration Infected samples")
    main_integration(obji, c("MVA-rep1", "MVA-rep2"), dir_i, result)
    print(paste0("Saved infected samples integration at ", dir_i))
    write.table(as.data.frame(Arguments), file.path(dir_i, "list_arguments_integration.tsv"), sep="\t")
  }
  if("both" %in% analysis){
    obja <- get(load(path_a))
    # Take variable features of the uninfected object
    # Depending of the value of the "intersection" argument,
    # Those are aither the computed variable features of all samples
    # Or the intersection of the vf of infected and uninfected
    VariableFeatures(obja) <- VariableFeatures(l[[1]])
    print("Integration All Samples samples")
    main_integration(obja, c(c("Uninf-rep1", "MVA-rep1"), c("Uninf-rep2", "MVA-rep2")), dir_a, result)
    print(paste0("Saved all samples integration at ", dir_a))
    write.table(as.data.frame(Arguments), file.path(dir_a, "list_arguments_integration.tsv"), sep="\t")
  }
}

## --------- RUN ---------- ## 

if (sys.nframe() == 0){
  integrate_all(result, path_i, path_u, path_a,
                dir_i, dir_u, dir_a)
}
