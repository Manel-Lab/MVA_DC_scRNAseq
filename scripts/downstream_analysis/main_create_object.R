#!/usr/bin/Rscript
## ---------------------------
##
## This script creates the R object used in the analysis
##
## Author: Kévin De Azevedo (kdeazevedo@github, kdeazeve@gitlab.curie.fr)
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-10-24

## ---------------------------

if (sys.nframe() == 0){
  .libPaths(c("~/R_packages",.libPaths()))
  
  ## Arguments
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 8) {
    print("Not enough arguments provided")
    print("Usage : Rscript inf uninf1 uninf2 mva_tsv class_tsv normalize path_out")
    q()
  } else if (length(args) > 8) {
    print("Too many arguments provided")
    print("Usage : Rscript inf uninf1 uninf2 mva_tsv class_tsv normalize path_out")  
    q()
  }
  
  inf1 <- args[1]          # Path of Seurat Robj of merged infected
  inf2 <- args[2]          # Path of Seurat Robj of merged infected
  uninf1 <- args[3]       # Path of Seurat Robj of uninfected 1
  uninf2 <- args[4]       # Path of Seurat Robj of uninfected 2
  mva_csv <- args[5]      # Path to tsv detailing MVA genes
  class_tsv <- args[6]    # Path to tsv annotating MVA genes
  normalize <- as.logical(args[7])    # Wether or not to use the alt_norm function
  path_out <- args[8]     # Path of the outputed R object
}

## ---------------------------

## Packages

library(Seurat)
source("utils/sc_utils.R")

## ------ FUNCTIONS ------- ##

# Create a infected object from two separate replicates
# path1 = path to the first object
# path2 = path to the second object
create_inf <- function(path1, path2){
  inf1$orig.ident <- "MVA-rep1"
  inf2$orig.ident <- "MVA-rep2"
  # Add the donor ID 
  inf1 <- AddMetaData(inf1, col.name = "cell_id_sil", metadata = 1)
  inf2 <- AddMetaData(inf2, col.name = "cell_id_sil", metadata = 2)
  inf <- merge(inf1, inf2, add.cell.ids = c("MVA-rep1", "MVA-rep2"))
  Idents(inf) <- inf$orig.ident
  return(inf)
}

# Create a Uninfected object from two separate replicates
# path1 = path to the first object
# path2 = path to the second object
create_uninf <- function(path1, path2){
  uninf1$orig.ident <- "Uninf-1"
  uninf2$orig.ident <- "Uninf-2"
  # Add the donor ID 
  uninf1 <- AddMetaData(uninf1, col.name = "cell_id_sil", metadata = 1)
  uninf2 <- AddMetaData(uninf2, col.name = "cell_id_sil", metadata = 2)
  uninf <- merge(uninf1, uninf2, add.cell.ids = c("Uninf-rep1", "Uninf-rep2"))
  Idents(uninf) <- uninf$orig.ident
  return(uninf)
}

# Main function for creating the objects
# inf = path to the uninfected object
# uninf1 = path to the uninfected replicate 1
# uninf2 = path to the uninfected replicate 2
# mva_tsv = path to the tsv annotating the genes (viral/non viral)
# class_tsv = path to the tsv with the viral gene classes (early, late etc.)
main_create_object <- function(inf1, inf2, uninf1, uninf2, mva_csv, class_tsv, path_out){
  
  # Create the output path directory
  dir.create(path_out)
  
  # Infected
  print("Creating MVA R object")
  if(normalize){filename <- "MVA_sepnorm.Robj"}else{filename <- "MVA_seuratnorm.Robj"}
  # Add the nature of the gene (viral or not)
  infected <- create_inf(inf1, inf2)
  infected <- add_gene_nature(infected, mva_csv)
  infected <- add_viral_class(infected, class_tsv)
  # Use the separate normalization if specified
  if(normalize){infected[["RNA"]]@data <- alt_norm(infected, mva_csv)}
  save(infected, "infected", file = file.path(path_out, filename))
  print(paste0("MVA R object saved at ", file.path(path_out, filename)))
  
  # For Uninfected, need to be merged first
  print("Creating Uninfected R object")
  if(normalize){filename <- "Uninfected_sepnorm.Robj"}else{filename <- "Uninfected_seuratnorm.Robj"}
  uninfected <- create_uninf(uninf1, uninf2)
  uninfected <- add_gene_nature(uninfected, mva_csv)
  uninfected <- add_viral_class(uninfected, class_tsv)
  if(normalize){uninfected[["RNA"]]@data <- alt_norm(uninfected, mva_csv)}
  save(uninfected, "uninfected", file = file.path(path_out, filename))
  print(paste0("Uninfected R object saved at ", file.path(path_out, filename)))
  
  # Create one with all samples
  print("Creating R object with all samples")
  if(normalize){filename <- "All-samples_sepnorm.Robj"}else{filename <- "All-samples_seuratnorm.Robj"}
  all <- merge(infected, uninfected)
  all <- add_gene_nature(all, mva_csv)
  all <- add_viral_class(all, class_tsv)
  if(normalize){all[["RNA"]]@data <- alt_norm(all, mva_csv)}
  save(all, "all", file = file.path(path_out, filename))
  print(paste0("All Samples R object saved at ", file.path(path_out, filename)))

}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  inf1_obj <- get(load(inf1))
  inf2_obj <- get(load(inf2))
  uninf1_obj <- get(load(uninf1))
  uninf2_obj <- get(load(uninf2))
  main_create_object(inf1_obj, inf2_obj, uninf1_obj, uninf2_obj, 
                     mva_tsv, class_tsv)
}