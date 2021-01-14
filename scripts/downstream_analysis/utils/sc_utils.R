#!/usr/bin/Rscript
## ---------------------------
##
## This script is comprised of utility functions that can be used
## For single cell projects 
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-10-25

## ---------------------------

## Packages

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tibble)
library(readr)
library(dplyr)
library(tidyr)

## ------ FUNCTIONS ------- ##

# Alternate normalization in case of presence of viral UMIs
# robj = Seurat object
# tsv = tsv annotating the genes (viral / non viral)
alt_norm <- function(robj, tsv) {
  # Read tsv
  isviral <- robj[["RNA"]]@meta.features["is_viral"]
  # Extract counts based on the nature of the gene (cellular or viral)
  virals <- subset(isviral, is_viral == 1)
  cellular <- subset(isviral, is_viral == 0)
  viralcounts <- robj[["RNA"]]@counts[rownames(robj[["RNA"]]@counts) %in% rownames(virals),]
  cellularcounts <- robj[["RNA"]]@counts[rownames(robj[["RNA"]]@counts) %in% rownames(cellular),]
  
  # Create a Seurat object with these counts
  # To use the Seurat normalization function separately
  viral_o <- robj
  cellular_o <- robj
  viral_o[["RNA"]] <- CreateAssayObject(viralcounts)
  cellular_o[["RNA"]] <- CreateAssayObject(cellularcounts)
  viral_o <- NormalizeData(viral_o)
  cellular_o <- NormalizeData(cellular_o)
  viraldata <- viral_o[["RNA"]]@data
  cellulardata <- cellular_o[["RNA"]]@data
  
  # New normalized data matrix
  newdata <- rbind(viraldata, cellulardata)
  return(newdata)
}

# Return a subset of a Seurat object, usually for testing purposes
# seurat_o = Seurat object
# n_cells = number of cells in the subset
# n_genes = number of genes in the subset
createTestSeurat <- function(seurat_o, n_cells, n_genes) {
  cells <- sample(colnames(seurat_o@assays$RNA@counts), n_cells, replace=FALSE)
  genes <- sample(rownames(seurat_o@assays$RNA@counts), n_genes, replace=FALSE)
  new_seurat = subset(seurat_o, cells = cells, features = genes)
  return(new_seurat)
}

# Add a meta.feature to Seurat object
# A boolean vector indicating if the gene is viral
# seurat_o = Seurat object
# csv = annotation csv
add_gene_nature <- function(seurat_o, csv) {
  is_viral <- read.csv2(csv)
  colnames(is_viral) <- c("gene", "is_viral_gene")
  is_viral$is_viral_gene <- 1
  tibble_rownames <- as_tibble(rownames(seurat_o[["RNA"]]@counts))
  colnames(tibble_rownames) <- c("gene")
  # Join viral genes and all genes in our object
  metadata <- full_join(is_viral, tibble_rownames)
  # If a gene is not in the csv, it is not a viral gene, else viral_gene = 0
  metadata$is_viral_gene <- replace_na(metadata$is_viral_gene, 0)
  my_df <- as.data.frame(metadata)
  rownames(my_df) <- metadata$gene
  # Annotate mcherry as a viral gene
  my_df[my_df$gene == "mcherry",]$is_viral_gene <- 1
  # Reorder names so they're in the same order at the count matrix
  my_df <- my_df[match(tibble_rownames$gene, my_df$gene),]
  seurat_o[["RNA"]] <- AddMetaData(seurat_o[["RNA"]], metadata = my_df$is_viral_gene, col.name = "is_viral")
  return(seurat_o)
}

# Add a meat.feature to Seurat object
# A string to describe the type of the viral gene
# The Seurat object must have been annotated with add_gene_nature !
# seurat_o = Seurat object
# tsv = annotation tsv
add_viral_class <- function(seurat_o, tsv) {
  classes <- read_tsv(tsv)
  all_genes <- tibble::as_tibble(seurat_o[["RNA"]]@meta.features["is_viral"]) %>%
    add_column(Gene = rownames(seurat_o[["RNA"]]@meta.features), .before=1)
  annot <- left_join(all_genes, classes) %>%
    mutate(Class = ifelse(is_viral == 0, replace_na(Class, "cellular"), replace_na(Class, "unknown")))
  seurat_o[["RNA"]] <- AddMetaData(seurat_o[["RNA"]], annot$Class, col.name = "gene_class")
  return(seurat_o)
}

# Take a list of ggplot of unknown lenght
# and create the grid and the PDF
# l = list of ggplot objects
# text = title
# pdf = output pdf path
gridtoPDF <- function(l, text, pdf){
  n <- length(l)
  print(n)
  nCol <- floor(sqrt(n))
  print(nCol)
  pdf(pdf)
  do.call("grid.arrange", c(l, ncol=nCol, top=text))
  dev.off()
}

# Find human variable features
# AND add them to the Seurat object
# robj = Seurat object
# method = method to find variable features (vst, mean.var.plot or disp)
# nfeatures = number of variable features to find
humanVarFeatures <- function(robj, method = "vst", nfeatures = 2000) {
  # Find top variable features and assign them back to the object
  robj <- FindVariableFeatures(robj, selection.method = method, nfeatures = nfeatures)
  # Replace variable features with the new ones
  one <- subset(robj[["RNA"]]@meta.features, is_viral == 0)
  two <- subset(one, vst.variable == TRUE)
  VariableFeatures(robj) <- rownames(two)
  return(robj)
}
