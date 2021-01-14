#!/usr/bin/Rscript
## ---------------------------
##
## Create various diagnostic plots
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-11-19

## ---------------------------

## Arguments

if (sys.nframe() == 0){
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    print("Not enough arguments provided")
    print("Usage : Rscript various_plot.R seurat_o [WORKDIR]")
    q()
  } else if (length(args) > 2) {
    print("Too many arguments provided")
    print("Usage : Rscript integration_eval.R seurat_o [WORKDIR]")  
    q()
  }
  
  seurat_o <- args[1]  # Path to Seurat object
}

## ---------------------------

## Packages

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tibble)
library(dplyr)

## ------ FUNCTIONS ------- ##

# Creates a tibble with interesting informations
# seurat_o : A Seurat object
# is_norm : boolean, use normalized value ?
# is_log : boolean, use log-transformed value ?
# threshold : min raw nUMI one's consider as a threshold for detection
create_tibble <- function(seurat_o, list_genes, is_norm = FALSE, is_log = TRUE, threshold = 10){
  
  # Use normalized value or not
  if(is_norm){
    my_counts <- seurat_o@assays$RNA@data
  }
  else{
    my_counts <- seurat_o@assays$RNA@counts
  }
  
  # Dataframe
  counts <- as_tibble(t(as.matrix(my_counts)))
  counts <- dplyr::select(counts, list_genes)
  df <- tibble(nUMI_MVA = seurat_o$nUMI_MVA, 
               nUMI_Cellular = seurat_o$nUMI_cellular,
               nCount_RNA = seurat_o$nCount_RNA)
  df <- bind_cols(df, counts)

  # Add cells, donor, and presence/absence boolean
  tibble = df %>%
    add_column("Cell" = colnames(my_counts), .before=1) %>%
    add_column("orig.ident" = seurat_o$orig.ident) %>%
    mutate_if(is.numeric, funs(`is` = ifelse(.> threshold, "Presence", "Absence"))) %>%
    # Create a new column sorting a cell in our free observed groups
    mutate(group = ifelse(log10(nUMI_MVA) > 2 & log10(nUMI_Cellular) < 3, "Dying",
                          ifelse(log10(nUMI_MVA) > 2 & log10(nUMI_Cellular) > 3, "Resist", "Normal")))
  
  # Log transform the tibble
  if(is_log){
    tibble <- tibble %>% mutate_at(c("nUMI_MVA", "nUMI_Cellular"), ~log10(.+1))
  }
  
  # Add CCA and MNN reductions if they exists
  if(!(is.null(seurat_o@reductions$cca))){
    tibble <- as_tibble(cbind(tibble, seurat_o@reductions$cca@cell.embeddings))
  }
  if(!(is.null(seurat_o@reductions$mnn))){
    tibble <- as_tibble(cbind(tibble, seurat_o@reductions$mnn@cell.embeddings))
  }
  if(!(is.null(seurat_o@reductions$umap))){
    tibble <- as_tibble(cbind(tibble, seurat_o@reductions$umap@cell.embeddings))
  }
  
  return(tibble)
  
}

# Create plot of cellular UMI and viral UMI
# tibble : relevant informations
# title : title of the plot
create_nUMI_plot <- function(tibble, title){
  p <- ggplot(tibble, mapping=aes(x = nUMI_Cellular, y= nUMI_MVA)) + geom_point() +
    labs(title = title, x = "Cellular UMI count", y = "Viral UMI count") +
    scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10')
  return(p)
}

# Plot to vizualise the expression or absence/presence of a gene
# tibble : relevant informations
# title : title of the plot
# gene : name of the gene to vizualise
# xlab : title of x axis
# ylab : title of y axis
# shape : column to shape points by (optional)
# x : x axis (nUMI_cellular or a dimensionality reduction)
# y : y axis (nUMI_MVA or a dimensionality reduction)
create_genexpr_plot <- function(tibble, title, gene,
                                x = "nUMI_Cellular", xlab = "Cellular UMI count",
                                y = "nUMI_MVA", ylab = "Viral UMI count") {
  p <- ggplot(tibble, mapping=aes_string(x = x, y = y)) + geom_point(aes_string(color=gene, shape=NULL)) +
    labs(title = title, x = xlab, y = ylab, color=paste0("Expression of ", gene))
  
  if(x=="nUMI_Cellular"){
    p <- p + scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10')
  }
  return(p)
}

# Plot viral or cellular UMIs on reduction
# tibble : relevant informations
# title : title of the plot
# point : nUMI_cellular or nUMI_MVA
# xlab : title of x axis
# ylab : title of y axis
# x : x axis (dimensionality reduction)
# y : y axis (dimensionality reduction)
plot_UMIs <- function(tibble, title, point = "nUMI_MVA",
                                x, xlab, y, ylab) {
  p <- ggplot(tibble, mapping=aes_string(x = x, y = y)) + geom_point(aes_string(color=point)) +
    labs(title = title, x = xlab, y = ylab)
  
  return(p)
}

# Plot the expression for every gene in a list
# Then save to PDF format
# tibble : relevant informations 
# list_genes : list of genes to plot
# x : x axis
# y : y axis
# xlab : x axis legend
# ylab : y axis legend
# pdf : output pdf path 
# text : title
gene_plot_list <- function(tibble, list_genes, x, y, xlab, ylab, pdf, text){
  list_titles <- paste("Expression of", list_genes, sep = " ")
  list_geneplot <- mapply(create_genexpr_plot, title=list_titles, gene=list_genes, SIMPLIFY = FALSE,
                          MoreArgs = list(tibble = tibble, x = x, xlab = xlab,
                                          y = y, ylab = ylab))
  gridtoPDF(list_geneplot, text, pdf)
}

# Plot the presence/absence for every gene in a list
# Then save to PDF format
# tibble : relevant informations 
# list_genes : list of genes to plot
# x : x axis
# y : y axis
# xlab : x axis legend
# ylab : y axis legend
# pdf : output pdf path 
# text : title
presence_plot_list <- function(tibble, list_genes, x, y, xlab, ylab, pdf, text){  
  list_bool <- paste(list_genes, "is", sep="_")
  list_titles <- paste0("Presence of ", list_genes)
  list_geneplot <- mapply(create_genexpr_plot, title=list_titles, gene=list_bool, SIMPLIFY = FALSE,
                          MoreArgs = list(tibble = tibble, x = x, xlab = xlab,
                                          y = y, ylab = ylab))
  gridtoPDF(list_geneplot, text, pdf)
}

# Plot viral UMIs depending on cellular UMIs
# And save it to PDF format
# tibble : relevant informations 
# pdf : output pdf path 
UMI_plot <- function(tibble, pdf){
  pdf(pdf)
  print(create_nUMI_plot(tibble, "UMIs repartition"))
  dev.off()
}


## --------- RUN ---------- ##  

if (sys.nframe() == 0){
  seurat_o <- get(load(seurat_o))
  listRDS <- readRDS("/Users/kevin/Documents/Bioinfo/MVA-infected DCs/inputs/MNN-objs-2/MNN_20_40_4.rds")
  seurat_o <- listRDS[[2]]
  my_tibble <- create_tibble(seurat_o, is_log=TRUE)
  list_genes <- c("B13R", "F1L", "N1L", "E3L", "PMAIP1")
  list_titles <- paste("Expression of", list_genes, sep = " ")
  list_geneplot <- mapply(create_genexpr_plot, title=list_titles, gene=list_genes, SIMPLIFY = FALSE,
                          MoreArgs = list(tibble = my_tibble))
  
  n <- length(list_geneplot)
  nCol <- floor(sqrt(n))
  do.call("grid.arrange", c(list_geneplot, ncol=nCol))
}
