## ---------------------------
##
## Purpose of script: Evaluate the relevance
## of a single cell clustering solution
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-01-20
##
## ---------------------------

## Arguments

if (sys.nframe() == 0){
  args <- commandArgs(trailingOnly = TRUE)
  
  seurat_o <- args[1]        # Path to R RDS Seurat object of interest
  jfile <- args[2]           # Json argument file
  path <- args[3]            # Path where to write the output
  
  # Read from json files
  result <- rjson::fromJSON(file = jfile)
}


## ---------------------------

## Packages

library(Seurat)
source("utils/various_plots.R")
source("utils/sc_utils.R")

## ------ FUNCTIONS ------- ##

# Create a tibble with relevant information for the plots
# seurat_o = Seurat object
clust_tibble <- function(seurat_o){
  #Barplots
  barplot_t <- create_tibble(seurat_o, c("mcherry"), TRUE, FALSE, 0) %>% 
    add_column(sample = seurat_o$orig.ident) %>% 
    add_column(cluster = seurat_o$clusters) %>%
    dplyr::select(Cell, sample, cluster, group) %>% 
    mutate(group = ifelse(group == "Dying", "HVLC",
                          ifelse(group == "Resist", "HVHC", "LVHC"))) %>%
    group_by(sample, cluster, group) %>%
    summarise(ncells = n())
    return(barplot_t)
}

# Barplots of the number of cells for each cluster
# with the possibility of representing the proportion 
# of a variable of interest (ex : donor distribution)
# tibble = result from clust_tibble()
# fill = what to represent on the plot
barplot_clusters <- function(tibble, fill){
  g <- ggplot(tibble, aes_string(x="cluster", y="ncells", fill=fill)) + 
    geom_bar(stat="identity", width= 0.3)
  return(g)
}

# VlnPlot for each cluster
# seurat_o = Seurat object
# features = list of features to plut
violin <- function(seurat_o, features){
  seurat_o$ratioMVA <- seurat_o$nUMI_MVA/seurat_o$nCount_RNA
  for(f in features){
    if(f %in% c("nUMI_MVA", "nUMI_cellular")){log = TRUE}else{log = FALSE}
    v <- VlnPlot(seurat_o, features = f, group.by = "clusters", pt.size = 0.5, log = log)
    print(v)
  }
}

# Create the violin plots and barplots
# seurat_o = Seurat object
# result = from json args file
# path = output path
main_cluster_eval <- function(seurat_o, result, path){
  
  # Key of the cluster configuration of interest
  seurat_cluster = result$clustplots$seurat_cluster
  # List of genes of interest
  features = result$clustplots$features
  
  dir.create(file.path(path))
  seurat_o <- readRDS(seurat_o)
  seurat_o$clusters <- seurat_o@meta.data[[seurat_cluster]]
  my_tibble <- clust_tibble(seurat_o)
  barplot_cluster <- barplot_clusters(my_tibble, "group")
  barplot_sample <- barplot_clusters(my_tibble, "sample")
  pdf(file.path(path, "cluster_barplots.pdf"))
  print(barplot_cluster)
  print(barplot_sample)
  dev.off()
  pdf(file.path(path, "violinplots.pdf"))
  l <- violin(seurat_o, features)
  dev.off()
}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  main_cluster_eval(seurat_o, seurat_cluster, features, path)
}