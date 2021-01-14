#!/usr/bin/Rscript
## ---------------------------
##
## Cluster the chosen reduction solution
## Compute Silhouette score of the clusters
## And plot various informations
## All to help in the decision of which reduction method to use
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-11-21

## ---------------------------

## Arguments

if (sys.nframe() == 0){
  
  .libPaths(c("~/R_packages", .libPaths()))
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    print("Not enough arguments provided")
    print("Usage : Rscript ")
    q()
  } else if (length(args) > 3) {
    print("Too many arguments provided")
    print("Usage : Rscript ")  
    q()
  }
  
  robj <- args[1]         # Seurat object
  outpath <- args[2]         # Path to write the output
  jfile <- args[3]        # Argument json file
  
  ## Global variables
  
  # Read arguments from json files
  result <- rjson::fromJSON(file = jfile)
}


  
## ---------------------------

## Packages

source("utils/clustering.R")
source("utils/various_plots.R")
source("utils/sc_utils.R")
source("utils/integration_eval.R")

## ------ FUNCTIONS ------- ##

# This function cluster our data and compute Silhouette score
# robj = path to Seurat object
# outpath = output path
# result = from json args file
cluster_and_sil <- function(robj, outpath, result){
  
  # Key of the reduction of interest
  red = result$eval$red
  # Vector of different k neighbors to test
  res = result$eval$res
  # Number of dimensions
  d = result$eval$d
  
  # Clustering with default k but varying level of resolutions
  print("Clustering")
  robj <- Neighborparams(robj, 30, red, d)
  listcluster <- lapply(res, Clusterparams, robj=robj)
  
  # Only keep indexes of the list where the Seurat
  # object has more clusters than previously
  a <- 0
  i <- 1
  list_index <- c()
  print("Keeping only relevant indexes")
  for(object in listcluster){
    n_clusters <- levels(object@meta.data$seurat_clusters)
    if(length(n_clusters) > a){
      list_index <- append(list_index, i)
      a <- length(levels(object@meta.data$seurat_clusters))
    }
    i <- i + 1
  }
  
  # Reverse order of list
  list_index <- rev(list_index)
  # Keep objects of interest
  list_seurat <- c(listcluster[list_index])
  dir.create(file.path(outpath, "Clusters"))
  
  # Apply UMAP on the reduction and plot the clusters
  print("UMAP reduction")
  listumap <- lapply(list_seurat, UMAPparams, red=red, d=d)
  listplotumap <- lapply(listumap, UMAPplot, group = "seurat_clusters", title = paste0("UMAP plot, ", red))

  # Plot the reduction but grouping by sample of origin
  plotumap <- UMAPplot(robj, group = "orig.ident", title = paste0("UMAP plot, ", red))
  
  # Save UMAP plots to PDF format
  gridtoPDF(listplotumap, 
            paste0("UMAP plot, reduction = ", red), 
            file.path(outpath, "Clusters", paste0("UMAP_clusters_", red, ".pdf")))
  
  pdf(file.path(outpath, "Clusters", paste0("UMAP_clusters_orig_", red, ".pdf")))
  print(plotumap)
  dev.off()

  # Compute Silhouette scores for clusters
  listsilscore <- lapply(list_seurat, computeSilhouette, red)
  listsilscore <- listsilscore[!is.na(listsilscore)]
  # Plot Silhouette scores
  listsilplot <- lapply(listsilscore, plotSilhouette)

  # Silhouette plots to PDF
  gridtoPDF(listsilplot, 
            paste0("Silhouette plots, reduction = ", red), 
            file.path(outpath, "Clusters", paste0("Silhouette_clusters_", red, ".pdf")))
} 

# Main function to evaluate the validity of a reduction
# robj = Seurat object path
# outpath = output path
# result = from jsona args file
main_red_eval <- function(robj, outpath, result){
  
  # Key of the reduction of interest
  red = result$eval$red
  # Number of dimensions
  d = result$eval$d
  # Vector of different k neighbors to test
  list_genes = result$eval$list_genes
  # Use normalized data for plot ?
  is_norm = result$eval$is_norm
  # Log10-transform data for plot ?
  is_log = result$eval$is_log
  # Detection threshold
  threshold = result$eval$threshold
  # Wether or not do the clustering part
  clustering = result$eval$clustering
  # Wether or not do the gene plots part
  geneplot = result$eval$geneplot
  
  # Create the tibble
  dir.create(outpath)
  robj <- readRDS(robj)
  robj <- UMAPparams(robj, red=red, d=d)
  
  if(clustering){
    # Clustering and Silhouette
    cluster_and_sil(robj, outpath, result)
  }
  
  my_tibble <- create_tibble(robj, list_genes, is_norm, is_log, threshold)
  
  if(geneplot){  
    # Gene plot on rawdata
    gene_plot_list(my_tibble, list_genes, 
                   "nUMI_Cellular", "nUMI_MVA", "Cellular UMI count", "MVA UMI count", 
                   file.path(outpath, "Clusters", paste0("Cell_Death_genes_rawdata_", red, ".pdf")),
                   "Cell death inhibiting MVA genes expression on normalized data"
    )
    # Presence plot on rawdata
    presence_plot_list(my_tibble, list_genes, 
                       "nUMI_Cellular", "nUMI_MVA", "Cellular UMI count", "MVA UMI count", 
                       file.path(outpath, "Clusters", paste0("Presence_Cell_Death_genes_rawdata", red, ".pdf")),
                       "Cell death inhibiting MVA genes prense on normalized data"
    )
    # Presence plot on rawdata
    UMI_plot(my_tibble, file.path(outpath, "Clusters", paste0("UMI_repartition_", red, ".pdf")))
    # UMI counts on reductions
    pdf(file.path(outpath, "Clusters", paste0("UMI counts in reduction", red, ".pdf")))
    p1 <- plot_UMIs(my_tibble, "UMI Cellular", "nUMI_Cellular", 
                    "UMAP_1", "Component 1", "UMAP_2", "Component 2")
    p2 <- plot_UMIs(my_tibble, "UMI MVA", "nUMI_MVA", 
                    "UMAP_1", "Component 1", "UMAP_2", "Component 2")
    grid.arrange(p1, p2,
                 ncol=2, nrow=1, top="Normalized UMI count in both replicates")
    dev.off()
  }
}

## --------- RUN ---------- ##


if (sys.nframe() == 0){
  main_red_eval(robj, outpath, result)
}


