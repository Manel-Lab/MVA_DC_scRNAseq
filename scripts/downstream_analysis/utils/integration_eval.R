#!/usr/bin/Rscript
## ---------------------------
##
## Functions for evaluating the efficency of integration
## methods and batch effect removal methods
## Tools used : Silhouette, kBET
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Based on Tarek Gharsalli's scripts
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-11-14

## ---------------------------

## Packages

library(Seurat)
library(kBET)
library(cluster)
library(ggplot2)
library(factoextra)
library(tibble)
library(readr)

## ------ FUNCTIONS ------- ##

# Function that calculate batch effect estimation 
# Using the kBET methods (see paper)
# data = a dataframe
# batch = a vector of batches
kBET_estimate <- function(data, batch) {
  batch.estimate <- kBET(data, batch, plot=FALSE, verbose=TRUE)
  # Prepare data for plotting
  kBETdata <- data.frame(class=rep(c('observed', 'expected'), 
                                   each=length(batch.estimate$stats$kBET.observed)), 
                         data =  c(batch.estimate$stats$kBET.observed,
                                   batch.estimate$stats$kBET.expected))
  return(kBETdata)
}

# Function that computes the Silhouette score for
# our batches of interest
# data = a dataframe
# batch = a vector of batches
sil_estimate <- function(data, batch) {
  # Compute distance on data matrix format
  dist_data <- as.matrix(dist(data))
  # Calculate the silhouette score
  clust <- as.numeric(batch)
  silscore <- silhouette(x = clust, dist = dist_data)
  return(silscore)
}

# Plot kBET results
# kBET_data = results of kBET_estimate()
# name_K = plot name
kBET_plotting <- function(kBET_data, name_K) {
  kBETplot <- ggplot(kBET_data, aes(class, data)) + geom_boxplot() + 
    labs(x='Test', y='Rejection rate', title=name_K) +
    theme_bw() +  
    scale_y_continuous(limits=c(0,1))
  return(kBETplot)
}

# Plot Silhouette scores for our batches
# sil_score = results of sil_estimate()
# name_S = plot name
sil_plotting <- function(sil_score, name_S) {
  silplot <- fviz_silhouette(sil_score, print.summary = FALSE, main = name_S +
                               scale_fill_brewer(palette = "Dark2") +
                               scale_color_brewer(palette = "Dark2") +
                               theme_minimal() +
                               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  return(silplot)
}


# Create a diagnostic tibble for MNN methods
# With average silhouette and kBET score aswell as
# Percentage of variance lost at each merge
# redlist = list of Seurat objects
# outpath = output path
# (defunct)
diagnostic_tibble <- function(redlist, outpath){
  print("Creating diagnostic tibble")
  
  lost_var1 <- c()
  lost_var2 <- c()
  message <- c()
  silhouette <- c()
  kBET <- c()
  i <- 1
  
  for(i in seq_along(redlist)){
    obj <- redlist[[i]]
    mnn <- obj[[8]]
    k[i] <- obj[[2]]$k
    d[i] <- obj[[2]]$d
    ndist[i] <- obj[[2]]$ndist
    lost_var1[i] <- round(mnn@tools$RunFastMNN@metadata$merge.info$lost.var[1,1], digits = 2)
    lost_var2[i] <- round(mnn@tools$RunFastMNN@metadata$merge.info$lost.var[1,2], digits = 2)
    silhouette[i] <- round(mean(obj[[4]][,3]), digits = 2)
    kBET[i] <- round(mean(obj[[3]][obj[[3]][,"class"] == "observed","data"]), digits = 2)
    i <- i + 1
    gc()
  }
  
  print("Writing tibble")
  my_tibble <- as_tibble(list(k=k, d=d, ndist=ndist, 
                              variance_lost_donor1=lost_var1, 
                              variance_lost_donor2=lost_var2,
                              silhouette=silhouette, kBET=kBET))
  write_tsv(my_tibble, file.path(outpath, "diagnostic_table.tsv"))
}
