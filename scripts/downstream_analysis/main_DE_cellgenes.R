## ---------------------------
##
## Purpose of script: Perform differential expression
## analysis between infected and uninfected cells
## in a Seurat object
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-02-14
##
## ---------------------------

## Arguments
if (sys.nframe() == 0){
  rm(list=ls())
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 4) {
    print("Not enough arguments provided")
    print("Usage : Rscript diff_expr infected uninfected jfile is_mva path")
    q()
  } else if (length(args) > 4) {
    print("Too many arguments provided")
    print("Usage : Rscript diff_expr infected uninfected jfile is_mva path")  
    q()
  }
  
  infected <- args[1]       # Seurat object with clusters
  uninfected <- args[2]     # Seurat object with clusters
  jfile <- args[3]          # Argument json file
  path <- args[4]           # Output path
  
  # Read from json files
  result <- rjson::fromJSON(file = jfile)

}


## ---------------------------

## Packages

source("utils/diff_expr.R")

## ------- FUNCTIONS ------- ##

# Main differential expression 
# infected = Seurat object of infected samples
# uninfected = Seurat object of uninfected samples
# path = output path 
# result = from json args file
main_DE_cell <- function(infected, uninfected, path, result) {
  
  set.seed("20020")
  #Log FC Threshold
  logFC = result$diffexp$logFC
  #Adjusted pvalue threshold
  pvalue = result$diffexp$pvalue
  #Minimum percentage of cells that must express the DEG
  min.pct = result$diffexp$min.pct
  #Test to use for differential expression
  test = result$diffexp$test
  #Key of the clustering solution to use
  cluster_sol = result$diffexp$cluster_sol
  #If the data should be scaled or not
  scaled = result$diffexp$scaled
  # Use the uninfected cluster or not as reference
  uninf = as.logical(result$diffexp$uninf)
  
  # Argument used
  Arguments <- c(unlist(result$diffexp), infected, uninfected)
  names(Arguments)[length(Arguments)-1] <- "infected"
  names(Arguments)[length(Arguments)] <- "uninfected"
  
  # Merge infected and uninfected
  seurat_o <- merge_inf_uninf(infected, uninfected, cluster_sol)
  
  # Create the color vector for the heatmap
  colors <- rev(RColorBrewer::brewer.pal(8, "RdBu"))
  # Create annotation of clusters
  clusters <- levels(Idents(seurat_o))[-length(levels(Idents(seurat_o)))]
  clustercol <- c(scales::hue_pal()(length(clusters)), "#000000")
  annotation_row <- cluster_annotation(seurat_o, "clustcomp")
  
  # Colors for annotation
  ann_color <- list(
    # Seurat colors
    Cluster = setNames(clustercol, levels(as.factor(Idents(seurat_o))))
  )
  
  # Matrix of heatmap
  if(scaled){
    seurat_o <- ScaleData(seurat_o, features=rownames(seurat_o))
    mat <- t(as.matrix(seurat_o[["RNA"]]@scale.data))
    
  } else {
    mat <- t(as.matrix(seurat_o[["RNA"]]@data))
  }
  mat <- mat[rownames(annotation_row),]
  
  # List to store results
  list_table <- list()
  list_top_DEG <- list()
  list_bp_UP <- list()
  list_mf_UP <- list()
  list_cc_UP <- list()
  list_kegg_UP <- list()
  list_bp_DOWN <- list()
  list_mf_DOWN <- list()
  list_cc_DOWN <- list()
  list_kegg_DOWN <- list()
  
  for(i in seq_along(clusters)){
    print(clusters[i])
    ## Differential expression
    # Remove scientific annotation for easy comparisons
    options(scipen=999)
    dir.create(path)
    if(uninf){
      markers <- sc_diffexp(seurat_o, clusters[i], "Uninfected", logFC, pvalue, test, min.pct)
    } else {
      markers <- sc_diffexp(seurat_o, clusters[i], clusters[-i], logFC, pvalue, test, min.pct)
    }
    # Put scientific annotation back
    options(scipen=0)
    # Remove viral genes
    markers <- filter_viral(markers, seurat_o)
    # Sort by absolute logFC and write genes
    marker_tibble <- as_tibble(markers) %>%
      add_column(Gene = rownames(markers)) %>%
      arrange(desc(abs(avg_logFC)))
    write_tsv(marker_tibble, file.path(path, paste0(clusters[i], "_diff_genes.tsv")))
    print("Writing tibble")
    # Add them to list of differential expression
    marker_tibble <- marker_tibble %>% add_column(cluster = rep(clusters[i], dim(marker_tibble)[1]))
    list_table[[i]] <- marker_tibble
    # Add the genes to the list of top DEG
    list_top_DEG[[i]] <- marker_tibble[1:50, "Gene"]$Gene
    
    ## Heatmap
    
    # Heatmap
    if(scaled){myBreaks <- seq(-2,2,by=0.5)}else{myBreaks <- NULL}
    heatM <- heatmap_DEG(mat[, rownames(markers)], colors, annotation_row, ann_color, myBreaks = myBreaks)
    pdf(file.path(path, paste0(clusters[[i]], "_heatmap.pdf")), height=7, width=18)
    print(heatM)
    dev.off()
    
    ## Kegg
    
    # More expressed
    dir.create(file.path(path, "UP"))
    dir.create(file.path(path, "DOWN"))
    up <- markers[markers$avg_logFC > 0,]
    up_enr <- main_diff(up, 
                        file.path(path, "UP", paste0(clusters[i], "_kegg.tsv")),
                        file.path(path, "UP", paste0(clusters[i], "_kegg.pdf")),
                        do_KEGG)
    
    # More repressed
    down <- markers[markers$avg_logFC < 0,]
    down_enr <- main_diff(down, 
                          file.path(path, "DOWN", paste0(clusters[i], "_kegg.tsv")),
                          file.path(path, "DOWN", paste0(clusters[i], "_kegg.pdf")),
                          do_KEGG)
    list_kegg_UP[[i]] <- up_enr
    list_kegg_DOWN[[i]] <- down_enr
    
    ## GO BP
    
    # More expressed
    up_enr <- main_diff(up, 
                        file.path(path, "UP", paste0(clusters[i], "_GO-BP.tsv")),
                        file.path(path, "UP", paste0(clusters[i], "_GO-BP.pdf")),
                        do_GO, "BP")
    # More repressed
    down_enr <- main_diff(down, 
                          file.path(path, "DOWN", paste0(clusters[i], "_GO-BP.tsv")),
                          file.path(path, "DOWN", paste0(clusters[i], "_GO-BP.pdf")),
                          do_GO, "BP")
    list_bp_UP[[i]] <- up_enr
    list_bp_DOWN[[i]] <- down_enr
    
    ## GO MF
    
    # More expressed
    up_enr <- main_diff(up, 
                        file.path(path, "UP", paste0(clusters[i], "_GO-MF.tsv")),
                        file.path(path, "UP", paste0(clusters[i], "_GO-MF.pdf")),
                        do_GO, "MF")
    # More repressed
    down_enr <- main_diff(down, 
                          file.path(path, "DOWN", paste0(clusters[i], "_GO-MF.tsv")),
                          file.path(path, "DOWN", paste0(clusters[i], "_GO-MF.pdf")),
                          do_GO, "MF")
    list_mf_UP[[i]] <- up_enr
    list_mf_DOWN[[i]] <- down_enr
    
    ## GO CC
    
    # More expressed
    up_enr <- main_diff(up, 
                        file.path(path, "UP", paste0(clusters[i], "_GO-CC.tsv")),
                        file.path(path, "UP", paste0(clusters[i], "_GO-CC.pdf")),
                        do_GO, "CC")
    # More repressed
    down_enr <- main_diff(down, 
                          file.path(path, "DOWN", paste0(clusters[i], "_GO-CC.tsv")),
                          file.path(path, "DOWN", paste0(clusters[i], "_GO-CC.pdf")),
                          do_GO, "CC")
    list_cc_UP[[i]] <- up_enr
    list_cc_DOWN[[i]] <- down_enr
  }
  
  # Union of the pathways
  top_pathways(list_bp_UP, path, "Pathways_GO-BP_UP.tsv")
  top_pathways(list_bp_DOWN, path, "Pathways_GO-BP_DOWN.tsv")
  top_pathways(list_cc_UP, path, "Pathways_GO-CC_UP.tsv")
  top_pathways(list_cc_DOWN, path, "Pathways_GO-CC_DOWN.tsv")
  top_pathways(list_mf_UP, path, "Pathways_GO-MF_UP.tsv")
  top_pathways(list_mf_DOWN, path, "Pathways_GO-MF_DOWN.tsv")
  top_pathways(list_kegg_UP, path, "Pathways_Kegg_UP.tsv")
  top_pathways(list_kegg_DOWN, path, "Pathways_Kegg_DOWN.tsv")
  
  ## Create a union of the top 50 DEG in each case
  DEG <- unionDEG(list_top_DEG)
  heatM <- heatmap_DEG(mat[, DEG], colors, annotation_row, ann_color, myBreaks = myBreaks)
  pdf(file.path(path, "DEG_heatmap.pdf"), height=7, width=18)
  print(heatM)
  dev.off()
  
  ## Create a recap table with all DEG and
  ## in which cluster they are differentially expressed
  # Data for LogFC Recap heatmap
  my_table <- data_RecapHM(list_table, path, 2)
  myBreaks <- seq(-2.5,4.5,by=1)
  # Tweak colors to have white color between -0.5 and 0.5 logFC
  colors_logFC <- rev(RColorBrewer::brewer.pal(9, "RdBu"))[-1:-2]
  mat <- as.matrix(my_table[-1])
  rownames(mat) <- my_table$Gene
  recapHM <- recap_heatmap(t(mat), colors_logFC, myBreaks = myBreaks)
  pdf(file.path(path, "LogFC_Heatmap.pdf"), height=7, width=18)
  print(recapHM)
  dev.off()
  
  # Data for average expression matrix
  matavg_t <- data_AvgHM(seurat_o, DEG, scaled)
  heatM <- avgHM(matavg_t[,DEG], colors)
  pdf(file.path(path, "avg_expression_DEG.pdf"), height=7, width=18)
  print(heatM)
  dev.off()
  
  # Write arguments used
  write.table(as.data.frame(Arguments), file.path(path, "list_arguments_diffexp.tsv"), sep="\t") 
}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  main_DE_cell(infected, uninfected, path, result)
}
