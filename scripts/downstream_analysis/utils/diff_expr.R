## ---------------------------
##
## Purpose of script: Perform differential expression
## analysis between groups of cells in a Seurat object
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2019-12-10
##
## ---------------------------

## Packages

library(Seurat)
library(clusterProfiler)
library("org.Hs.eg.db")
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)


## ------ FUNCTIONS ------- ##

# Merge infected and uninfected together
# infected = RDS object of infected Seurat object with cluster solutions
# uninfected = uninfected R object
# cluster_sol = name of the cluster solution
merge_inf_uninf <- function(infected, uninfected, cluster_sol) {
  # Load infected and uninfected files, and merge together
  infected <- readRDS(infected)
  uninfected <- get(load(uninfected))
  seurat_o <- merge(infected, uninfected)
  seurat_o[["RNA"]]@meta.features <- infected[["RNA"]]@meta.features
  # Create a metadata with the clusters + a "Uninfected" cluster
  times <- length(seurat_o@meta.data[[cluster_sol]]) - length(infected@meta.data[[cluster_sol]])
  clustcomp <- c(as.vector(infected@meta.data[[cluster_sol]]), rep("Uninfected", times))
  seurat_o <- AddMetaData(seurat_o, clustcomp, "clustcomp")
  Idents(seurat_o) <- clustcomp
  return(seurat_o)
}

# Differential expression with Seurat's FindMarkers function
# seurat_o = Seurat object with clusters 
# comp = vector of clusters to find DE for
# ref = vector of clusters used for reference
# logfc = logFC threshold
# pvalue = pvalue threshold 
# test = statistical test (must be implemented in Seurat's FindMarkers)
# min.pct = threshold of the percent of cells expressing a DE gene
# outfile = file where to write results with all genes
sc_diffexp <- function(seurat_o, comp, ref, logfc, pvalue, test, min.pct) {
  markers <- FindMarkers(seurat_o, ident.1 = comp, ident.2 = ref, logfc.threshold	= logfc, min.pct = min.pct,
                         test.use = test, latent.vars	= c("cell_id_sil"), assay = "RNA")
  markers$p_val_adj <- sapply(markers$p_val_adj, round, digits=10)
  markers <- markers[markers$p_val_adj < pvalue, ]
  return(markers)
}

# Filter viral genes from rownames of dataframe
# df = dataframe of differential genes in rownames
# seurat_o = Seurat object
filter_viral <- function(df, seurat_o){
  metafeatures <- seurat_o[["RNA"]]@meta.features
  df <- df[(rownames(df) %in% rownames(metafeatures[metafeatures$is_viral == 0,])),]
  return(df)
}

# Filter cellular genes from rownames of dataframe
# df = dataframe of differential genes in rownames
# seurat_o = Seurat object
filter_cellular <- function(df, seurat_o){
  metafeatures <- seurat_o[["RNA"]]@meta.features
  df <- df[(rownames(df) %in% rownames(metafeatures[metafeatures$is_viral == 1,])),]
  return(df)
}

# Enrichment analysis with GO annotations
# ids = list of genes
# ont = type of ontologie (CC, BP, MF)
# keyType = annotation of genes in list
do_GO <- function(ids, ont, keyType="ENTREZID"){
  enrich <- enrichGO(ids, OrgDb = org.Hs.eg.db, ont = ont, keyType = keyType)
  enrich <- simplify(enrich, cutoff = 0.7, by ="p.adjust", select_fun = min)
  return(enrich)
}

# Enrichment analysis with KEGG annotations
# ids = list of genes
do_KEGG<- function(ids){
  return(enrichKEGG(ids))
}

# Enrichment analysis function
# df = dataframe with differential genes in rownames
# pdf_file = name of the pdf_file with output plots 
# FUN = which enrichment function to use (do_KEGG or do_GO)
# + additional arguments passed to FUN
sc_enrichment <- function(df, pdf_file, FUN, ...){
  # Transform from Symbol to Entrez
  entrez_ids <- bitr(rownames(df), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)
  # Do the enrichment analysis
  enrich <- FUN(entrez_ids$ENTREZID, ...)
  
  # If the enrichment analysis is conclusive
  if(!is.null(enrich)){
    for(i in seq_along(enrich@result$geneID)){
      # Extract individual genes in the pathways
      init <- ""
      ids <- unlist(strsplit(enrich@result$geneID[i], split = "/"))
      for (j in seq_along(ids)){
        j <- entrez_ids[(match(ids[j], entrez_ids$ENTREZID)),"SYMBOL"]
        init <- paste(init, j, sep = "/")
      }
      init <- sub('/', '', init)
      enrich@result$geneID[i] <- init
    }
    
    
    geneList <- df[rownames(df) %in% entrez_ids[!is.na(entrez_ids$ENTREZID),]$SYMBOL,]
    foldChange <- geneList$avg_logFC
    names(foldChange) <- rownames(geneList)
    if(dim(data.frame(enrich))[1] != 0){
      pdf(pdf_file)  
      print(dotplot(enrich, showCategory = 10, font.size = 8))
      print(cnetplot(enrich, showCategory = 10, font.size = 8, categorySize = "pvalue", foldChange = foldChange))
      dev.off()
    }
    # Return the enrichment tsv
    return(enrich)
  }
  # Or return NULL if analysis inconclusive
  return(NULL)
}

# Subset dataframe of diff exp with genes in enrichment
# enrich = results from enrichment analysis
# df = dataframe with differential genes in rownames
filter_enrich <- function(enrich, df){
  genes <- unique(unlist(lapply(enrich@result$geneID, strsplit, split = "/")))
  df <- df[rownames(df) %in% genes,]
  return(df)
}

# Main enrichment
# df = dataframe with DE in rownames
# path = output path 
# pdf = pdf name 
# FUN = do_KEGG or do_GO
#+additional args passed to FUN
main_diff <- function(df, path, pdf, FUN, ...){
  enrich <- sc_enrichment(df, pdf, FUN, ...)
  if(!is.null(enrich)){
    write_tsv(as.tibble(enrich), path)
  }
  return(as.tibble(enrich))
}

# Create an annotation for the genes for the heatmap
# Specifically, add the class of the gene (cellular, early, interm, late and unknown)
# seurat_o = Seurat object
gene_annotations <- function(seurat_o) {
  gene_annotations <- as_tibble(seurat_o[["RNA"]]@meta.features["gene_class"]) %>% 
    add_column(Gene = rownames(seurat_o[["RNA"]]@meta.features), .before=1)
  colnames(gene_annotations) <- c("Gene", "Class")
  gene_annotations = as.data.frame(gene_annotations)
  # Transform to dataframe for pheatmap
  rownames(gene_annotations) <- gene_annotations$Gene
  gene_annotations <- gene_annotations[-1]
  # Convert to vector for pheatmap
  gene_annotations$Class <- as.factor(gene_annotations$Class)
  return(gene_annotations)
}

# Create an annotation for the clusters for the heatmap
# seurat_o = Seurat_object
# cluster_sol = key of the solution cluster
cluster_annotation <- function(seurat, cluster_sol) {
  cluster_ann <- as.data.frame(seurat@meta.data[[cluster_sol]])
  colnames(cluster_ann) <- "Cluster"
  # Order cells by cluster number
  rownames(cluster_ann) <- colnames(seurat[["RNA"]])
  cluster_ann <- cluster_ann[order(cluster_ann$Cluster), , drop = FALSE]
  # Convert to factor for pheatmap
  cluster_ann$Cluster <- as.factor(cluster_ann$Cluster)
  return(cluster_ann)
}


# Create a heatmap with genes of interests (e.g., DEG)
# mat = matrix with genes in columns, cells in rows
# colors = vector of color
# row_annot = dataframe of cell annotation 
# ann_color = colors for annotations 
# cols_annotation = dataframe of gene annotation 
# cluster_rows = wether to cluster cells together
# cluster_cols = wether to cluster genes together 
# myBeaks = numerical vector of breaks
heatmap_DEG <- function(mat, colors, row_annot, ann_color, cols_annotation = NULL, show_rownames = FALSE,
                        show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = TRUE, myBreaks = NULL,
                        gaps_col = NULL, gaps_row = NULL,
                        main = "Gene expression values") {
  
  
  HM <- pheatmap::pheatmap(mat, breaks = myBreaks, gaps_col = gaps_col, gaps_row = gaps_row,
                           cluster_rows = cluster_rows, cluster_cols = cluster_cols, show_rownames = show_rownames, 
                           show_colnames = show_colnames,
                           annotation_col = cols_annotation, annotation_row = row_annot, annotation_colors = ann_color,
                           color=colors, fontsize = 7, scale="none", clustering_method="complete", main = main)
  return(HM)
}


# Create a dataframe recaping which genes are differential
# in which cluster, and with which logFC
# Used for the recap_heatmap function
# list_tables = list of DE tables
# path = output path 
# threshold = threshold of logFC, to restrict to top hits
data_RecapHM <- function(list_table, path, threshold = NULL){
  # Create a tibble with all marker_tibbles
  my_table <- dplyr::bind_rows(list_table)
  for(i in unique(my_table$cluster)){
    # Column for cluster number i
    column1 <- rlang::parse_expr(paste0("cluster", (i)))
    # Column for logFC value in cluster number i
    column2 <- rlang::parse_expr(paste0("logFC_c", (i)))
    my_table <- my_table %>% 
      # Boolean column : 1 if differential in this cluster
      mutate(!!column1 := ifelse(cluster == (i), 1, 0)) %>%
      # Value is LogFC if differential, else 0
      mutate(!!column2 := ifelse(cluster == (i), avg_logFC, 0))
  }
  
  # Collapse table
  my_table <- my_table %>% 
    dplyr::select(-p_val, -avg_logFC, -pct.1, -pct.2, -p_val_adj) %>% 
    group_by(Gene) %>% 
    # There is i value of logFC and cluster for each gene
    # all in i different columns 
    # therefore we can sum while collapsing to keep one line per gene
    # Without modyfing the logFC for a cluster
    summarize_if(is.numeric, sum)
  
  # Write the table
  write_tsv(my_table, file.path(path, "resume_diffgenes.tsv"))
  
  # Apply a threshold of LogFC (optional)
  if(!is.null(threshold)){
    my_table <- dplyr::filter_if(my_table, is.numeric, any_vars(abs(.) > threshold))
  }
  
  # Rename some variable for Heatmap
  my_table <- dplyr::select(my_table, -starts_with("cluster"))
  names(my_table) <- names(my_table) %>%
    gsub("logFC_c", "Cluster ", .) 
  
  return(my_table)
}


# Heatmap with logFC values of genes across clusters
# mat = matrix with genes in columns, cells in rows
# colors = vector of color
# ann_color = colors for annotations 
# annotation_col = dataframe of gene annotation 
# myBeaks = numerical vector of breaks
recap_heatmap <- function(mat, colors, 
                          annotation_col = NULL, ann_color = NULL, myBreaks = NULL) {
  print(mat)
  if(is.null(myBreaks)){
    myBreaks <- c(seq(min(mat), 0, length.out=ceiling(11/2) + 1), 
                  seq(max(mat)/11, max(mat), 
                      length.out=floor(11/2)
                  )
    )
  }
  print(myBreaks)

  logFCheatmap <- pheatmap::pheatmap(mat, color=colors, 
                                     cluster_rows = FALSE, cluster_cols = TRUE,
                                     display_numbers = TRUE, breaks=myBreaks, 
                                     annotation_col = annotation_col, ann_color = ann_color,
                                     cellwidth = 22, fontsize = 12, cellheight = 22, fontsize_number = 8)
  
  return(logFCheatmap)
}

# Create a matrix of average expression across clusters 
# seurat_o = Seurat object
# genes = genes of interest 
# scaled = wether to use scaled data or log-normalized
data_AvgHM <- function(seurat_o, genes, scaled){
  if(scaled){slot="scale.data"}else{slot="data"}
  avg <- AverageExpression(seurat_o, "RNA", genes, slot=slot)
  matavg <- as.matrix(avg$RNA)
  matavg_t <- t(matavg)
  rownames(matavg_t) <- colnames(matavg)
  #Rownames sorting
  clusters <- sort(rownames(matavg_t))
  if(clusters[length(clusters)] == "Uninfected"){
    # Put Uninfected first
    clusters <- clusters[c(length(clusters),1:(length(clusters)-1))]
  }
  matavg_t <- matavg_t[clusters,]
  return(matavg_t)
}

# Create a heatmap with average expression of genes across clusters
# mat = matrix with genes in columns, cells in rows
# colors = vector of color
# ann_color = colors for annotations 
# cols_annotation = dataframe of gene annotation 
# scaled = wether to used scaled value or log normalized 
# myBeaks = numerical vector of breaks
# ... = additional arguments to pheatmap
avgHM <- function(mat, colors, row_annot = NULL,
                  cols_annotation = NULL, ann_color = NULL,
                  scaled = TRUE, myBreaks = NULL,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  gaps_col = NULL, gaps_row = NULL,
                  title = "Average expression values",
                  ...) {
  
  quantile_breaks <- function(xs, n = 10){
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    return(breaks[!duplicated(breaks)])
  }
  
  if(is.null(myBreaks)){
    if(scaled){
      myBreaks <- seq(-2,2,by=0.5)
    }
    else{
      myBreaks <- quantile_breaks(matavg, n = 11)
    }
  }
  
  avgHM <- pheatmap::pheatmap(mat, color=colors, breaks = myBreaks, gaps_col = gaps_col, gaps_row = gaps_row,
                              cluster_rows = cluster_rows, cluster_cols = cluster_cols, fontsize = 8, 
                              annotation_col = cols_annotation, annotation_row = row_annot, annotation_colors = ann_color,
                              main = title, ...)
  return(avgHM)
}

# Return the union of top pathways
# list_enrich = list of results of pathways enrichment analysis 
# filepath = output path 
# filename = name of the file
top_pathways <- function(list_enrich, filepath, filename){
  for(i in seq_along(list_enrich)){
    if(i == 1){
      my_tibble <- list_enrich[[1]] %>% top_n(10)
    }
    else{
      my_tibble2 <- list_enrich[[i]] %>% top_n(10)
      my_tibble <- full_join(my_tibble, my_tibble2, 
                             by=c("ID", "Description"),
                             suffix=c(as.character(i-1), as.character(i)))
    }
  }
  my_tibble <- my_tibble %>% dplyr::select(ID, Description, contains("p.adjust"), contains("Count"))
  write_tsv(my_tibble, file.path(filepath, filename))
}

# Return the union of a list of differential genes
# list_top_DEG = list of differential genes
unionDEG <- function(list_top_DEG) {
  for(i in seq_along(list_top_DEG)[-length(list_top_DEG)]){
    un <- union(list_top_DEG[[i]], list_top_DEG[[i+1]])
    if(i < (length(list_top_DEG)-1)){
      list_top_DEG[[i+1]] <- un
    }
    else{
      return(un)
    }
  }
}

