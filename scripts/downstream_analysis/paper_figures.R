#!/usr/bin/Rscript
## ---------------------------
##
## Figures for the paper
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
## Date Created: 2019-12-02
##
## ---------------------------

if (sys.nframe() == 0){
  .libPaths(c("~/R_packages",.libPaths()))
  
  ## Arguments
  
  args <- commandArgs(trailingOnly = TRUE)
  
  robj_file <- args[1]       # Path to R RDS object of interest
  uninfected_file <- args[2] # R object for uninfected samples
  jfile <- args[3]           # Json argument file
  outpath <- args[4]         # Path where to write the output
  indir <- args[5]           # Directory with the DE and pathway tables
  
  # Read from json files
  result <- rjson::fromJSON(file = jfile)
}

## ---------------------------

## Packages
source("utils/various_plots.R")
source("utils/clustering.R")
source("utils/sc_utils.R")
source("utils/diff_expr.R")
library(tibble)
library(tidyr)
library(readr)
library(gridExtra)
library(ggpmisc)

## ------ FUNCTIONS ------- ##

# Personal png function for this session
png2 <- function(path){
  png(path, width = 1000, height = 1000)
}

# FeaturePlot of genes of interest
# robj = Seurat object
# red = key of reduction of interest
# d = number of dimensions
# g = gene of interest
# outpath = output path for PDF
expr_plots <- function(robj, red, d, g, outpath) {

  umap <- UMAPplot(robj, reduction = red,
                   group = "orig.ident", title = paste0("UMAP plot of MNN reduction on infected samples")
                   ) 
  umap <- umap + ylim(-7.5,5)
  pdf(file.path(outpath, paste(red, "umap_plot.pdf", sep="_")))
  print(umap)
  dev.off()
  png2(file.path(outpath, paste(red, "umap_plot.png", sep="_")))
  print(umap)
  dev.off()
  
  robj$mcherry <- robj$mcherry_n
  if(g %in% c("mcherry", "mcherry_n", "nUMI_MVA_norm")){
    cols = c("#CFCFCF", "#FF0000")
  }
  else{
    cols = c("#CFCFCF", "#0000FF")
  }
  
  plot <- FeaturePlot(robj, g, pt.size = 1,
                      cols = cols,
                      reduction=red) + ylim(-7.5,5)
  path1 <- file.path(outpath, paste(g, "expr", red, "umap.pdf", sep = "_"))
  path2 <- file.path(outpath, paste(g, "expr", red, "umap.png", sep = "_"))
  
  # Expression plot
  pdf(path1)
  print(plot + 
          labs(title = paste0("Expression of ", g)))
  dev.off()
  
  png2(path2)
  print(plot + 
          labs(title = paste0("Expression of ", g)))
  dev.off()
}

# Create table with proportion of a gene class among expressed gene
# robj = Seurat R object of infected cells
viral_counts_tibble_gene <- function(robj) {
  # Get Assay Counts and class of MVA genes
  mycounts <- as.matrix(robj[["RNA"]]@counts)

  # Keep counts of viral genes only
  is_mva <- robj@assays$RNA@meta.features$is_viral
  features <- as.vector(rownames(robj@assays$RNA@counts))
  keep <- features[which(is_mva == 1)]
  my_counts <- as_tibble(t(mycounts[keep,])) %>% 
    add_column("Cell" = colnames(robj[["RNA"]]@counts), .before = 1)
  # Determine the class of each gene
  early_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "early",])
  int_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "interm",])
  late_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "later",])
  unk_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "unknown",])
  
  # Create sub tibble of early, int, late and unknown genes
  cells <- my_counts$Cell
  my_counts_e <- my_counts[c(early_g)] %>% mutate(n_e = rowSums(. != 0)) %>% add_column(Cell = cells)
  my_counts_i <- my_counts[c(int_g)] %>% mutate(n_i = rowSums(. != 0)) %>% add_column(Cell = cells)
  my_counts_l <- my_counts[c(late_g)] %>% mutate(n_l = rowSums(. != 0)) %>% add_column(Cell = cells)
  my_counts_u <- my_counts[c(unk_g)] %>% mutate(n_u = rowSums(. != 0)) %>% add_column(Cell = cells)
  
  # Join all these four tibbles
  my_counts <- left_join(
    left_join(my_counts_e, my_counts_i, by = "Cell"),
    left_join(my_counts_l, my_counts_u, by = "Cell"),
    by = "Cell"
  ) %>%
    dplyr::select(Cell, n_e, n_i, n_l, n_u)
  
  return(my_counts)
}

# Create counts of viral UMIs
# robj = Seurat R object of infected cells
viral_counts_tibble_umis <- function(robj) {
  # Get Assay Counts and class of MVA genes
  mycounts <- as.matrix(robj[["RNA"]]@counts)

  # Keep counts of viral genes only
  is_mva <- robj@assays$RNA@meta.features$is_viral
  features <- as.vector(rownames(robj@assays$RNA@counts))
  keep <- features[which(is_mva == 1)]
  my_counts <- as_tibble(t(mycounts[keep,])) %>% 
    add_column("Cell" = colnames(robj[["RNA"]]@counts), .before = 1)
  # Determine the class of each gene
  early_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "early",])
  int_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "interm",])
  late_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "later",])
  unk_g <- rownames(robj[["RNA"]]@meta.features[robj[["RNA"]]@meta.features$gene_class == "unknown",])
  
  # Create sub tibble of early, int, late and unknown genes
  cells <- my_counts$Cell
  my_counts_e <- my_counts[c(early_g)] %>% mutate(n_e = rowSums(.)) %>% add_column(Cell = cells)
  my_counts_i <- my_counts[c(int_g)] %>% mutate(n_i = rowSums(.)) %>% add_column(Cell = cells)
  my_counts_l <- my_counts[c(late_g)] %>% mutate(n_l = rowSums(.)) %>% add_column(Cell = cells)
  my_counts_u <- my_counts[c(unk_g)] %>% mutate(n_u = rowSums(.)) %>% add_column(Cell = cells)
  
  # Join all these four tibbles
  my_counts <- left_join(
    left_join(my_counts_e, my_counts_i, by = "Cell"),
    left_join(my_counts_l, my_counts_u, by = "Cell"),
    by = "Cell"
  ) %>%
    dplyr::select(Cell, n_e, n_i, n_l, n_u)
  
  return(my_counts)
}


# Create tibble for barplot
# my_tibble = tibble with genes of interest
# arrange_col = variable to arrange x axis by
# my_counts = viral counts
# robj = seurat Object
# all = add cellular genes or not
# add_unk = add unknown genes or not
# threshold = detection threshold
barplot_tibble <- function(my_tibble, arrange_col, my_counts, robj, all = FALSE, add_unk = FALSE, threshold = 2) {
  # The tibble is different if we choose to include in the
  # plot the unanottated genes or not
  if(all){
    denominator <- rlang::parse_expr("nGene_MVA+nGene_cellular")
  }
  else if(add_unk){
    denominator <- rlang::parse_expr("n_l+n_e+n_i+n_u")
  }
  else{
    denominator <- rlang::parse_expr("n_l+n_e+n_i")
  }
  
  # Create the tibble with % of "class" gene 
  # among all (with or without unannotated) viral genes
  
  arrange_col <- rlang::parse_expr(arrange_col)
  barplot_tibble <- left_join(my_tibble, my_counts, by="Cell") %>%
    add_column(nGene_MVA = robj$nGene_MVA, .after=3) %>%
    add_column(nGene_cellular = robj$nGene_cellular, .after=3) %>%
    add_column(log10_nUMI_Cellular = log10(robj$nUMI_cellular), .after=3) %>%
    add_column(log10_nUMI_MVA = log10(robj$nUMI_MVA), .after=3) %>%
    add_column(ratioMVA = robj$nUMI_MVA_norm, .after=3) %>%
    add_column(mcherry_n = robj$mcherry_n, .after=3) %>%
    dplyr::filter(log10_nUMI_MVA > threshold) %>%
    mutate(denominator = !!denominator) %>%
    mutate(percent_late = round(n_l*100/(denominator))) %>%
    mutate(percent_inter = round(n_i*100/(denominator))) %>%
    mutate(percent_early = round(n_e*100/(denominator))) %>%
    mutate(percent_unk = round(n_u*100/(denominator))) %>%
    dplyr::select(Cell, percent_early, percent_inter, percent_late, percent_unk, nGene_cellular,
                  nGene_MVA, log10_nUMI_Cellular, log10_nUMI_MVA, ratioMVA, mcherry_n) %>% 
    mutate_all(funs(replace(., is.na(.), 0))) %>%
    arrange(!!arrange_col) %>%
    filter(percent_early != 0)

  if(all){
    barplot_tibble <- barplot_tibble %>%
      mutate(denominator = !!denominator) %>%
      mutate(percent_cellular = round(nGene_cellular*100/(denominator))) %>%
      dplyr::select(-percent_unk) %>%
      dplyr::select(-denominator)
  }
  else if(!add_unk){
    barplot_tibble <- barplot_tibble %>% dplyr::select(-percent_unk)
  }
  
  # Put the tibble in long format so the plot has expected aspect
  barplot_tibble <- pivot_longer(barplot_tibble, 
                                 -c("Cell", "nGene_cellular",
                                    "nGene_MVA", "log10_nUMI_Cellular", 
                                    "log10_nUMI_MVA", "ratioMVA", "mcherry_n")
                                 )
  barplot_tibble$x_axis <- factor(barplot_tibble$Cell, levels=unique(barplot_tibble$Cell))
  
  return(barplot_tibble)
}

# Barplot representing the proportion of each class of MVA gene
# barplot_tibble = tibble with informations for the barplot
# outpath = output path 
# filename = name of the output file
# sort = sort x axis by
# title = title of the plot 
# filetype = png or pdf
class_mva_barplot <- function(barplot_tibble, outpath, filename, sort, title, filetype) {
  p <- ggplot(barplot_tibble, aes(x_axis, value, fill = name)) +
    geom_bar(position = "fill", stat = "identity", width = 2) +
    scale_fill_manual(values=c("#D55E00", "#CC79A7", "#0072B2", "#009E73", "#FFFFFF")) +
    scale_y_reverse() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    labs(y = "Proportion of the gene class in viral UMIs", 
         x = paste0("Cells sorted by ", sort),
         title = title)
  
  if(filetype == "png"){
    png2(file.path(outpath, paste0(filename, ".png")))
    print(p)
    dev.off()
  }
  if(filetype == "pdf"){
    pdf(file.path(outpath, paste0(filename, ".pdf")))
    print(p)
    dev.off()
  }
  
}

# Geom_line with ratio MVA
# barplot_tibble = tibble with infos
# outpath = output path
# filename = name of the file
# filetype = png or pdf
# y_axis = variable in y
# ylab = name of the y axis
# x_axis = variable in x
# xlab = name of the x axis
# title = title of the plot
# stat = regress out or not (boolean)
# color = color of the line
# points = show data point (boolean)
mva_line <- function(barplot_tibble, outpath, filename, filetype, 
                     y_axis, ylab, x_axis = "x_axis", xlab, title, stat = FALSE, color = "black", points = TRUE) {
  barplot_tibble$group = 1
  p <- ggplot(barplot_tibble, aes_string(x_axis, y_axis, group = "group")) +
    scale_y_continuous(position = "right") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") +
    labs(y = ylab, x = xlab, 
         title = title)
  
  if(filetype == "png"){
    png2(file.path(outpath, paste0(filename, ".png")))
  }
  if(filetype == "pdf"){
    pdf(file.path(outpath, paste0(filename, ".pdf")))
  }
  
  # Fit a smoothing line
  if(stat){
    p1 <- p +
      stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, label.x.npc="left") +
      geom_smooth(method="loess", color=color)
    if(points){
      p1 <- p1 + geom_point()
    }
    print(p1)
  }
  # Or just print a line
  else{
    p <- p +
      geom_line(color=color)
    print(p)
  }
  dev.off()
}

# Regression mcherry with a variable of interest
# my_tibble = tibble with genes of interest
# x_axis, y_axis = variable in x, y
# xlab, ylab = name of the x and y axis
# title = title of the plot
# reg = method for the regression
# logscale = put the plot in log scale
reg_mcherry <- function(my_tibble, x_axis, y_axis, xlab, ylab, title, reg = "lm", logscale = FALSE){
  p <- ggplot(my_tibble, aes_string(x_axis, y_axis)) +
    geom_point() + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, label.x.npc="left") +
    geom_smooth(method=reg, color="black") +
    labs(title = title, xlab=xlab, ylab=ylab)
  
  if(logscale){
    p <- p + scale_y_continuous(trans="log") + scale_x_continuous(trans="log")
  }
  
  print(p)
}

# Return correlation table of mcherry with the genes
# seurat_o = Seurat object
# list_genes = list of genes of interest
# corr = feature to correlate to other features (gene or other feature)
# method = correlation method
corr_table <- function(seurat_o, list_genes, corr = "nUMI_MVA_norm", method = "spearman"){
  seurat_o[["RNA"]]@data["mcherry",] <- seurat_o$mcherry_n
  seurat_o$nUMI_cellular <- log(seurat_o$nUMI_cellular + 1)*10000
  cor_matrix <- apply(GetAssayData(seurat_o), 1, function(m) cor(x = m, y = seurat_o@meta.data[[corr]], method = method))
  return(cor_matrix[list_genes])
}

# Violin Plots with features of interest
# seurat_o = Seurat object
# features = list of features of interest
# cols = color vector
violin <- function(seurat_o, features, cols = NULL){
  seurat_o$ratioMVA <- seurat_o$nUMI_MVA/seurat_o$nCount_RNA
  for(f in features){
    if(f %in% c("nUMI_MVA", "nUMI_cellular")){log = TRUE}else{log = FALSE}
    v <- VlnPlot(seurat_o, features = f, group.by = "clusters", pt.size = 0.5, cols = cols, log = log)
    print(v)
  }
}

# Create a reduced expression matrix
# mat = expression matrix
# DE_file = dataframe, result from differential analysis
# annotation_cluster = dataframe of cluster annotation
# top = top number of genes to display
# genes = take only up-regulated genes into account (up), down-regulated (down) or 
# filter based on absolute logFC (all)
matrix_expression <- function(mat, DE_file, annotation_cluster, top = 50, genes = "all") {
  if(genes == "up"){
    tab <- read_tsv(DE_file) %>% filter(avg_logFC > 0) %>% top_n(top, avg_logFC)
  } else if(genes == "down"){
    tab <- read_tsv(DE_file) %>% filter(avg_logFC < 0) %>% top_n(top, desc(avg_logFC))
  } else if(genes == "all"){
    tab <- read_tsv(DE_file) %>% top_n(top, abs(avg_logFC))
  } else{
    stop("Invalid values of genes arguments, must be all, up or down")
  }
  mat <- mat[tab$Gene, rownames(annotation_cluster)]
  return(mat)
}

# Apply hierarchical clustering on cells of a matrix
# mat = expression matrix
# annotation_cluster = dataframe, cluster annotation
cluster_cells_HM <- function(mat, annotation_cluster) {
  # Cluster cells INSIDE function clusters (Uninfected, 0, 1, 2)
  duninf <- dist(t(mat[, rownames(annotation_cluster[annotation_cluster$Cluster == "Uninfected", ,drop = FALSE])]))
  clustuninf <- hclust(duninf)
  d0 <- dist(t(mat[, rownames(annotation_cluster[annotation_cluster$Cluster == 0, ,drop = FALSE])]))
  clust0 <- hclust(d0)
  d1 <- dist(t(mat[, rownames(annotation_cluster[annotation_cluster$Cluster == 1, ,drop = FALSE])]))
  clust1 <- hclust(d1)
  d2 <- dist(t(mat[, rownames(annotation_cluster[annotation_cluster$Cluster == 2, ,drop = FALSE])]))
  clust2 <- hclust(d2)
  # Take order of these clusters and create one vector for ordering the matrix
  oduninf <- clustuninf$order
  od0 <- clust0$order
  od1 <- clust1$order
  od2 <- clust2$order
  od = c(which(annotation_cluster == "Uninfected")[oduninf], which(annotation_cluster == 0)[od0], 
         which(annotation_cluster == 1)[od1], which(annotation_cluster == 2)[od2])
  return(od)
}

# List the pathway files of interest
# filenames = name of the files
# end = suffix
list_pathways <- function(filenames, end) {
  # List all files satistying the "end" condition
  filelist <- filenames[unlist(lapply(filenames, function(x) endsWith(x, end)))]
  # Names for easier handling late
  names(filelist) <- c("C0", "C1", "C2")
  print(filelist)
  return(filelist)
}

# Create the tibble for plotting
# inpath = directories with pathwayf iles
# filelist = list of comparisons names
pathway_tibble <- function(inpath, filelist, padjust = 1e-2, count = 5){
  for(i in seq_along(filelist)){
    if(i == 1){
      # Read the first pathway table
      df <- read_tsv(file.path(inpath, filelist[i]))
      df$set = names(filelist)[i]
      df <- df %>% filter(p.adjust < padjust, Count >= count)
    }
    else{
      # Join the other pathway table with the previous one
      tab <- read_tsv(file.path(inpath, filelist[i]))
      tab$set = names(filelist)[i]
      tab <- tab %>% filter(p.adjust < padjust, Count >= count)
      df <- bind_rows(df, tab)
    }
  }
  return(df)
}

# Plot the recap pathway
# df = pathways tibble
# title = ttle of the plot
recap_pathways <- function(df, title) {
  g <- ggplot(df, aes(set, forcats::fct_inorder(Pathway))) + 
    #facet_grid(cols = vars(grid)) +
    geom_point(aes(size=Count, color=`p.adjust`)) +
    scale_color_gradient(low="black", high="gray") +
    scale_size_continuous(range = c(0.1, 3)) +
    labs(title = title) +
    xlab("Gene sets") +
    ylab("Pathways") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          text = element_text(size = 12),
          plot.title = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(g)
}

# Plot the pathways graph
# title = title of the graph
# prefix = MF, BP, CC or Kegg
plot_pathway <- function(indir, title, prefix = "MF", padjust = 1e-2, count = 5) {
  filenames <- list.files(file.path(indir, "DOWN"))
  filelist <- list_pathways(filenames, paste0(prefix, ".tsv"))
  down <- pathway_tibble(file.path(indir, "DOWN"), filelist, padjust, count)
  down$set <- unlist(lapply(down$set, function(x) paste0("DOWN ", x)))
  filenames <- list.files(file.path(indir, "UP"))
  filelist <- list_pathways(filenames, paste0(prefix, ".tsv"))
  up <- pathway_tibble(file.path(indir, "UP"), filelist, padjust, count)
  up$set <- unlist(lapply(up$set, function(x) paste0("UP ", x)))
  df <- bind_rows(down, up)
  # Shorten the name of some pathways
  df$Description <- replace(df$Description, 
                            df$Description=="exonuclease activity, active with either ribo- or deoxyribonucleic acids and producing 5'-phosphomonoesters",
                            c("exonuclease activity"))
  df$Description <- replace(df$Description, 
                            df$Description=="negative regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
                            c("negative regulation of adaptive immune response based on recombination of immune receptors"))
  df$Description <- replace(df$Description, 
                            df$Description=="positive regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
                            c("positive regulation of adaptive immune response based on recombination of immune receptors"))
  df$Description <- replace(df$Description, 
                            df$Description=="regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
                            c("regulation of adaptive immune response based on recombination of immune receptors"))
  df$Description <- replace(df$Description, 
                            df$Description=="adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
                            c("adaptive immune response based on recombination of immune receptors"))
  # Sort by set and Count for better vizualisation
  df <- df %>%
    group_by(set) %>%
    top_n(100, Count) %>%
    arrange(desc(set), desc(Count)) %>%
    ungroup() %>%
    mutate(Pathway = factor(Description))
  dotplot <- recap_pathways(df, title)
  return(dotplot)
}

# Gene set over representation analysis
# list_DE_files <- list of tables of DE results
# list_names = list of cluster names
# file_gmt = path to gmt file
# title = plot title
enrich_geneset <- function(list_DE_files, list_names, file_gmt, title,
                           pvalue = 1e-2, count = 5) {
  list_names = c("C0", "C1", "C2")
  geneset <- read.gmt(file_gmt)
  df <- tibble()
  for(i in seq_along(list_DE_files)){
    degs <- read_tsv(list_DE_files[[i]])
    up_degs <- degs %>% filter(avg_logFC > 0)
    down_degs <- degs %>% filter(avg_logFC < 0)
    up_enrich <- enricher(up_degs$Gene, TERM2GENE = geneset, pvalueCutoff = pvalue)
    up_enrich <- as_tibble(up_enrich) %>% filter(Count >= count)
    down_enrich <- enricher(down_degs$Gene, TERM2GENE = geneset, pvalueCutoff = pvalue)
    down_enrich <- as_tibble(down_enrich) %>% filter(Count >= count)
    up_enrich$set <- paste0("UP_", list_names[[i]])
    down_enrich$set <- paste0("DOWN_", list_names[[i]])
    enrich <- bind_rows(up_enrich, down_enrich)
    df <- bind_rows(df, enrich)
  }
  df <- df %>% dplyr::rename(Pathway = ID)
  print(df)
  return(df)
}

# Differential expression of one cluster against each other one individually
# used to put stars in the violin plot
# infected = Seurat RDS object of infected cells
# uninfected = Seurat R object of uninfecred cells
# outpath = output path 
# result = list of arguments from json file
pairwise_DE <- function(infected, uninfected, outpath, result){
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
  #List of genes of interest for the stars
  list_genes = result$figures$list_vln
  
  # Merge infected and uninfected
  seurat_o <- merge_inf_uninf(infected, uninfected, cluster_sol)
  
  # Create the color vector for the heatmap
  colors <- rev(RColorBrewer::brewer.pal(8, "RdBu"))
  # Create annotation of clusters
  clusters <- levels(Idents(seurat_o))
  # Create output path
  dir.create(file.path(outpath, "Pairwise_DE"))
  
  # Add nUMI_cellular, nUMI_MVA and ratioMVA to the matrix
  matcell <- matrix(log(seurat_o$nUMI_cellular/sum(seurat_o$nUMI_cellular)*10000+1))
  rownames(matcell) <- names(seurat_o$nUMI_cellular)
  colnames(matcell) <- c("nUMI_cellular")
  seurat_o@assays$RNA@data <- rbind(seurat_o@assays$RNA@data, t(matcell))
  matcell <- matrix(log(seurat_o$nUMI_MVA/sum(seurat_o$nUMI_MVA)*10000+1))
  rownames(matcell) <- names(seurat_o$nUMI_MVA)
  colnames(matcell) <- c("nUMI_MVA")
  seurat_o@assays$RNA@data <- rbind(seurat_o@assays$RNA@data, t(matcell))
  matcell <- matrix(seurat_o$nUMI_MVA_norm)
  rownames(matcell) <- names(seurat_o$nUMI_MVA_norm)
  colnames(matcell) <- c("ratioMVA")
  seurat_o@assays$RNA@data <- rbind(seurat_o@assays$RNA@data, t(matcell))
  
  ref <- clusters
  recap_tib <- tibble()
  for(i in seq_along(clusters)[-length(clusters)]){
    ## Differential expression
    # Remove scientific annotation for easy comparisons
    options(scipen=999)
    ref <- ref[-1]
    for(j in seq_along(ref)){
      markers <- sc_diffexp(seurat_o, clusters[i], ref[j], logFC, pvalue, test, min.pct)
      # Put scientific annotation back
      options(scipen=0)
      # Sort by absolute logFC and write genes
      marker_tibble <- as_tibble(markers) %>%
        add_column(Gene = rownames(markers)) %>%
        arrange(desc(abs(avg_logFC))) %>%
        filter(Gene %in% list_genes) %>%
        add_column(comparison = paste0(clusters[i], "_VS_", ref[j]))
      # Write individual results
      write_tsv(marker_tibble, file.path(outpath, "Pairwise_DE", paste0(clusters[i], "_VS_", ref[j], ".tsv")))
      recap_tib <- bind_rows(recap_tib, marker_tibble)
    }
  }
  # Write recap of DE tables
  write_tsv(recap_tib, file.path(outpath, "Pairwise_DE", "Pairwise_DE.tsv"))
}


# Create the barplots, FeaturePlot and VlnPlot of the paper
# robj = Seurat infected object
# uninfected = Seurat uninfected object
# outpath = output path
# result = from args json file
main_paper_figures <- function(robj_file, uninfected_file, outpath, indir, result) {
  
  # Key of the reduction of interest
  red = result$figures$red
  # Number of dimensions to plot
  d = result$figures$d
  # List of genes of interest
  list_genes = result$figures$list_genes
  # Threshold for filtering tibble
  threshold = result$figures$threshold
  # List of features of interest for VlnPlot
  list_vln = result$figures$list_vln
  # Cluster of interest 
  seurat_cluster = result$figures$seurat_cluster
  # List of DE files
  DE_files <- result$figures$DE_files
  # Top number of genes to display in heatmaps
  top <- result$figures$top
  # wether to use scaled value or not
  scaled <- as.logical(result$figures$scaled)
  # list of genes of interest for the heatmap
  genes_heatmap <- result$figures$list_HM
  # List of genes for the correlation barplot
  list_corr <- result$figures$list_corr

  dir.create(outpath)
  robj <- readRDS(robj_file)
  rep1 <- SubsetData(robj, subset.name = "orig.ident", accept.value = "MVA-rep1")
  rep2 <- SubsetData(robj, subset.name = "orig.ident", accept.value = "MVA-rep2")
  
  list_genes <- c("TNF", "IL6")
  
  my_tibble_both <- create_tibble(robj, list_genes, TRUE, FALSE, 1)
  my_tibble_1 <- create_tibble(rep1, list_genes, TRUE, FALSE, 1)
  my_tibble_2 <- create_tibble(rep2, list_genes, TRUE, FALSE, 1)
  
  #################
  #### 1) UMAP ####
  #################
  
  for(genes in c(list_genes, "nUMI_MVA_norm", "nUMI_cellular")){
    expr_plots(robj, red, d, genes, outpath)
  }
  
  #####################
  #### 2) BARPLOTS ####
  #####################
  
  my_counts_both <- viral_counts_tibble_umis(robj)
  my_counts_1 <- viral_counts_tibble_umis(rep1)
  my_counts_2 <- viral_counts_tibble_umis(rep2)
  
  # BOTH SAMPLES
  
  # Sort by mcherry norm
  before1 <- dim(filter(my_counts_both, startsWith(Cell, "MVA-rep1")))[1]
  before2 <- dim(filter(my_counts_both, startsWith(Cell, "MVA-rep2")))[1]
  final_tibble <- barplot_tibble(my_tibble_both, "mcherry_n", my_counts_both, robj, FALSE, FALSE, threshold)
  final_tibble <- dplyr::filter(final_tibble, mcherry_n != 0)
  after1 <- dim(filter(final_tibble, startsWith(Cell, "MVA-rep1")))[1]/3
  after2 <- dim(filter(final_tibble, startsWith(Cell, "MVA-rep2")))[1]/3
  
  # Titles
  title_both <- paste("Proportion of viral classes in viral UMIs in mcherry positive infected cells \nexpressing more than",
                      threshold, "log 10 UMI MVAs.\n", 
                      after1, "out of", before1, "cells represented for replicate 1.\n",
                      after2, "out of", before2, "cells represented for replicate 2.\n",
                      sep = " "
  )
  
  # Barplot of classes of genes
  class_mva_barplot(final_tibble, outpath, "MVA-all_sorted_by_mcherry", "normalized value of mcherry", 
                    title_both, "pdf")
  
  # Regression (threshold = 0 cause we want to keep all mcherry positive cells)
  my_tibble_both <- create_tibble(robj, list_genes, TRUE, FALSE, 0)
  final_tibble <- barplot_tibble(my_tibble_both, "mcherry_n", my_counts_both, robj, FALSE, FALSE, 0)
  final_tibble <- dplyr::filter(final_tibble, mcherry_n != 0)
  mva_line(final_tibble, outpath,
           "MVA-all_mcherry_sorted_by_mcherry_RED", "pdf",
           "mcherry_n", "Normalized value of mcherry", 
           "x_axis", "Cells sorted by value of mcherry", "Mcherry distribution", FALSE, "red")
  mva_line(final_tibble, outpath,
           "MVA-all_mcherry_sorted_by_mcherry_BLACK", "pdf"
           , "mcherry_n", "Normalized value of mcherry", 
           "x_axis", "Cells sorted by value of mcherry", "Mcherry distribution", FALSE)
  mva_line(final_tibble, outpath,
           "MVA-all_ratioMVA_sorted_by_mcherry_REGRESSION", "pdf",
           "ratioMVA", "Ratio of MVA in the cell", 
           "x_axis", "Cells sorted by value of mcherry", "Ratio regression", TRUE, "black", FALSE)
  mva_line(final_tibble, outpath,
           "MVA-all_ratioMVA_sorted_by_mcherry_REGRESSION_WITH_POINTS", "pdf",
           "ratioMVA", "Ratio of MVA in the cell", 
           "x_axis", "Cells sorted by value of mcherry", "Ratio regression", TRUE, "red")


  ###########################
  ### 3) REGRESSION TESTS ###
  ###########################
  
  my_tibble_both$mcherry_raw <- robj[["RNA"]]@counts["mcherry",]
  my_tibble_both <- my_tibble_both %>% mutate(ratioMVA = log((nUMI_MVA/(nUMI_Cellular+nUMI_MVA))*10000+1)) %>%
    mutate(ratioMVA_raw = nUMI_MVA/(nUMI_Cellular+nUMI_MVA))
  my_tibble_both$mcherry_seurat <- robj$mcherry_n
  
  reg_mcherry(my_tibble_both, "mcherry_raw", "nUMI_MVA", 
              "Raw values of mcherry", "Number of UMI MVA", 
              "Linear regression of viral UMIs, log-scale", "lm", logscale = TRUE)
  reg_mcherry(my_tibble_both, "mcherry_seurat", "ratioMVA", 
              "Normalized values of mcherry", "Ratio of MVA, log-transformed", 
              "Linear regression of viral ratio", "lm")
  
  
  ###########################
  ##### 4) VIOLIN PLOTS #####
  ###########################
  
  uninfected <- get(load(uninfected_file))
  merged <- merge(robj, uninfected)
  # Create a metadata with the clusters + a "Uninfected" cluster
  times <- length(merged@meta.data[[seurat_cluster]]) - length(robj@meta.data[[seurat_cluster]])
  clusters <- c(as.vector(robj@meta.data[[seurat_cluster]]), rep("Uninfected", times))
  merged <- AddMetaData(merged, clusters, "clusters")
  
  cols = c(scales::hue_pal()(3), "#BCBCBC")
  pdf(file.path(outpath, "violinplots.pdf"))
  violin(merged, list_vln, cols)
  dev.off()
  
  #############################
  ### 5) CORRELATION TABLES ###
  #############################
  
  ## Correlations with ratio of MVA in infected cells
  mat <- corr_table(robj, list_corr, method = "pearson")
  tib <- tibble(Gene=names(mat), Pearson=mat)
  tib <- mutate(tib, x_axis = factor(Gene, levels=list_corr))
  write_tsv(tib, file.path(outpath, "correlations_pearsons.tsv"))
  g <- ggplot(tib, aes(x=x_axis, y=Pearson)) +
    geom_bar(stat="identity", fill="black") + 
    theme(axis.text.x = element_text(angle = 60),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") +
    ylab("Pearson correlation") +
    xlab("Feature") +
    ylim(c(-0.5,1))
  pdf(file.path(outpath, "barplots_correaltion_pearson.pdf"))
  print(g)
  dev.off()
  
  #############################
  ### 6) EXPRESSION HEATMAP ###
  #############################
  
  # Merge infected and uninfected
  seurat_o <- merge_inf_uninf(robj_file, uninfected_file, seurat_cluster)
  clustcomp <- Idents(seurat_o)
  # Create color vector for the heatmap
  clusters <- levels(Idents(seurat_o))[-length(levels(Idents(seurat_o)))]
  # Create annotation of clusters
  colors <- rev(RColorBrewer::brewer.pal(10, "RdBu"))
  clustercol <- c("#BCBCBC", scales::hue_pal()(length(clusters)))
  annotation_cluster <- cluster_annotation(seurat_o, "clustcomp")
  annotation_cluster$Cluster <- factor(annotation_cluster$Cluster, levels=c("Uninfected",0,1,2))
  annotation_cluster <- annotation_cluster[order(annotation_cluster$Cluster), , drop=FALSE]
  # Colors for annotation
  clustcomp <- factor(clustcomp, levels=c("Uninfected",0,1,2))
  ann_color <- list(
    # Seurat colors
    Cluster = setNames(clustercol, levels(as.factor(clustcomp)))
  )
  
  # Expression heatmap
  if(scaled){
    seurat_o <- ScaleData(seurat_o, features=rownames(seurat_o))
    mat <- as.matrix(seurat_o[["RNA"]]@scale.data)
  } else {
    mat <- as.matrix(seurat_o[["RNA"]]@data)
  }
  
  # Get the list of DE results files
  list_DE_files <- unlist(lapply(DE_files, function(x) file.path(indir, x)))
  
  # Create gaps for heatmap
  gaps <- as_tibble(annotation_cluster) %>% group_by(Cluster) %>% tally()
  my_gaps <- c(gaps$n[1]+1,gaps$n[2]+gaps$n[1]+1,gaps$n[3]+gaps$n[2]+gaps$n[1]+1)
  myBreaks = seq(-2.5,2.5,0.5)
  
  # Heatmap of Top 25 of each gene
  mat_c0 <- matrix_expression(mat, list_DE_files[[1]], annotation_cluster, top, "up")
  mat_c1 <- matrix_expression(mat, list_DE_files[[2]], annotation_cluster, top, "up")
  mat_c2 <- matrix_expression(mat, list_DE_files[[3]], annotation_cluster, top, "up")
  mat_top25 <- rbind(mat_c0, mat_c1, mat_c2)
  cell_order <- cluster_cells_HM(mat_top25, annotation_cluster)
  topHM <-  pheatmap::pheatmap(mat_top25[,cell_order], color = colors, annotation_colors = ann_color, 
                     annotation_col = annotation_cluster, fontsize_row = 9,
                     show_rownames = TRUE, show_colnames = FALSE, gaps_col = my_gaps,
                     cluster_rows = FALSE, cluster_cols = FALSE, breaks = myBreaks,
                     main = paste0("Scaled expression values"))
  pdf(file.path(outpath, "Top_25_UP_DEG.pdf"), width=6, height=10)
  print(topHM)
  dev.off()
  
  png(file.path(outpath, "Top_25_UP_DEG.png"), height=900, width=600)
  print(topHM)
  dev.off()
  
  
  # Heatmap for features of interest (the one in the VlnPlots)
  vUMI <- as.matrix(seurat_o$nUMI_MVA)
  vUMI <- t(vUMI)
  rownames(vUMI) <- "vUMI"
  cUMI <- as.matrix(seurat_o$nUMI_cellular)
  cUMI <- t(cUMI)
  rownames(cUMI) <- "cUMI"
  ratio <- as.matrix(seurat_o$nUMI_MVA_norm)
  ratio <- t(ratio)
  rownames(ratio) <- "vUMI/tUMI"
  
  mat <- as.matrix(seurat_o[["RNA"]]@data)
  mat2 <- rbind(mat, vUMI, cUMI, ratio)
  mat2 <- mat2[genes_heatmap, rownames(annotation_cluster)]
  cell_order <- cluster_cells_HM(mat2, annotation_cluster)

  vlnHM <- pheatmap::pheatmap(mat2[,cell_order], color = colors, annotation_colors = ann_color, 
                              annotation_col = annotation_cluster, scale = "row",
                              show_rownames = TRUE, show_colnames = FALSE, gaps_col = my_gaps,
                              gaps_row = c(4, 8, 12, 16),
                              cluster_rows = FALSE, cluster_cols = FALSE, breaks = myBreaks,
                              main = paste0("Scaled expression values"))
  pdf(file.path(outpath, "HM_Gene_interest.pdf"), height=6, width=6)
  print(vlnHM)
  dev.off()
  
  png(file.path(outpath, "HM_Gene_interest.png"), height=600, width=600)
  print(vlnHM)
  dev.off()
  
  ###########################
  ### 7) PATHWAYS DOTPLOT ###
  ###########################

  # Pathway dotplot for GO and Kegg
  pdf(file.path(outpath, "Pathways_dotplot_MF.pdf"), height=4, width=7.5)
  print(plot_pathway(indir, "GO MF Pathways in our clusters", prefix = "MF"))
  dev.off()
  pdf(file.path(outpath, "Pathways_dotplot_BP.pdf"), height=13, width=10)
  print(plot_pathway(indir, "GO BP Pathways in our clusters", prefix = "BP"))
  dev.off()
  pdf(file.path(outpath, "Pathways_dotplot_CC.pdf"), height=6, width=7)
  print(plot_pathway(indir, "GO CC Pathways in our clusters", prefix = "CC"))
  dev.off()
  pdf(file.path(outpath, "Pathways_dotplot_Kegg.pdf"), height=8, width=7)
  print(plot_pathway(indir, "Kegg Pathways in our clusters", prefix = "kegg"))
  dev.off()
  
  # Hallmark gene set enrichment
  hallmark <- enrich_geneset(list_DE_files, c("C0", "C1", "C2"), 
                             "../../data/GSEA_GeneSets/hallmark.all.v7.1.symbols.gmt", 
                             "Hallmark gene sets")
  
  hallmark$Pathway <- gsub("nfkb", "NFKB", gsub("tnfa", "TNFA", gsub("_", " ", R.utils::capitalize(tolower(hallmark$Pathway)))))
  hallmark_plot <- recap_pathways(hallmark, "Hallmark gene sets")
  
  # C7 gene set enrichment
  imm <- enrich_geneset(list_DE_files, c("C0", "C1", "C2"), 
                        "../../data/GSEA_GeneSets/c7.immunological.all.v7.1.symbols.gmt", 
                        "Immunological gene sets", count = 15)
  imm_plot <- recap_pathways(imm, "Immunological gene sets")
  
  # Plot dotplots to PDF
  pdf(file.path(outpath, "Dotplot_Hallmark.pdf"), height=4.5, width=5.5)
  print(hallmark_plot)
  dev.off()
  
  pdf(file.path(outpath, "Dotplot_C7_Immunological.pdf"), height=11.5, width=9.5)
  print(imm_plot)
  dev.off()
  
  ######################
  ### 8) PAIRWISE DE ###
  ######################
  
  pairwise_DE(robj_file, uninfected_file, outpath, result)
  
  ######################
  ### 9) SILHOUETTE  ###
  ###    AND BARPLOT ###
  ######################
  
  robj$clusters <- robj@meta.data[[seurat_cluster]]
  clust_t <- clust_tibble(robj)
  barplot <- barplot_clusters(clust_t, "sample") +
    theme_bw()
  pdf(file.path(outpath, "Sample_barplot.pdf"))
  print(barplot)
  dev.off()
  
  # Compute Silhouette score
  robj$seurat_clusters <- robj$clusters
  silscore <- computeSilhouette(robj, "mnn1_")
  
  # Plot Silhouette score
  Silhouetteplot <- plotSilhouette(silscore)
  
  pdf(file.path(outpath, "Silhouettescore.pdf"))
  print(Silhouetteplot)
  dev.off()
  
}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  main_paper_figures(robj_file, uninfected_file, outpath, indir, result)
}


