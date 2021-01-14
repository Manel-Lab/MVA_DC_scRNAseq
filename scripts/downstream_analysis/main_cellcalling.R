## ---------------------------
##
## Purpose of script: Cell calling using EmptyDrops
##
## Author: Tarek Gharsalli
##
## Date Created: 2020-07-03
##
## ---------------------------

## Arguments

if (sys.nframe() == 0){
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 9) {
    cat("Usage: Rscript cellCalling.R <matrix_dir_raw>
  <matrix_dir_filt> <out_dir> <max-umi_value>")
    q()
  }
  
  MATRIX_DIR_RAW <- args[1]
  MATRIX_DIR_FILT <- args[2]
  OUT_DIR <- args[3]
  MAX_UMI_MODEL<-as.numeric(args[4])
  FDR <- as.numeric(args[5])
  MIN_UMI_COUNT <- as.numeric(args[6])
  scale_factor <- as.numeric(args[7])
  MVA_genes_file <- args[8]
  MVA_class <- args[9]
  OUTFILE <- args[10]
}

## ---------------------------

## Packages

library(Seurat)
library(DropletUtils)

## ------ FUNCTIONS ------- ##

# Cell calling using EmptyDrops
# MATRIX_DIR_RAW : path to CellRanger raw matrix directory
# MATRIX_DIR_FILT : path to CellRanger filtered matrix directory
# OUT_DIR : output directory
# MAX_UMI_MODEL : lower bound of total UMI count for EmptyDrops
# FDR : false disovery rate threshold for EmptyDrops results
# MIN_UMI_COUNT : minimum number of UMIs to retain a cell
# scale_factor : scaling factor for log-normalization
# MVA_genes_file : csv indicating the names of the genes of VACV virus
# MVA_class : tsv indicating the class of the viral genes (early, intermediary, late)
# OUTFILE : name of the output file
cellCalling <- function(MATRIX_DIR_RAW, MATRIX_DIR_FILT, OUT_DIR, MAX_UMI_MODEL, 
                        FDR, MIN_UMI_COUNT, scale_factor, MVA_genes_file, MVA_class, OUTFILE) {
  set.seed(269)
  dir.create(OUT_DIR)
  
  # Create a Seurat object of raw matrix
  object_raw.data <- Read10X(data.dir = MATRIX_DIR_RAW)
  object_raw <- CreateSeuratObject(counts = object_raw.data, project = "DC_MVA")
  barcodes_cr_all <- colnames(GetAssayData(object = object_raw))
  
  # Create a Seurat object of filtered matrix
  object_filt.data <- Read10X(data.dir = MATRIX_DIR_FILT)
  object_filt <- CreateSeuratObject(counts = object_filt.data, project = "DC_MVA")
  
  # Create an emptyDrop object of raw matrix
  sce <- read10xCounts(MATRIX_DIR_RAW)
  ed_raw <- emptyDrops(counts(sce), lower = MAX_UMI_MODEL)
  barcodes_ed_all <- sce@colData$Barcode
  
  # Reading the number of UMI and genes
  nUMI <-  object_raw$nCount_RNA 
  nGene <- object_raw$nFeature_RNA
  
  # Determine gene expression
  
  # Creating a matrix of all gene names (cell + MVA)
  genes <- GetAssayData(object = object_raw)
  genes <- genes@Dimnames
  genes <- genes[[1]]
  
  # Determine viral gene expression
  # Create an array of all MVA gene names
  MVA_genes <- read.csv2(MVA_genes_file)
  
  # Replace viral genes nomenclature in Seurat objects
  genes <- rownames(GetAssayData(object = object_raw, slot = "counts"))
  for(i in seq_along(genes)){
    if(genes[[i]] %in% MVA_genes$V){
      genes[[i]] <- as.character(MVA_genes[MVA_genes$V == genes[[i]],1])
    }
  }
  rownames(object_raw[["RNA"]]@counts) <- genes
  rownames(object_raw[["RNA"]]@data) <- genes

  # Extract expression matrix of MVA genes
  matrix_MVA_genes <- GetAssayData(object = object_raw, slot = "counts")[MVA_genes$MVA_genes, , drop = FALSE]
  # UMI count for MVA genes
  nUMI_MVA <- Matrix::colSums(matrix_MVA_genes) 
  # Number of detected MVA genes
  nGene_MVA <- Matrix::colSums(matrix_MVA_genes>0)
  # Assign cell names
  names(nUMI_MVA) <- names(nGene_MVA) <- colnames(object_raw@assays$RNA@data)
  # add nUMI_MVA to the raw_object
  object_raw$nUMI_MVA <- nUMI_MVA
  object_raw$nGene_MVA <- nGene_MVA
  # Reading and computing all barcodes
  # Read the barcode names from 10X and cellRanger
  barcodes <- colnames(GetAssayData(object = object_raw))
  barcodes_cr <- colnames(GetAssayData(object = object_filt))
  
  # Read the barcodes names from emptyDrop
  x <- (barcodes_ed_all %in% barcodes_cr_all ) & !(is.na(ed_raw$FDR))
  barcodes_ed <- barcodes_ed_all[x]
  
  # Determine the match between nUMI_MVA and barcodes_ed
  barcodes_ed <- barcodes_ed_all[x]
  # Create a vector that contain the same barcodes_ed and nUMI_MVA
  nUMI_MVA_and_ed <- nUMI_MVA[names(nUMI_MVA) %in% barcodes_ed]
  # Check how we recover cells with nUMI_MVA and FDR
  barcodes_ed_nUMI_MVA <- barcodes_ed[nUMI_MVA_and_ed > 100 & ed_raw$FDR[x] < FDR]
  # Check how we recover cells with nUMI_MVA, MIN_UMI_COUNT and FDR
  barcodes_ed_nUMI_MVA_MIN_UMI_COUNT <- barcodes_ed[((ed_raw$Total[x] > MIN_UMI_COUNT)|(nUMI_MVA_and_ed > 100)) & ed_raw$FDR[x] < FDR]
  barcodes_ed <- barcodes_ed_nUMI_MVA_MIN_UMI_COUNT
  # Compute the sets
  barcodes_inters <- barcodes_cr[barcodes_cr %in% barcodes_ed]
  barcodes_cr_not_ed <- barcodes_cr[!(barcodes_cr %in% barcodes_ed)]
  
  # Determine cellular gene expression
  nUMI_cellular <- nUMI - nUMI_MVA
  object_raw$nUMI_cellular <- nUMI_cellular
  nGene_cellular <- nGene - nGene_MVA
  object_raw$nGene_cellular <- nGene_cellular

  # Compute the sets
  barcodes_inters <- barcodes_cr[barcodes_cr %in% barcodes_ed]
  barcodes_cr_not_ed <- barcodes_cr[!(barcodes_cr %in% barcodes_ed)]
  
  # Determine mcherry expression
  mcherry <- GetAssayData(object = object_raw, slot = "counts")[genes %in% c("mcherry"),]  # extract mcherry gene
  names(mcherry) <- colnames(object_raw@assays$RNA@data)  # assign cell names
  object_raw$mcherry <- mcherry # add mcherry to object_raw
  
  GetAssayData(object = object_raw, slot = "counts")[genes %in% c("mcherry"),]
  
  # Determine the number of mcherry positive and mcherry negative in barcodes_cr and barcodes_ed
  # Extract barcodes with mcherry negative from the raw matrix
  mcherry_neg1 <- mcherry[mcherry < 1]
  # Select cells called by cr with mcherry negative
  mcherry_neg2 <- mcherry_neg1[names(mcherry_neg1) %in% barcodes_cr]
  # Select cells called by cr with mcherry negative
  mcherry_neg <- mcherry_neg1[names(mcherry_neg1) %in% barcodes_cr | names(mcherry_neg1) %in% barcodes_ed]
  
  # Extract barcodes with mcherry positive from the raw matrix
  mcherry_pos_cells1 <- mcherry[mcherry > 0]
  # Select cells called by cr with mcherry positive
  mcherry_pos_cells2 <- mcherry_pos_cells1[names(mcherry_pos_cells1) %in% barcodes_cr]
  mcherry_pos_cells <- mcherry_pos_cells1[names(mcherry_pos_cells1) %in% barcodes_cr | names(mcherry_pos_cells1) %in% barcodes_ed]
  
  # Determine the union between barcodes_cr and barcodes_ed
  barcodes_cr_union_ed <- union(barcodes_cr, barcodes_ed)
  
  # Create a new object containing the union between cr and ed
  object <- subset(object_raw, cells = barcodes_cr_union_ed) 
  # Limit based on UMI and total RNA threshold
  object <- subset(object, 
                   cells = names(object$nUMI_MVA[object$nUMI_MVA > 100 | object$nCount_RNA > MIN_UMI_COUNT]))
  
  
  # Normalize object_cr_and_ed 
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = scale_factor)
  
  # Normalize nUMI_MVA
  nUMI_MVA_norm <- log((object$nUMI_MVA/object$nCount_RNA)*scale_factor + 1)
  
  # Save in nUMI_MVA in metadata
  names(nUMI_MVA_norm) <- colnames(object@assays$RNA@data)
  object$nUMI_MVA_norm <- nUMI_MVA_norm
  
  # Normalize mcherry
  mcherry_n <- log((object$mcherry/object$nCount_RNA)*scale_factor + 1)
  
  # Save mcherry_n in metadata 
  names(mcherry_n) <- colnames(object@assays$RNA@data)
  object$mcherry_n <- mcherry_n

  # Get the class of MVA gene (early, intermediate, late)
  class_mva <- readr::read_tsv(MVA_class)
  
  # An array for each type of MVA gene
  # Early genes
  MVA_gene_ID_E <- dplyr::filter(class_mva, Class == "early")$Gene # extract early genes
  ## Intermediate  genes
  MVA_gene_ID_I <- dplyr::filter(class_mva, Class == "interm")$Gene # Intermediate genes
  ## Late genes
  MVA_gene_ID_L <- dplyr::filter(class_mva, Class == "later")$Gene # later genes
  
  # Create metadata with number of UMI counts for these classes in each cell
  object$early_genes <- colSums(GetAssayData(object, slot = "counts")[MVA_gene_ID_E,])
  object$interm_genes <- colSums(GetAssayData(object, slot = "counts")[MVA_gene_ID_I,])
  object$late_genes <- colSums(GetAssayData(object, slot = "counts")[MVA_gene_ID_L,])
  
  # Normalize gene class
  late_genes_n <- log((object$late_genes/object$nCount_RNA)*scale_factor + 1)
  object$late_genes_n <- late_genes_n
  early_genes_n <- log((object$early_genes/object$nCount_RNA)*scale_factor + 1)
  object$early_genes_n <- early_genes_n
  interm_genes_n <- log((object$interm_genes/object$nCount_RNA)*scale_factor + 1)
  object$interm_genes_n <- interm_genes_n
  
  # Save object_merged
  save(object, list = c("object"), file = file.path(OUT_DIR, OUTFILE))
  return(object)
}

cellCalling2 <- function(MATRIX_DIR_FILT, OUT_DIR, scale_factor, MIN_UMI_CELL, MIN_UMI_MVA, MVA_genes_file, MVA_class, OUTFILE) {
  
  set.seed(269)
  dir.create(OUT_DIR)
  
  # Create a Seurat object of filtered matrix
  object_filt.data <- Read10X(data.dir = MATRIX_DIR_FILT)
  object_filt <- CreateSeuratObject(counts = object_filt.data, project = "DC_MVA")
  
  # Reading the number of UMI and genes
  nUMI <-  object_filt$nCount_RNA 
  nGene <- object_filt$nFeature_RNA
  
  # Determine gene expression
  
  # Creating a matrix of all gene names (cell + MVA)
  genes <- GetAssayData(object = object_filt)
  genes <- genes@Dimnames
  genes <- genes[[1]]
  
  # Determine viral gene expression
  # Create an array of all MVA gene names
  MVA_genes <- read.csv2(MVA_genes_file)
  
  # Replace viral genes nomenclature in Seurat objects
  genes <- rownames(GetAssayData(object = object_filt, slot = "counts"))
  for(i in seq_along(genes)){
    if(genes[[i]] %in% MVA_genes$V){
      genes[[i]] <- as.character(MVA_genes[MVA_genes$V == genes[[i]],1])
    }
  }
  rownames(object_filt[["RNA"]]@counts) <- genes
  rownames(object_filt[["RNA"]]@data) <- genes
  
  # Extract expression matrix of MVA genes
  matrix_MVA_genes <- GetAssayData(object = object_filt, slot = "counts")[MVA_genes$MVA_genes, , drop = FALSE]
  # UMI count for MVA genes
  nUMI_MVA <- Matrix::colSums(matrix_MVA_genes) 
  # Number of detected MVA genes
  nGene_MVA <- Matrix::colSums(matrix_MVA_genes>0)
  # Assign cell names
  names(nUMI_MVA) <- names(nGene_MVA) <- colnames(object_filt@assays$RNA@data)
  # add nUMI_MVA to the filtered cell calling object
  object_filt$nUMI_MVA <- nUMI_MVA
  object_filt$nGene_MVA <- nGene_MVA
  
  # Determine cellular gene expression
  nUMI_cellular <- nUMI - nUMI_MVA
  object_filt$nUMI_cellular <- nUMI_cellular
  nGene_cellular <- nGene - nGene_MVA
  object_filt$nGene_cellular <- nGene_cellular
  
  ### script 3 start ###
  
  # Determine mcherry expression
  mcherry <- GetAssayData(object = object_filt, slot = "counts")[genes %in% c("mcherry"),]  # extract mcherry gene
  names(mcherry) <- colnames(object_filt@assays$RNA@data)  # assign cell names
  object_filt$mcherry <- mcherry # add mcherry to object_raw
  # Apply UMI threshold 
  # Keep cells with at least 1000 total RNA or 100 UMI MVA
  object_filt <- SubsetData(object_filt, 
                            cells = names(object_filt$nUMI_MVA[object_filt$nUMI_MVA > MIN_UMI_MVA | object_filt$nCount_RNA > MIN_UMI_CELL]))
  
  ### script 4 start ### ---> Normalize + nUMI_mva_norm and other stuff
  
  object <- object_filt
  # Normalize object_cr_and_ed 
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = scale_factor)
  
  # Normalize nUMI_MVA
  nUMI_MVA_norm <- log((object$nUMI_MVA/object@meta.data$nCount_RNA)*scale_factor +1)
  
  # Save in nUMI_MVA in metadata
  names(nUMI_MVA_norm) <- colnames(object@assays$RNA@data)
  object$nUMI_MVA_norm <- nUMI_MVA_norm
  
  # Normalize mcherry
  mcherry_n <- log((object$mcherry/object$nCount_RNA)*scale_factor +1)
  
  # Save mcherry_n in metadata 
  names(mcherry_n) <- colnames(object@assays$RNA@data)
  object$mcherry_n <- mcherry_n
  
  ### script 6_1 start ###
  
  # Get the class of MVA gene (early, intermediate, late)
  class_mva <- readr::read_tsv(MVA_class)
  
  # An array for each type of MVA gene
  # Early genes
  MVA_gene_ID_E <- dplyr::filter(class_mva, Class == "early")$Gene # extract early genes
  ## Intermediate  genes
  MVA_gene_ID_I <- dplyr::filter(class_mva, Class == "interm")$Gene # Intermediate genes
  ## Late genes
  MVA_gene_ID_L <- dplyr::filter(class_mva, Class == "later")$Gene # later genes
  
  # Create metadata with number of UMI counts for these classes in each cell
  object$early_genes <- colSums(GetAssayData(object, slot = "counts")[MVA_gene_ID_E,])
  object$interm_genes <- colSums(GetAssayData(object, slot = "counts")[MVA_gene_ID_I,])
  object$late_genes <- colSums(GetAssayData(object, slot = "counts")[MVA_gene_ID_L,])
  
  # Normalize gene class
  late_genes_n <- log((object$late_genes/object$nCount_RNA)*scale_factor + 1)
  object$late_genes_n <- late_genes_n
  early_genes_n <- log((object$early_genes/object$nCount_RNA)*scale_factor + 1)
  object$early_genes_n <- early_genes_n
  interm_genes_n <- log((object$interm_genes/object$nCount_RNA)*scale_factor + 1)
  object$interm_genes_n <- interm_genes_n
  
  # Save object_merged 
  save(object, list = c("object"), file = file.path(OUT_DIR, OUTFILE))
  return(object)
}

## --------- RUN ---------- ##
  
if (sys.nframe() == 0){
  cellCalling(MATRIX_DIR_RAW, MATRIX_DIR_FILT, OUT_DIR, MAX_UMI_MODEL,
               FDR, MIN_UMI_COUNT, scale_factor, MVA_genes_file, OUTFILE)
}
