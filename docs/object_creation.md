# Annotating Seurat object

## Rationale
The goal is to create a Seurat object from CellRanger count tables, and add to a Seurat object annotations relating to the nature of the genes (viral or cellular), as well as applying a different normalization. The normalization in question is the same as the one from the _Normalize_ function of the Seurat package, with the denominator being different : for viral/cellular genes, it is the total number of UMIs for viral/cellular genes.

## Steps

There are two parts to this, first is the use of the script ``cellcalling.R``. For each sample, it recover additional infected cells through EmptyDrops, apply additional filtering based on UMI threshold, and add some metadata to the Seurat object. It has some command line parameters :
1. Path to the directory with the raw count table
2. Path to the directory with the filtered count table
3. Directory output
4. Lower bound of total UMI count for EmptyDrops : cells below this UMI threshold will be discarded by EmptyDrops
5. FDR threshold used to call cells for EmptyDrops
6. Minimum UMI count in a cell to be retained
7. scale factor for the _Normalize_ data
8. semi-colon separated file giving the "common" gene name equivalent to the official MVA gene names
9. tsv separated file indicating the "class" of the MVA genes (early, intermediary or late)
10. Name of the output R object file

The script is followed by ``create_Seurat_object.R``, to allow the normalization of feature depending of their nature, and to add additional metadata. It also merge samples together to create infected and uninfected Seurat objects. It has some command line parameters :
1. Path of the Seurat object of the uninfected replicate 1
2. Path of the Seurat object of the uninfected replicate 2
3. Path of the Seurat object of the uninfected replicate 1
4. Path of the Seurat object of the uninfected replicate 2
5. semi-colon separated file giving the "common" gene name equivalent to the official MVA gene names
6. tsv separated file indicating the "class" of the MVA genes (early, intermediary or late)
7. boolean, wether to normalize viral and cellular counts separately or not
8. Output path where the objects will be created

This will create in path_out three R objects : one with the infected data, one with the uninfected data, and one with all the replicates.
