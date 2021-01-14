# Clustering

## Rationale
These scripts will create most of the figures of the paper.

## Steps

The script ``paper_figures.R`` will create _FeaturePlot_ with genes of interest (____Figure 1G.__), barplots with repartition of viral genes classes among viral UMIs (__Figure 3A.__ ), heatmap of top 25 genes in each cluster (__Figure S11B.__), barplots with the correlation (__Figure S1D__), aswell violin plots (__Figure 6B__) comparing level expressions of genes across our clusters. The command line arguments are :

1. Path of infected Seurat object
2. Path to uninfected Seurat object
3. Path to argument json file
4. Output path
5. Directory where the differential tables are

There are also arguments contained inside the json argument file :

- __red__ : key of the reduction to use for FeaturePlot
- __d__ : number of dimensions for FeaturePlot
- __list_genes__ : list of genes of interest for FeaturePlot
- __threshold__ : threshold of number of MVA UMIs for the barplot : below these threhsold, cells are discarded.
- __seurat_cluster__ : key of the cluster solution of interest
- __scaled__ : for the heatmap, use scaled value or not (boolean)
- __DE_files__ : list of differential files for the heatmap
- __top__ : top number of DE genes of each cluster to plot in the heatmap
- __list_HM__ : list of genes to plot in a heatmap, different from the one with the top 25 genes
- __list_vln__ : list of features for Violin Plot
- __list_cor__ : list of feature to plot for the correlation with the ratio of MVA.
