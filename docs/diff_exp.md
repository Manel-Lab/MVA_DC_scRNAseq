# Clustering

## Rationale
These scripts will take a Seurat object with a clustering solutions aswell as a Seurat object with Uninfected samples, and do a differential expression between the clusters and the Uninfected samples.

## Steps

### 1) Cellular genes

The script is ``main_DE_cellgenes.R`` will do a differential expression analysis on the infected clusters against all Uninfected cells or other clusters, aswell as the pathways analysis of the down-regulated and up-regulated genes of these clusters separately. Command lines arguments are :

1. Seurat object of infected samples
2. Seurat object of uninfected samples
3. path to args.json file
4. output path

Some arguments are passed through the json argument file

- __logFC__ : The logFC threshold (absolute value)
- __pvalue__ : the pvalue threshold
- __min.pct__ : minimum percent of cells expressing the differential gene
- __test__ : which test from Seurat's _FinderMarkers_ to use
- __cluster_sol__ : name of the cluster solution to use
- __scaled__ : wether to use scaled data ("TRUE") or not ("FALSE") for heatmaps
- __uninf__ : if TRUE, compare each cluster to the uninfected cells. Otherwise, compare each cluster to the other ones combined.

In the output path, there will be :
- heatmaps representing the top 50 hits in each cluster
- a heatmap with the union of the top 50
- a recapitulating tsv with all differential genes
- a heatmap with the logFC of top hits across cluster
- DE genes table for each clusters
- pathway analysis (tables, plots) for each down and up regulated genes in each cluster.
