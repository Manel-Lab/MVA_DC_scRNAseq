## ---------------------------
##
## Purpose of script: Pipeline to recreate our analysis
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-02-20
##
## ---------------------------

## Arguments
.libPaths(c("~/R_packages",.libPaths()))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  print("Not enough arguments provided")
  print("Usage : Rscript inf1 inf2 uninf1 uninf2 mva_tsv class_tsv normalize path_out jfile")
  q()
} else if (length(args) > 7) {
  print("Too many arguments provided")
  print("Usage : Rscript inf1 inf2 uninf1 uninf2 mva_tsv class_tsv normalize path_out jfile")
  q()
}


indir <- args[1]        # Input directory with cellranger alignment results
mva_map <- args[2]      # csv file mapping VACWR terminology of genes to gene names
class_tsv <- args[3]    # Path to tsv annotating MVA genes
normalize <- args[4]    # Wether or not to use the alt_norm function
path_out <- args[5]     # Path of the outputed R object
jfile <- args[6]        # Json argument file
working_dir <- args[7]  # Script directory

## Global variables

# Read from json files
result <- rjson::fromJSON(file = jfile)

# CD to worind_dir
setwd(working_dir)

## ---------------------------

## 

source("main_cellcalling.R")
source("main_create_object.R")
source("main_integration.R")
source("main_cluster_red_solution.R")
source("main_clustering.R")
source("main_clust_eval.R")
source("main_DE_cellgenes.R")
source("paper_figures.R")

## --------- RUN ---------- ##

print("1) Create objects from cellranger data")
inf1 <- cellCalling(file.path(indir, "MVA-rep1", "outs", "raw_feature_bc_matrix"), 
                    file.path(indir, "MVA-rep1", "outs", "filtered_feature_bc_matrix"), 
                    path_out, 100, 0.1, 1000, 10000, mva_map, class_tsv, "MVA-rep1.Robj")
inf2 <- cellCalling(file.path(indir, "MVA-rep2", "outs", "raw_feature_bc_matrix"), 
                    file.path(indir, "MVA-rep2", "outs", "filtered_feature_bc_matrix"), 
                    path_out, 100, 0.1, 1000, 10000, mva_map, class_tsv, "MVA-rep2.Robj")
uninf1 <- cellCalling(file.path(indir, "Uninfected-rep1", "outs", "raw_feature_bc_matrix"), 
                      file.path(indir, "Uninfected-rep1", "outs", "filtered_feature_bc_matrix"), 
                      path_out, 100, 0.1, 1000, 10000, mva_map, class_tsv, "Uninf-rep1.Robj")
uninf2 <- cellCalling(file.path(indir, "Uninfected-rep2", "outs", "raw_feature_bc_matrix"), 
                      file.path(indir, "Uninfected-rep2", "outs", "filtered_feature_bc_matrix"), 
                      path_out, 100, 0.1, 1000, 10000, mva_map, class_tsv, "Uninf-rep2.Robj")
#inf1 <- cellCalling2(file.path(indir, "MVA-rep1", "filtered_feature_bc_matrix"), 
#                    path_out, 10000, 1000, 100, mva_map, class_tsv, "MVA-rep1.Robj")
#inf2 <- cellCalling2(file.path(indir, "MVA-rep2", "filtered_feature_bc_matrix"), 
#                     path_out, 10000, 1000, 100, mva_map, class_tsv, "MVA-rep2.Robj")
#uninf1 <- cellCalling2(file.path(indir, "Uninfected-rep1", "filtered_feature_bc_matrix"), 
#                       path_out, 10000, 1000, 100, mva_map, class_tsv, "Uninfected-rep1.Robj")
#uninf2 <- cellCalling2(file.path(indir, "Uninfected-rep2", "filtered_feature_bc_matrix"), 
#                       path_out, 10000, 1000, 100, mva_map, class_tsv, "Uninfected-rep2.Robj")
print("2) Creating Seurat objects with relevant annotations and correct normalization")
main_create_object(inf1, inf2, uninf1, uninf2, mva_map, class_tsv, path_out)
print("3) Data reduction and removal of the donor effect")
integrate_all(result, file.path(path_out, "MVA_sepnorm.Robj"), 
              file.path(path_out, "Uninfected_sepnorm.Robj"), 
              file.path(path_out, "All-samples_sepnorm.Robj"),
              file.path(path_out, "donor_effect_removal"), 
              file.path(path_out, "donor_effect_removal"), 
              file.path(path_out, "donor_effect_removal"))
print("4) Evaluation of the reduction solution")
main_red_eval(file.path(path_out, "donor_effect_removal", "New_MNN.rds"), 
              file.path(path_out, "donor_effect_removal", "cluster_eval"),
              result)
print("5) Clustering of the reduction solution")
main_clustering(file.path(path_out, "donor_effect_removal", "New_MNN.rds"), 
                file.path(path_out, "donor_effect_removal", "clustering"),
                result)
print("6) Evaluating clustering solution")
main_cluster_eval(file.path(path_out, "donor_effect_removal", "clustering", "10k_Clusters_All-samples.rds"), 
                  result,
                  file.path(path_out, "donor_effect_removal", "clustering", "Cluster_eval"))
print("7) Cellular differential expression and pathway analysis")
main_DE_cell(file.path(path_out, "donor_effect_removal", "clustering", "10k_Clusters_All-samples.rds"), 
             file.path(path_out, "Uninfected_sepnorm.Robj"), 
             file.path(path_out, "donor_effect_removal", "clustering", "cellgenes_diffexpr"),
             result)
print("8) Various figures for the paper")
main_paper_figures(file.path(path_out, "donor_effect_removal", "clustering", "10k_Clusters_All-samples.rds"), 
                   file.path(path_out, "Uninfected_sepnorm.Robj"), 
                   file.path(path_out, "paper_figures"),
                   file.path(path_out, "donor_effect_removal", "clustering", "cellgenes_diffexpr"),
                   result)
