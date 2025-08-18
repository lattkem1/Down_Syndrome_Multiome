message("\n\n##########################################################################\n",
        "# Start C01: Subsetting and re-integration ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E14_241219_DS_foetal_brain_grafts_for_man_v02_low_string/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#specify script/output index as prefix for file names
script_ind = "C01_"

#specify output directory
out_dir = paste0(main_dir,"C_subsetting_exc_lin_from_all_non_cx_excl/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}


#load group and file info; keep only samples
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered_non_cx_excl.csv")

#load cleaned/integrated dataset
load(file = paste0(main_dir,"C_subsetting_all_cells_non_cx_excl/C03_seur_integr_labelled.rda")) 

#define subset, load cluster labels for subsetting
clust_tab = read_csv(paste0(main_dir,"C_subsetting_all_cells_non_cx_excl/C02_subset_cluster_assignment_exc_lin.csv"))




###########################################################
# subset and re-integrate dataset
###########################################################

options(future.globals.maxSize = 30 * 1024^3) #SCTransform exceeds default memory limit for parallelisation with futures

message("\n\n          *** Subsetting and re-normalising dataset... ", Sys.time(),"\n\n")

seur <- subset(
  x = seur,
  subset = cluster_name %in% clust_tab$cluster_name & 
    sample %in% gr_tab$sample
  )

seur[["RNA"]] <- split(seur[["RNA"]], f = seur$sample)

#normalise full dataset (SCT transform)
gc()
seur <- SCTransform(seur, ncells = 5000,variable.features.n = 3000, 
                    assay = "RNA",new.assay.name = "SCT", 
                    conserve.memory = TRUE,
                    verbose = TRUE)

message("\n\n          *** Saving subsetted dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir,script_ind, "seur_subset.rda")) 


# re-integrate dataset

message("\n\n          *** Re-integrating subsetted dataset... ", Sys.time(),"\n\n")

DefaultAssay(seur) <- "SCT"

seur <- RunPCA(seur,  assay = "SCT", reduction.name = "pca")

seur <- IntegrateLayers(
  object = seur, method = HarmonyIntegration, assay = "SCT",
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

message("\n\n          *** Saving re-integrated dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir,script_ind, "seur_subset.rda")) 


#############################################
#joining RNA layers, re-normalising and dimension reduction
#############################################

message("\n\n          *** Re-normalising re-integrated dataset... ", Sys.time(),"\n\n")

seur <- JoinLayers(seur, assay = "RNA")

seur <- SCTransform(seur, ncells = 3000,variable.features.n = 2000, 
                    assay = "RNA",new.assay.name = "SCT", 
                    conserve.memory = TRUE,
                    verbose = TRUE)

seur <- RunPCA(seur,  assay = "SCT", reduction.name = "pca")
seur <- RunUMAP(seur, reduction = "harmony", dims = 1:30, return.model = TRUE)
seur = FindNeighbors(seur, reduction = "harmony", dims = 1:30)

message("\n\n          *** Saving Re-normalised re-integrated dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir,script_ind, "seur_subset.rda")) 




message("\n\n##########################################################################\n",
        "# Completed C01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


