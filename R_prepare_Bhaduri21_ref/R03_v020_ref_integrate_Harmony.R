### load and merge reference and analysis datasets
# perform basic QC plots and preliminary clustering for initial assessment of dataset

message("###################\n####### Start R03b: ", Sys.time(),"\n########################")

setwd("/rds/general/user/mlattke/projects/dsseq23/live/E02_230608_foetal_brain_ref_atlas/")
set.seed(1234)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(colorRamps)

# define output folder
out_dir = "./R_out_Seur5/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

### load normalised Seurat reference dataset
load(file = paste0(out_dir,"R02_ref_seur_integr.rda")) 

gr_tab_ref = read_csv(file = "R_in/gr_tab_ref.csv")

#get marker gene panels
marker_tab = read_csv("R_in/cell type markers 22-02-23.csv")

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###data loaded  ", Sys.time(),"\n")


###########################################################
# integrate data with Harmony
###########################################################

seur_ref = merge(seur_list_merged[[1]],seur_list_merged[-1])

#checkpoint
rm(seur_list_merged)
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###data merged  ", Sys.time(),"\n")

seur_ref <- SCTransform(seur_ref, ncells = 3000,variable.features.n = 2000, 
                    assay = "RNA",new.assay.name = "SCT", 
                    conserve.memory = TRUE,
                    verbose = TRUE)

gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###merged data SCT transformed  ", Sys.time(),"\n")

seur_ref <- RunPCA(seur_ref,  assay = "SCT", reduction.name = "pca")

seur_ref <- IntegrateLayers(
  object = seur_ref, method = HarmonyIntegration, assay = "SCT",
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

save(seur_ref, file = paste0(out_dir,"R03_ref_seur_integr.rda")) 

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###Seurat integrated  ", Sys.time(),"\n")

#joining RNA layers, re-normalising and initial clustering
seur_ref <- JoinLayers(seur_ref, assay = "RNA")
seur_ref <- SCTransform(seur_ref, ncells = 3000,variable.features.n = 2000, 
                        assay = "RNA",new.assay.name = "SCT", 
                        conserve.memory = TRUE,
                        verbose = TRUE)

seur_ref <- RunPCA(seur_ref,  assay = "SCT", reduction.name = "pca")
seur_ref <- RunUMAP(seur_ref, reduction = "harmony", dims = 1:30, return.model = TRUE)
seur_ref = FindNeighbors(seur_ref, reduction = "harmony", dims = 1:30)

gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###save integrated dataset ", Sys.time(),"\n")

save(seur_ref, file = paste0(out_dir,"R03_ref_seur_integr.rda")) 

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###data saved  ", Sys.time(),"\n")


message("/n####### End R03b: ", Sys.time(),"\n")


