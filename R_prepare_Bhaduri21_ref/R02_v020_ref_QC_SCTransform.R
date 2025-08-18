### load normalise reference atlas dataset
# step 1: generate group table, split and downsample dataset, SCT normalise

message("###################\n####### Start R02: ", Sys.time(),"\n########################")

setwd("/rds/general/user/mlattke/projects/dsseq23/live/E02_230608_foetal_brain_ref_atlas/")
set.seed(1234)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(colorRamps)

# define output folder
out_dir = "./R_out_Seur5/"

### load merged Seurat reference dataset
load (file = paste0(out_dir,"R01_seur_ref.rda"))

#create gr_tab and sample/region overview from Seurat metadata
meta = seur_ref@meta.data
names(meta)
t1 = tibble(meta[,c("orig.ident", "individual", "age", "structure", "area", "subarea", "lamina")])
t2 = distinct(t1)
t2$group = paste0("GW", t2$age)
t2$sample = t2$individual
write_csv(t2, file = paste0(out_dir, "02_sample_region_overview.csv"))
gr_tab_ref = distinct(t2[,c("group", "sample", "individual", "age")])
write_csv(gr_tab_ref, file = "R_in/gr_tab_ref.csv")

#add group/sample information and umap from UCSC cell browser
seur_ref$sample = seur_ref$individual
seur_ref$group = paste0("GW", seur_ref$age)
seur_ref$dataset = "Ref_Bhadhuri21"
# calculate mitochondrial content 
seur_ref[["percentMito"]] = PercentageFeatureSet(seur_ref, pattern = "^MT-", assay = "RNA")


#get marker gene panels
marker_tab = read_csv("R_in/cell type markers 22-02-23.csv")

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###data loaded  ", Sys.time(),"\n")

####################################################################
# remove low quality cells, split reference dataset into list with samples, downsample to 500 cells/cluster and sample
####################################################################

#remove low quality cells
seur_ref <- subset(
  x = seur_ref,
  subset = 
    percentMito <= 10 &
    nFeature_RNA >= 750 
)

Idents(seur_ref) = "cell_cluster" #required for downsampling

seur_list_merged = SplitObject(seur_ref, split.by = "sample")

#free up space
rm(seur_ref)
gc()


###########################################################
# SCTransform
###########################################################

seur_list_merged <- lapply(X = seur_list_merged, FUN = SCTransform, ncells = 5000, variable.features.n = 3000)

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###data normalised  ", Sys.time(),"\n")


save(seur_list_merged, file = paste0(out_dir,"R02_ref_seur_integr.rda")) 

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n###data saved  ", Sys.time(),"\n")

message("/n####### End R02b: ", Sys.time(),"\n")

