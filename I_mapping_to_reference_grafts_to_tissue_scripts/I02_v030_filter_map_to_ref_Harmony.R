message("\n\n##########################################################################\n",
        "# Start I02: Filter dataset and map to reference ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Remove low quality cells, normalise dataset, map to reference \n",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#specify script/output index as prefix for file names
script_ind = "I02_"

#specify output directory
out_dir = paste0(main_dir,"I_mapping_to_reference_grafts_to_tissue/")

#load reference dataset
load(file = paste0(main_dir,"B_basic_analysis/B04_seur_integr_labelled.rda") )
seur_ref = seur

gr_tab_ref = read_csv(file =paste0(main_dir, "B_basic_analysis/B02_gr_tab_filtered.csv"))


#load data to map on reference

gr_tab_data = read_csv(paste0(main_dir, "A_input/group_tab_grafts.csv"))

load(file = paste0(out_dir,"I01_seur_merged.rda") )



#############################################
#filter low quality cells
#############################################

message("\n\n          *** Removing low quality cells... ", Sys.time(),"\n\n")

set.seed(1234)
DefaultAssay(seur) <- "RNA"

filter_stats = seur@meta.data %>% group_by(sample) %>% summarise(N_cells_unfiltered = n())

set.seed(123)

seur = subset(x = seur,
              subset = nCount_RNA < 30000 &
                nCount_RNA > 500 &
                percent.mt < 2
)

t1 = seur@meta.data %>% group_by(sample) %>% summarise(N_cells_filtered = n())

filter_stats$N_cells_filtered = t1$N_cells_filtered[match(filter_stats$sample, t1$sample)]

filter_stats$fract_removed = (filter_stats$N_cells_unfiltered - filter_stats$N_cells_filtered)/filter_stats$N_cells_unfiltered

write_csv(filter_stats, file = paste0(out_dir,script_ind, "filter_stats.csv"))



###########################################################
# normalise dataset
###########################################################

message("\n\n          *** Normalising dataset... ", Sys.time(),"\n\n")

seur <- SCTransform(seur, ncells = 3000,variable.features.n = 2000, 
                    assay = "RNA",new.assay.name = "SCT", 
                    conserve.memory = TRUE,
                    verbose = TRUE)

message("\n\n          *** Normalisation finished. Saving dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir,script_ind,"seur_filtered_mapped.rda")) 



###########################################################
# map own dataset (seur_list) on reference
###########################################################

message("\n\n          *** Dataset saved. Mapping to reference... ", Sys.time(),"\n\n")

anchors <- FindTransferAnchors(reference = seur_ref, query = seur,
                               dims = 1:30, reference.reduction = "pca")

seur <- MapQuery(anchorset = anchors, reference = seur_ref, query = seur,
                 refdata = list(cluster_name = "cluster_name", cell_type = "cell_type",
                                cell_class = "cell_class",
                                stage = "stage", dev_PCW = "dev_PCW",
                                group = "group"), 
                 reference.reduction = "pca", reduction.model = "umap")

message("\n\n          ### Mapping completed. Saving dataset... ", Sys.time(),"\n\n")


save(seur, file = paste0(out_dir,script_ind,"seur_filtered_mapped.rda"))




#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed I02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


