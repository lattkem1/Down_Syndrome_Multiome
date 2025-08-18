message("\n\n##########################################################################\n",
        "# Start B02: Filter normalise and integrate dataset ", Sys.time(),
        "\n##########################################################################\n",
        "\n   \n", 
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
library(colorRamps)
library(GenomicRanges)

#specify script/output index as prefix for file names
script_ind = "B02_"

#specify output directory
out_dir = paste0(main_dir,"B_basic_analysis/")

#load group and file info
gr_tab = read_csv("A_input/group_tab_tissue.csv")

#load raw dataset
load(file = paste0(out_dir,"B01_seur_merged.rda")) 


#############################################
#filter low quality cells
#############################################


message("\n\n          *** Removing low quality cells... ", Sys.time(),"\n\n")

set.seed(1234)
DefaultAssay(seur) <- "RNA"

filter_stats = seur@meta.data %>% group_by(library) %>% summarise(N_cells_unfiltered = n())

set.seed(123)

seur = subset(x = seur,
              nCount_ATAC < 25000 &
                nCount_RNA < 30000 &
                nCount_ATAC > 100 &
                nCount_RNA > 500 &
                percent.mt < 2 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1.1
)


t1 = seur@meta.data %>% group_by(library) %>% summarise(N_cells_filtered = n())

filter_stats$N_cells_filtered = t1$N_cells_filtered[match(filter_stats$library, t1$library)]

filter_stats$fract_removed = (filter_stats$N_cells_unfiltered - filter_stats$N_cells_filtered)/filter_stats$N_cells_unfiltered

# remove libraries with >50% cells removed due to low quality

libraries_retain = filter_stats$library[filter_stats$fract_removed<=0.5 & 
                                       filter_stats$N_cells_filtered>=500 &
                                       filter_stats$library %in% gr_tab$library]
filter_stats$cells_retained = 0
filter_stats$cells_retained[filter_stats$library %in% libraries_retain] = 
  filter_stats$N_cells_filtered[filter_stats$library %in% libraries_retain]

write_csv(filter_stats, file = paste0(out_dir,script_ind, "filter_stats.csv"))

seur = subset(x = seur, subset = library %in% libraries_retain)

gr_tab = gr_tab[gr_tab$library %in% libraries_retain,]
write_csv(gr_tab, file = paste0(out_dir,script_ind,"gr_tab_filtered.csv"))

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)


###########################################################
# integrate data with Harmony
###########################################################

message("\n\n          *** Normalising dataset... ", Sys.time(),"\n\n")

options(future.globals.maxSize = 30 * 1024^3) #SCTransform exceeds default memory limit for parallelisation with futures

seur[["RNA"]] <- split(seur[["RNA"]], f = seur$library)

set.seed(123)

seur <- SCTransform(seur, ncells = 3000,variable.features.n = 2000, 
                        assay = "RNA",new.assay.name = "SCT", 
                        conserve.memory = TRUE,
                        verbose = TRUE)

message("\n\n          *** Dataset normalised. Save dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir,script_ind, "seur_integr.rda")) 

message("\n\n          *** Integrate dataset... ", Sys.time(),"\n\n")

set.seed(123)

seur <- RunPCA(seur,  assay = "SCT", reduction.name = "pca")

set.seed(123)

seur <- IntegrateLayers(
  object = seur, method = HarmonyIntegration, assay = "SCT",
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

message("\n\n          *** Integration done. Saving dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir,script_ind, "seur_integr.rda")) 

#############################################
#joining RNA layers, re-normalising and initial clustering
#############################################

message("\n\n          *** Re-normalising integrated dataset... ", Sys.time(),"\n\n")

seur <- JoinLayers(seur, assay = "RNA")

set.seed(123)

seur <- SCTransform(seur, ncells = 3000,variable.features.n = 2000, 
                    assay = "RNA",new.assay.name = "SCT", 
                    conserve.memory = TRUE,
                    verbose = TRUE)

set.seed(123)
seur <- RunPCA(seur,  assay = "SCT", reduction.name = "pca")

set.seed(123)
seur <- RunUMAP(seur, reduction = "harmony", dims = 1:30, return.model = TRUE)
set.seed(123)
seur = FindNeighbors(seur, reduction = "harmony", dims = 1:30)

message("\n\n          *** Dataset normalised. Save dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir,script_ind, "seur_integr.rda")) 



message("\n\n##########################################################################\n",
        "# Completed B02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


