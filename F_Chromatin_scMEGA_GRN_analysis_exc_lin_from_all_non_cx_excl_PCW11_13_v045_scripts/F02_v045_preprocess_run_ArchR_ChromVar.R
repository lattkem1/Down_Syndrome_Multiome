message("\n\n##########################################################################\n",
        "# Start H01: scMEGA setup ", Sys.time(),
        "\n##########################################################################\n",
        "\n   downsample Seurat, get TF motifs, run ArchR analysis",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E14_241219_DS_foetal_brain_grafts_for_man_v02_low_string/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(ArchR)
library(scMEGA)
library(GenomeInfoDb)
library(JASPAR2024)
library(TFBSTools)
library(ggraph)
library(BiocParallel)

#specify script/output index as prefix for file names
script_ind = "F02_"

#specify output directory
out_dir = paste0(main_dir,"F_Chromatin_scMEGA_GRN_analysis_exc_lin_from_all_non_cx_excl_PCW11_13_v045/")

#load group and file info
#gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")


#clust_tab = read_csv(paste0(main_dir,"B_basic_analysis/B03_cluster_assignment_exc_lin.csv"))

load(file = paste0(out_dir, "F01_seur_w_peaks_by_cluster_quant.rda")) 


#load RNA pseudobulk dataset (to extract cell populations for GRN analysis)

load(file = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_PCW11_13_v045/",
                   "E03_bulk_data_w_expr_z_scores.rda"))


###focus on main excitatory lineage (remove NEU_RELN and AST_OPC)
cluster_names = unique(bulk_data$meta$cluster_name)
cluster_names = cluster_names[!cluster_names %in% c("AST_s9")]


##############################################################
# save peaks with annotations
##############################################################

annot = as_tibble(Annotation(seur[["peaks_by_cluster"]]))

# identify nearest TSS

DefaultAssay(seur) = "peaks_by_cluster"

t1 = as_tibble(distanceToNearest(seur, subject = Annotation(seur)))

t2 = as_tibble(seur@assays$peaks_by_cluster@ranges)
t2$peak_id = paste0(t2$seqnames,"-",t2$start, "-", t2$end)
t2$distanceToNearest = t1$distance
t2$gene_name = annot$gene_name[t1$subjectHits]
peak_annot = t2

write_csv(peak_annot, file = paste0(out_dir, script_ind, "peaks_annotated.csv"))


#############################################################################
# subset dataset, retain only selected clusters
#     use only clusters also used for DESeq2 analysis (at least 2 pseudobulks per cluster_stage)
#############################################################################

t1 = seur@meta.data
t1$cluster_sample = paste0(t1$cluster_name,"_",t1$sample)

t2 = t1 %>% group_by(cluster_sample, cluster_name) %>%  dplyr::summarise(N_cells = n())

t2 = t2[t2$cluster_sample %in% bulk_data$meta$cluster_sample & t2$cluster_name %in% cluster_names,]

t3 = t1[t1$cluster_sample %in% t2$cluster_sample,]

# v3 = NULL
# 
# for (cl in t2$cluster_sample){
#   v1 = rownames(t1)[t1$cluster_sample==cl]
#   if (length(v1) >100){
#     set.seed(123)
#     v2 = sample(v1, size = 100)
#   } else {v2 = v1}
#   v3 = c(v3,v2)
# }

seur = subset(seur, cells = rownames(t3))

gc()


#############################################################################
# infer trajectory
#############################################################################
cluster_names = unique(bulk_data$meta$cluster_name)

seur <- AddTrajectory(object = seur, 
                      trajectory = cluster_names, #include all cells
                      group.by = "cluster_name", 
                      reduction = "umap",
                      dims = 1:2, 
                      use.all = TRUE)

# we only plot the cells that are in this trajectory
seur <- seur[, !is.na(seur$Trajectory)]


#visualise trajectory
p1 <- DimPlot(object = seur, 
              group.by = "cluster_name", 
              reduction = "umap",
              label = TRUE) + NoLegend()

p2 <- TrajectoryPlot(object = seur, 
                     reduction = "umap",
                     continuousSet = "blueYellow",
                     size = 1,
                     addArrow = TRUE)



pdf(file = paste0(out_dir,script_ind, "Trajectory.pdf"), 
    width = 10, height = 5)
{
  
  p1 + p2
  
}

dev.off()




#############################################################################
# initiate network, get TF motifs, run chromVar (for TF activity)
#############################################################################

message("\n\n          *** Get motif database... ", Sys.time(),"\n\n")

JASPAR2024 <- JASPAR2024()

# Get a list of motif position frequency matrices from the JASPAR database
pfm = getMatrixSet(x = JASPAR2024@db, 
                    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
                    )

message("\n\n          *** Add motifs to peak set... ", Sys.time(),"\n\n")


# add motif information

seur <- AddMotifs(
  object = seur,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "peaks_by_cluster"
)

message("\n\n          *** Run ChromVar... ", Sys.time(),"\n\n")

gc()


#limit number of cores for ChromVar, else extremely high memory requirement

register(MulticoreParam(4))

# run chromVAR

seur <- RunChromVAR(
  object = seur,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "peaks_by_cluster"
)

gc()

message("\n\n          *** Save dataset... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir, script_ind, "seur_with_ArchR.rda"))

#load(file = paste0(out_dir, script_ind, "seur_with_ArchR.rda"))

#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed H01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


