message("\n\n##########################################################################\n",
        "# Start F03: scMEGA network generation ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
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
library(igraph)
library(ggraph)
library(tidygraph)

#specify script/output index as prefix for file names
script_ind = "F03_"

#specify output directory
out_dir = paste0(main_dir,"F_Chromatin_scMEGA_GRN_analysis_exc_lin_from_all_non_cx_excl_PCW16_20_v045/")

#load group and file info
# gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")

#load Seurat dataset with ArchR analysis

# clust_tab = read_csv(paste0(main_dir,"B_basic_analysis/B03_cluster_assignment_exc_lin.csv"))

load( file = paste0(out_dir, "F02_seur_with_ArchR.rda"))


#load pseudobulk dataset

load(file = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_PCW16_20_v045/",
                   "E03_bulk_data_w_expr_z_scores.rda"))



#########################################
# select TFs and genes for GRN inference
#########################################

### select network TFs and plot activity


message("\n\n          *** Select network TFs... ", Sys.time(),"\n\n")
# 
# sel.tfs <- SelectTFs(object = seur, 
#                      tf.assay = "chromvar",
#                      rna.assay = "SCT",
#                      atac.assay = "peaks_by_cluster",
#                      trajectory.name = "Trajectory",
#                      groupEvery = 1,
#                      return.heatmap = TRUE,
#                      cor.cutoff = -1,
#                      p.cutoff = 1)

# SelectTFs error when expression of TF is 0 throughout trajectory => 
# alternative: manually get all TFs with expression and in chromvar assay

trajMM <- GetTrajectory(seur, assay = "chromvar", trajectory.name = "Trajectory", 
                        groupEvery = 1, slot = "data", smoothWindow = 7, 
                        log2Norm = FALSE)
t2 = trajMM@assays@data$mat
rownames(trajMM) <- seur@assays[["peaks_by_cluster"]]@motifs@motif.names
trajRNA <- GetTrajectory(seur, assay = "SCT", trajectory.name = "Trajectory", 
                         groupEvery = 1, slot = "data", smoothWindow = 7, 
                         log2Norm = TRUE)

mat1 <- assay(trajMM)
mat1 = mat1[apply(abs(mat1), 1, sum)>0,]
mat2 <- assay(trajRNA)
mat2 = mat2[apply(abs(mat2), 1, sum)>0,]

TFs.use <- intersect(rownames(mat1), rownames(mat2))


### select network genes and plot activity

message("\n\n          *** Select network genes... ", Sys.time(),"\n\n")

sel.genes <- SelectGenes(object = seur,
                          rna.assay = "SCT",
                          atac.assay = "peaks_by_cluster",
                         labelTop1 = 0,
                         labelTop2 = 0)

df.p2g <- sel.genes$p2g
ht <- sel.genes$heatmap

pdf(file = paste0(out_dir,script_ind, "Network_genes_activity_heatmap.pdf"), 
    width = 6, height = 6)
{
  draw(ht)
}

dev.off()



#############################################################################
# infer network (calculate TF-gene correlation)
#############################################################################

message("\n\n          *** Calculate GRN... ", Sys.time(),"\n\n")

tf.gene.cor <- GetTFGeneCorrelation(object = seur, 
                                    tf.use = TFs.use, 
                                    gene.use = unique(intersect(df.p2g$gene, 
                                                                c(unlist(bulk_data$DEGs), TFs.use))),
                                    tf.assay = "chromvar", 
                                    gene.assay = "SCT",
                                    atac.assay = "peaks_by_cluster",
                                    trajectory.name = "Trajectory")


save(seur, tf.gene.cor, sel.genes, TFs.use, file = paste0(out_dir, script_ind, "seur_with_GRN.rda"))



#############################################################################
# extract network
#############################################################################

message("\n\n          *** Extract GRN data... ", Sys.time(),"\n\n")


###matching motifs with regulatory elements to identify direct targets

df.p2g <- sel.genes$p2g

motif.matching <- seur@assays$peaks_by_cluster@motifs@data
colnames(motif.matching) <- seur@assays$peaks_by_cluster@motifs@motif.names
motif.matching <-
  motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]


###extract network

df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = tf.gene.cor, 
                 df.p2g = df.p2g)



message("\n\n          *** Saving GRN dataset... ", Sys.time(),"\n\n")

save(seur, tf.gene.cor, sel.genes, df.p2g, df.grn, motif.matching,  TFs.use,
     file = paste0(out_dir, script_ind, "seur_with_GRN.rda"))



#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed F03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


