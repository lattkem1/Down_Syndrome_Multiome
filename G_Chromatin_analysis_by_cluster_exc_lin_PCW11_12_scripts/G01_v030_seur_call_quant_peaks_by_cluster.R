message("\n\n##########################################################################\n",
        "# Start G01: Call peaks by cluster  ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
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
script_ind = "G01_"

#specify output directory
out_dir = paste0(main_dir,"G_Chromatin_analysis_by_cluster_exc_lin_PCW11_12/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#load dataset

load(file = paste0(main_dir,"C_subsetting_exc_lin_PCW11_12/C03_seur_integr_labelled.rda")) 




###########################################################
# call peaks by cluster_name
###########################################################

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

DefaultAssay(seur) = "ATAC"

peaks <- CallPeaks(
  object = seur,
  group.by = "cluster_name"
)


# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks  <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

#save peakset as bed
export(peaks, paste0(out_dir, script_ind, "peaks_by_cluster.bed"))

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seur),
  features = peaks,
  cells = colnames(seur)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
seur[["peaks_by_cluster"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(seur),
  annotation = annotation
)


###access peaks normalisation

DefaultAssay(seur) <- "peaks_by_cluster"
seur <- FindTopFeatures(seur, min.cutoff = 5)
seur <- RunTFIDF(seur)
seur <- RunSVD(seur)


save(seur, peaks, file = paste0(out_dir, script_ind,"seur_w_peaks_by_cluster_quant.rda"))



#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed G01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


