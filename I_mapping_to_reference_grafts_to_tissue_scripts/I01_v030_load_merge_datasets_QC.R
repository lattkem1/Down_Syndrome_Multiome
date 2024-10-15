message("\n\n##########################################################################\n",
        "# Start I01: Load dataset into Seurat ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Plot cellranger QC, load count data, merge datasets, plot Seurat QC, \n",
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
library(patchwork)

#specify script/output index as prefix for file names
script_ind = "I01_"

#specify output directory
out_dir = paste0(main_dir,"I_mapping_to_reference_grafts_to_tissue/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#load group and file info for data to map on reference
gr_tab = read_csv(paste0(main_dir, "A_input/group_tab_grafts.csv"))



#######################################
#load cellranger QC data
#######################################

message("\n\n     *** Loading/plotting Cellranger QC data (", Sys.time(), ") \n\n")

set.seed(1234)

#load cellranger QC data

cellranger_QC = NULL

for (i in 1:nrow(gr_tab)){
  
  t1 = read_csv(paste0(gr_tab$dataset_folder[i],gr_tab$sample[i],"/outs/summary.csv"))
  cellranger_QC = rbind(cellranger_QC, t1)
  
}

write_csv(cellranger_QC, paste0(out_dir,script_ind,"cellranger_QC.csv"))


### plot cellranger QC

t1 = cellranger_QC

plot_cols = c('Estimated number of cells', 'GEX Mean raw reads per cell', 
              'GEX Median UMI counts per cell','GEX Median genes per cell',
              'GEX Fraction of transcriptomic reads in cells',
              'GEX Reads mapped confidently to genome', 
              'GEX Reads mapped confidently to transcriptome',
              'ATAC Mean raw read pairs per cell', 
              'ATAC Median high-quality fragments per cell',
              'ATAC Fraction of high-quality fragments in cells',
              'ATAC Fraction of genome in peaks', 
              'ATAC Fraction of high-quality fragments overlapping peaks',
              'ATAC Fraction of high-quality fragments overlapping TSS',
              'ATAC Number of peaks')

pl =  lapply(plot_cols, function(plot_col){
  t2 = t1[,c('Sample ID', plot_col)]
  names(t2) = c("sample", "pl_col")
  g1 = ggplot(t2, aes(x = sample, y = pl_col))+geom_col()+
    theme_classic()+ labs(title = plot_col, y = plot_col)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
})           

pdf(file = paste0(out_dir,script_ind,"cellranger_QC_barplots.pdf"), 
    width = nrow(gr_tab)/7+2, height = 40)
wrap_plots(pl, ncol=1)
dev.off()




#######################################
#load RNA, ATAC data, create Seurat objects and load cellranger QC data
#######################################

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

#load RNA and ATAC data
count_list = list()
seur_list = list()
fragpath_list = list()

for (i in 1:nrow(gr_tab)){
  
  message("\n\n     ***Loading count data for ", gr_tab$sample[i], 
          " (sample ",i, " of ", nrow(gr_tab)," )\n\n")
  
  counts = 
    Read10X_h5(paste0(gr_tab$dataset_folder[i],gr_tab$sample[i],"/outs/filtered_feature_bc_matrix.h5"))
  
  # create a Seurat object containing the RNA adata
  seur_temp <- CreateSeuratObject(
    counts = CreateAssayObject(counts = counts$`Gene Expression`),
    assay = "RNA"
  )
  
  #add sample metadata
  for (meta1 in colnames(gr_tab)){
    seur_temp[[meta1]] = gr_tab[[meta1]][i]
  }
  
  #add percent.mt QC
  v1 = seur_temp@assays$RNA@counts@Dimnames[[1]]
  v2 = v1[grepl(v1, pattern = "^MT-")]
  seur_temp[["percent.mt"]] <- PercentageFeatureSet(seur_temp, features = v2, assay = "RNA")
  
  
  ### add ATAC data
  
  fragpath = paste0(gr_tab$dataset_folder[i],
                    gr_tab$sample[i],"/outs/atac_fragments.tsv.gz")
  
  # create ATAC assay and add it to the object
  seur_temp[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  
  #add ATAC QC
  DefaultAssay(seur_temp) <- "ATAC"
  seur_temp <- NucleosomeSignal(seur_temp)
  seur_temp <- TSSEnrichment(seur_temp)
  
  #rename cells to make names unique to allow merging
  seur_temp <- RenameCells(object = seur_temp, add.cell.id = gr_tab$sample[i])
  
  seur_list[[ gr_tab$sample[i] ]] = seur_temp
  
}



#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n\n     *** Count data loaded. Saving Seurat objects... (", Sys.time(), ") \n\n")

save(seur_list, 
     file = paste0(out_dir,script_ind,"seur_list.rda")) 

gc(verbose = TRUE, reset = FALSE, full = TRUE)



#########################################
# Merge datasets 
#########################################

message("\n\n     *** Seurat objects saved. Merging Seurat objects... (", Sys.time(), ") \n\n")

#merge datasets preliminarily for QC plots

seur <- merge(seur_list[[1]], y = seur_list[-1], merge.data = TRUE)

rm(seur_list) #free up space
gc()

message("\n\n     *** Seurat objects merged. Saving combined dataset... (", Sys.time(), ") \n\n")

save(seur, file = paste0(out_dir,script_ind,"seur_merged.rda")) 



#########################################
# QC plots for merged dataset 
#########################################

message("\n\n     *** Combined dataset saved. Plotting Seurat QC... (", Sys.time(), ") \n\n")


###QC plots

#for large datasets, use subsample for plotting (else plots get  too large)
if (nrow(seur@meta.data)>50000){
  seur = seur[, sample(colnames(seur), size =50000, replace=F)]
} 

meta = seur@meta.data

pl = list()

pl[["nCount_RNA_max50k"]] = VlnPlot(seur, features = "nCount_RNA", group.by = "sample", pt.size = 0.01) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend()+ ylim(0, 50000)+ labs(title = "nCount_RNA_max50k") 
pl[["nCount_RNA_max3000"]] = VlnPlot(seur, features = "nCount_RNA", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  +
  NoLegend()+ ylim(0, 3000)+labs(title = "nCount_RNA_max3k")
pl[["nFeature_RNA"]] = VlnPlot(seur, features = "nFeature_RNA", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + 
  NoLegend()+ labs(title = "nFeature_RNA")
pl[["percent.mt"]] = VlnPlot(seur, features = "percent.mt", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + 
  NoLegend()+ labs(title = "percent.mt")
pl[["percent.mt_max5"]] = VlnPlot(seur, features = "percent.mt", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  NoLegend()+ylim(0,5)+ labs(title = "percent.mt_max5%")
pl[["nCount_ATAC_max100k"]] = VlnPlot(seur, features = "nCount_ATAC", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend()+ ylim(0, 100000)+ labs(title = "nCount_ATAC_max100k")
pl[["nCount_ATAC_max3000"]] = VlnPlot(seur, features = "nCount_ATAC", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend()+ ylim(0, 3000)+labs(title = "nCount_ATAC_max3k")
pl[["TSS.enrichment_max20"]] = VlnPlot(seur, features = "TSS.enrichment", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend()+ ylim(0, 20)+labs(title = "TSS.enrichment_max20")
pl[["nucleosome_signal"]] = VlnPlot(seur, features = "nucleosome_signal", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend()+labs(title = "nucleosome_signal")


pdf(file = paste0(out_dir,script_ind,"QC_plots.pdf"), width = 5, height = 6)
lapply(pl, function(x){x} )
dev.off()


#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed I01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


