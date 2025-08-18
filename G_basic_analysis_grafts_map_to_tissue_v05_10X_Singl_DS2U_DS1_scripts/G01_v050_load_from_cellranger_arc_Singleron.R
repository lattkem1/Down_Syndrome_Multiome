message("\n\n##########################################################################\n",
        "# Start B01: Load dataset into Seurat (10X Multiome RNA and Singleron RNA)", Sys.time(),
        "\n##########################################################################\n",
        "\n   load count data, merge datasets, plot Seurat QC, \n",
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
library(patchwork)


#specify script/output index as prefix for file names
script_ind = "G01_"

#specify output directory
out_dir = paste0(main_dir,"G_basic_analysis_grafts_map_to_tissue_v05_10X_Singl_DS2U_DS1/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#load group and file info for analysis dataset
gr_tab = read_csv("A_input/group_tab_grafts_DS2U_DS1_SCZ.csv")





#######################################
#load RNA data, create Seurat objects 
#######################################

seur_list = list()


for (i in 1:nrow(gr_tab)){
  
  message("\n\n     ***Loading count data for ", gr_tab$library[i], 
          " (library ",i, " of ", nrow(gr_tab)," )\n\n")
  
  ### load RNA counts for 10X multiome data 

  counts_rna = NULL
  
  if (gr_tab$seq_tech[i] == "Multiome_10X"){
    
    l1 = Read10X_h5(paste0(gr_tab$dataset_folder[i],gr_tab$library[i],"/outs/filtered_feature_bc_matrix.h5"))
    counts_rna = l1$`Gene Expression`
    
  }
 
  ### load RNA counts for Singleron data 
  
  if (gr_tab$seq_tech[i] == "snRNA_Singleron"){
    
    #unzip tar folder in temporary folder
    
    untar(paste0(gr_tab$dataset_folder[i], gr_tab$library[i],"/",
                 gr_tab$library[i],"_counts_matrix.tar"), 
          exdir = paste0(out_dir, "temp/"))
    
    counts_rna = ReadMtx(mtx = paste0(out_dir, "temp/matrix.mtx.gz"),
                     cells = paste0(out_dir, "temp/barcodes.tsv.gz"),
                     features = paste0(out_dir, "temp/features.tsv.gz"))
    
    #delete unzipped temporary folder
    unlink(paste0(out_dir, "temp"), recursive = TRUE) 
    
  }
  
  
  ### create a Seurat object containing the RNA counts
  
  seur_temp <- CreateSeuratObject(
    counts = CreateAssayObject(counts = counts_rna),
    assay = "RNA"
  )
  
  #add library metadata
  for (meta1 in colnames(gr_tab)){
    seur_temp[[meta1]] = gr_tab[[meta1]][i]
  }

  #add percent.mt QC
  v1 = seur_temp@assays$RNA@counts@Dimnames[[1]]
  v2 = v1[grepl(v1, pattern = "^MT-")]
  seur_temp[["percent.mt"]] <- PercentageFeatureSet(seur_temp, features = v2, assay = "RNA")
  
  #rename cells to make names unique to allow merging
  seur_temp <- RenameCells(object = seur_temp, add.cell.id = gr_tab$library[i])
  
  seur_list[[ gr_tab$library[i] ]] = seur_temp

}


gc(verbose = TRUE, reset = FALSE, full = TRUE)



#########################################
# Merge datasets 
#########################################

message("\n\n     *** Merging Seurat objects... (", Sys.time(), ") \n\n")

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

pl[["nCount_RNA_max50k"]] = VlnPlot(seur, features = "nCount_RNA", group.by = "library", pt.size = 0.01) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend()+ ylim(0, 50000)+ labs(title = "nCount_RNA_max50k") 
pl[["nCount_RNA_max5000"]] = VlnPlot(seur, features = "nCount_RNA", group.by = "library", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  +
  NoLegend()+ ylim(0, 5000)+labs(title = "nCount_RNA_max5k")
pl[["nFeature_RNA"]] = VlnPlot(seur, features = "nFeature_RNA", group.by = "library", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + 
  NoLegend()+ labs(title = "nFeature_RNA")
pl[["percent.mt"]] = VlnPlot(seur, features = "percent.mt", group.by = "library", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + 
  NoLegend()+ labs(title = "percent.mt")
pl[["percent.mt_max5"]] = VlnPlot(seur, features = "percent.mt", group.by = "library", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  NoLegend()+ylim(0,5)+ labs(title = "percent.mt_max5%")


pdf(file = paste0(out_dir,script_ind,"QC_plots.pdf"), width = 10, height = 6)
lapply(pl, function(x){x} )
dev.off()


message("\n\n##########################################################################\n",
        "# Completed B01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

