message("\n\n##########################################################################\n",
        "# Start B03: Clustering tests and basic cell population analysis ", Sys.time(),
        "\n##########################################################################\n",
        "\n   \n", 
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
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
script_ind = "B03_"

#specify output directory
out_dir = paste0(main_dir,"B_basic_analysis/")

#load group and file info
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")

#load raw dataset
load(file = paste0(out_dir,"B02_seur_integr.rda")) 

#get marker gene panels
GOI = list()
t1 = read_csv("A_input/cell_type_markers_240806_consolidated.csv")
GOI$cell_type_markers = t1$gene[t1$level %in% c("cell_types", "neuronal_lineage")]
GOI$subtype_markers = t1$gene[t1$level %in% c("NPC_patterning","cortical_layers", "Interneuron_subtypes")]

#classify stages into early (<= PCW12, till end layer 6), mid (PCW13-16, end of layer 2/3), late (>=PCW17)
t1 = seur@meta.data
t1$stage = "PCW11_12"
t1$stage[t1$dev_PCW > 12] = "PCW13_16"
t1$stage[t1$dev_PCW > 16] = "PCW17_20"
seur@meta.data = t1


###########################################################
#functions
###########################################################

#custom colour palette for variable values defined in vector v
pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey", "blue")
  } else if (v2 ==3){
    p2 = c("blue", "grey", "orange")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}


#############################################
#Plot UMAP and marker expression (test different clustering resolutions)
#############################################

set.seed(1234)

#define resolutions for clustering tests
test_res = c(0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5)


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


#create lists for plot objects
pl_umap = list()
pl_dotplots_cell_type_markers = list()
pl_dotplots_subtype_markers = list()
pl_umap_expr = list()


### generate umap plot by group and sample


#define grouping variables and colors for plotting
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
dev_PCW_stages = unique(gr_tab$dev_PCW[order(gr_tab$dev_PCW)])
batches = unique(gr_tab$batch_seq[order(gr_tab$batch_seq)])
stages = unique(seur$stage[order(seur$stage)])


pl_umap[["group"]] = DimPlot(seur, group.by = "group", shuffle = TRUE,
                             label = FALSE, reduction = "umap", 
                             pt.size = 0.01)+
  scale_color_manual(limits = gr, values = pal(gr))+
  labs(title = "UMAP by group")

pl_umap[["sample"]] = DimPlot(seur, group.by = "sample", shuffle = TRUE,
                              label = FALSE, reduction = "umap", 
                              pt.size = 0.01)+
  scale_color_manual(limits = samples, values = pal(samples))+
  labs(title = "UMAP by sample")

pl_umap[["dev_PCW"]] = DimPlot(seur, group.by = "dev_PCW", shuffle = TRUE,
                              label = FALSE, reduction = "umap", 
                              pt.size = 0.01)+
  scale_color_manual(limits = factor(dev_PCW_stages), values = pal(dev_PCW_stages))+
  labs(title = "UMAP by PCW")

pl_umap[["batch_seq"]] = DimPlot(seur, group.by = "batch_seq", shuffle = TRUE,
                               label = FALSE, reduction = "umap", 
                               pt.size = 0.01)+
  scale_color_manual(limits = batches, values = pal(batches))+
  labs(title = "UMAP by batch")

pl_umap[["stage"]] = DimPlot(seur, group.by = "stage", shuffle = TRUE,
                               label = FALSE, reduction = "umap", 
                               pt.size = 0.01)+
  scale_color_manual(limits = stages, values = pal(stages))+
  labs(title = "UMAP by stage")


### create plots for clusters and marker dot plot for different resolutions (each test_res)
j = test_res[1]
for (j in test_res){
  
  set.seed(1234)
  
  seur0 <- FindClusters(
    object = seur0,
    graph.name = "SCT_snn",
    algorithm = 1,
    resolution = j,
    verbose = FALSE
  )
  
  #subset for plotting
  
  seur = seur0
  
  cells = colnames(seur@assays$SCT@scale.data)
  
  if (length(cells)>100000){
    set.seed(1234)
    seur = seur[, sample(cells, size =100000, replace=F)]
  } 
  
  cl = unique(seur$seurat_clusters)
  pl_umap[[paste0("resol_", j)]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE,
                                           label = TRUE, reduction = "umap", pt.size = 0.01)+
    scale_color_manual(limits = cl, values = pal(cl))+
    NoLegend()+labs(title = paste0("clusters resol_", j))
  
  
  #create dotplot with gene expression of marker genes
  
  DefaultAssay(seur0) = "SCT"
  
  pl_dotplots_cell_type_markers[[paste0("resol_", j)]] = 
    DotPlot(seur0, features = intersect(GOI$cell_type_markers, rownames(seur0)), 
            scale.by = "size") +RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("resol_", j))
  
  pl_dotplots_subtype_markers[[paste0("resol_", j)]] = 
    DotPlot(seur0, features = intersect(GOI$subtype_markers, rownames(seur0)), 
            scale.by = "size") +RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("resol_", j))
  
}


seur = seur0

#save dataset with cluster assignment for different cluster resolutions
save(seur, file = paste0(out_dir, script_ind, "seur_integr.rda")) 


### save umap clustering and dotplots

pdf(file = paste0(out_dir, script_ind, "UMAP_clustering_test_integr_dataset.pdf"), width = 6, height = 5)
lapply(pl_umap, function(x){x})
dev.off()

pdf(file = paste0(out_dir, script_ind, "Dotplots_clustering_test_integr_dataset_cell_type_markers.pdf"), width = 12, height = 10)
lapply(pl_dotplots_cell_type_markers, function(x){x})
dev.off()

pdf(file = paste0(out_dir, script_ind, "Dotplots_clustering_test_integr_dataset_subtype_markers.pdf"), width = 12, height = 10)
lapply(pl_dotplots_subtype_markers, function(x){x})
dev.off()



############################################################
# UMAP plots split by group vs repl/batch_seq/dev_PCW/stage
############################################################

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


t1 = FetchData(seur, vars = c("umap_1", "umap_2", "group", "sample", "stage", "dev_PCW", "batch_seq"))

#plot group vs stage

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2))+geom_point(size = 0.1, alpha = 0.5)+
  ggplot2::facet_grid(cols = vars(stage), rows = vars(group))+
  theme_bw()+
  labs(title = "UMAP by group vs stage")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_stage.pdf"), width = 6, height = 3)
plot(p1)
dev.off()

#plot group vs batch

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2))+geom_point(size = 0.1, alpha = 0.5)+
  ggplot2::facet_grid(cols = vars(batch_seq), rows = vars(group))+
  theme_bw()+
  labs(title = "UMAP by group vs batch")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_batch.pdf"), width = 15, height = 4)
plot(p1)
dev.off()

#plot group vs PCW

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2))+geom_point(size = 0.1, alpha = 0.5)+
  ggplot2::facet_grid(cols = vars(dev_PCW), rows = vars(group))+
  theme_bw()+
  labs(title = "UMAP by group vs PCW")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_PCW.pdf"), width = 15, height = 4)
plot(p1)
dev.off()


#plot group vs sample

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2))+geom_point(size = 0.05, alpha = 0.5)+
  ggplot2::facet_wrap(facets = factor(t1$sample, levels = samples), nrow = 2)+
  theme_bw()+
  labs(title = "UMAP by group vs sample")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_sample.pdf"), width = 10, height = 2.5)
plot(p1)
dev.off()


######################################################
# plots for marker expression on UMAP
######################################################

#subsample seurat object for large datasets (else plots become too large (pdf with hundreds of MB))

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>10000){
  set.seed(1234)
  seur = seur[, sample(cells, size =10000, replace=F)]
} 

DefaultAssay(seur) = "SCT"

pl_umap_expr = FeaturePlot(seur, features = c(GOI$cell_type_markers, GOI$subtype_markers), order = TRUE,
                           reduction = "umap", pt.size = 0.01, ncol = 9)

pdf(file = paste0(out_dir, script_ind, "marker_expr_UMAP_integr_dataset.pdf"), width = 30, height = 15)
plot(pl_umap_expr)
dev.off()


###################################################################################
# save table for cluster assignment/labelling (with highest resolution) 
###################################################################################
# for downstream analyses assign cell types in csv; remove not needed clusters for lower resolutions

seur = seur0

v1 = unique(seur$seurat_clusters)
v2 = v1[order(v1)]
t1 = tibble(cluster = v2, cluster_name = paste0("s",v2), cell_type = "all", cell_class = "all")

if (file.exists(paste0(out_dir, script_ind, "cluster_assignment.csv"))){
  write_csv(t1, paste0(out_dir, script_ind, "cluster_assignment_new.csv"))
} else {write_csv(t1, paste0(out_dir, script_ind, "cluster_assignment.csv"))}




message("\n\n##########################################################################\n",
        "# Completed B03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")




