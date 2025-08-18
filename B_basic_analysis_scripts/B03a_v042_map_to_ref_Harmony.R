### load RNA data, map to reference dataset, add ATAC data

message("\n#################################\n",
        "####### Start B03a Map data to reference: ", Sys.time(),
        "\n##################################\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E14_241219_DS_foetal_brain_grafts_for_man_v02_low_string/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(Signac)
library(colorRamps)
library(viridis)
library(pheatmap)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#specify script/output index as prefix for file names
script_ind = "B03a_"

#specify output directory
out_dir = paste0(main_dir,"B_basic_analysis/")


### load normalised Seurat reference dataset

#define directory for reference dataset
ref_dir = "/rds/general/user/mlattke/projects/dsseq23/live/E02_230608_foetal_brain_ref_atlas/"
load(file = paste0(ref_dir,"/R_out_Seur5/R03_ref_seur_integr.rda") )

gr_tab_ref = read_csv(file =paste0(ref_dir, "R_in/gr_tab_ref.csv"))

#define order of areas and cell_types

t1 = seur_ref@meta.data

t2 = t1 %>% group_by(structure, area) %>% summarise(N = n())
t3 = t2[grepl("cortex", t2$structure),]

predicted.areas = c(unique(t3$area), unique(t2$area[!(t2$area %in% t3$area)]))

predicted.cell_types = unique(t1$cell_type)


### load own data

#load group and file info
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")

#load raw dataset
load(file = paste0(out_dir,"B02_seur_integr.rda")) 

#classify stages into early (<= PCW12, till end layer 6), mid (PCW13-16, end of layer 2/3), late (>=PCW17)
t1 = seur@meta.data
t1$stage = "PCW11_12"
t1$stage[t1$dev_PCW > 12] = "PCW13_16"
t1$stage[t1$dev_PCW > 16] = "PCW17_20"
seur@meta.data = t1



###########################################################
# functions
###########################################################

#custom colour palette for variable values defined in vector v
pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey20", "dodgerblue")
  } else if (v2 ==3){
    p2 = c("dodgerblue", "grey20", "orange")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}

### heatmap function with viridis (magma) colouring

heatmap_vir = function(m1, main = "", cluster_rows = FALSE, cluster_cols = FALSE){
  
  p1 = pheatmap::pheatmap(m1, show_rownames=TRUE, cluster_rows = cluster_rows,
                          cluster_cols = cluster_cols, show_colnames = TRUE, 
                          clustering_distance_rows = "euclidean",
                          clustering_method = "ward.D2",
                          treeheight_row = 50,
                          color = viridis_pal(option = "magma")(250),
                          breaks = seq(0, max(m1), length.out = 251),
                          border_color = NA, fontsize = 10,
                          cellwidth = 10, cellheight = 10,
                          main = main
  )
}



###########################################################
# map own dataset (seur_list) on reference
###########################################################

options(future.globals.maxSize = 30 * 1024^3) #SCTransform exceeds default memory limit for parallelisation with futures

#add pca model to reference for mapping query data
seur_ref <- RunUMAP(seur_ref, dims = 1:30, reduction = "harmony", return.model = TRUE)


###mapping query data and add ref.umap and add percent.mt to metadata 

message("\n\n          *** Finding transfer anchors... ", Sys.time(),"\n\n")

anchors <- FindTransferAnchors(reference = seur_ref, query = seur,
                               dims = 1:30, reference.reduction = "pca")


message("\n\n          *** Mapping query data... ", Sys.time(),"\n\n")


seur <- MapQuery(anchorset = anchors, reference = seur_ref, query = seur,
                 refdata = list(cell_cluster = "cell_cluster", cell_type = "cell_type",
                                area = "area", age = "age"), 
                 reference.reduction = "pca", reduction.model = "umap")


message("\n\n          *** Saving dataset with ref mapping... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir, script_ind,"seur_with_ref_mapping.rda")) 




#################################################
#plot predicted cell type and area on own dataset UMAP
#################################################

#load(file = paste0(out_dir, script_ind,"seur_with_ref_mapping.rda")) 

message("\n\n          *** plot predictions on UMAP... ", Sys.time(),"\n\n")


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


### generate umap plot by group and sample

#create lists for plot objects
pl_umap = list()

#define grouping variables and colors for plotting
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
libraries = unique(gr_tab$library)
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

pl_umap[["library"]] = DimPlot(seur, group.by = "library", shuffle = TRUE,
                               label = FALSE, reduction = "umap", 
                               pt.size = 0.01)+
  scale_color_manual(limits = libraries, values = pal(libraries))+
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

pl_umap[["predicted.area"]] = DimPlot(seur, group.by = "predicted.area", shuffle = TRUE,
                             label = FALSE, reduction = "umap", 
                             pt.size = 1)+
  scale_color_manual(limits = predicted.areas, values = pal(predicted.areas))+
  labs(title = "UMAP by predicted.area")

pl_umap[["predicted.area_noLegend"]] = DimPlot(seur, group.by = "predicted.area", shuffle = TRUE,
                                      label = FALSE, reduction = "umap", 
                                      pt.size = 0.01)+NoLegend()+
  scale_color_manual(limits = predicted.areas, values = pal(predicted.areas))+
  labs(title = "UMAP by predicted.area")

pl_umap[["predicted.cell_type"]] = DimPlot(seur, group.by = "predicted.cell_type", shuffle = TRUE,
                                      label = FALSE, reduction = "umap", 
                                      pt.size = 1)+
  scale_color_manual(limits = predicted.cell_types, values = pal(predicted.cell_types))+
  labs(title = "UMAP by predicted.cell_type")

pl_umap[["predicted.cell_type_NoLegend"]] = DimPlot(seur, group.by = "predicted.cell_type", shuffle = TRUE,
                                           label = FALSE, reduction = "umap", 
                                           pt.size = 0.01)+NoLegend()+
  scale_color_manual(limits = predicted.cell_types, values = pal(predicted.cell_types))+
  labs(title = "UMAP by predicted.cell_type")

pdf(file = paste0(out_dir,script_ind, "UMAP_covars_ref_predictions.pdf"), width = 6, height = 5)
lapply(pl_umap, function(x){x})
dev.off()

pdf(file = paste0(out_dir,script_ind, "UMAP_covars_ref_predictions_smaller.pdf"), width = 4, height = 3.5)
lapply(pl_umap, function(x){x})
dev.off()


#############################################
# create table to plot reference vs query data
#############################################

message("\n\n          *** Extract ref mapping... ", Sys.time(),"\n\n")

#load(file = paste0(out_dir, script_ind,"seur_with_ref_mapping.rda")) 

#subsample reference dataset for smaller plot sizes
if (nrow(seur_ref@meta.data)>50000){
  seur_ref_plot = seur_ref[, sample(colnames(seur_ref), size =50000, replace=F)]
} else {seur_ref_plot = seur_ref}

t1 = seur_ref_plot@meta.data
colnames(t1)

#extract plot data from reference dataset + add query dataset-specific metadata columns (set sample as "ref")

t1 = FetchData(seur_ref_plot, vars = c("umap_1", "umap_2", 
                                       "area", "age", "cell_type"))
pl_tab = cbind(t1, sample = "reference", group = "reference",  
               dev_PCW = t1$age-2, cell_type_ref = t1$cell_type, 
               predicted.cell_type.score = 1, predicted.cell_type = t1$cell_type,         
               predicted.area.score = 1, predicted.area = t1$area, dataset = "reference")



#add plot data from graft datasets

t1 = seur@meta.data
colnames(t1)

t2 = FetchData(seur, vars = c("refUMAP_1", "refUMAP_2"))
names(t2) = c("umap_1", "umap_2")
t3 = cbind(t2, area = "this_study", age = t1$dev_PCW+2, cell_type = "this_study", 
           sample = t1$sample, group = t1$group, dev_PCW = t1$dev_PCW, 
           cell_type_ref = "this_study", 
           predicted.cell_type.score = t1$predicted.cell_type.score, 
           predicted.cell_type = t1$predicted.cell_type,         
           predicted.area.score = t1$predicted.area.score, predicted.area = t1$predicted.area, 
           dataset = "this_study")

pl_tab = rbind(pl_tab, t3)


#############################################
# plot graft samples on reference UMAP
#############################################

message("\n\n          *** Plot ref UMAP mapping... ", Sys.time(),"\n\n")

#shuffle cells (except reference) to avoid unbalanced overplotting (keep only 50k cells from own data to reduce plot size)
t1 = pl_tab[pl_tab$sample != "reference", ]
set.seed(123)
t1 = t1[sample(nrow(t1), 50000),]
pl_tab_shuffled = rbind(pl_tab[pl_tab$sample == "reference",], t1)

#define grouping variables
datasets = unique(pl_tab$dataset)
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)



### plot graft cells by samples and covariates on reference UMAP

set.seed(123)
pl_tab_shuffled  = pl_tab_shuffled [sample(nrow(pl_tab_shuffled )),]

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_light()

pl[["cell_type_ref"]] = p1 + geom_point(aes(color = cell_type_ref), size = 2, alpha =0.7)+
  scale_color_manual(limits = c("this_study", predicted.cell_types), values = c("black", pal(predicted.cell_types)))+
  labs(title = "Dataset mapped on reference cell types")

pl[["cell_type_ref_no_legend"]] = p1 + geom_point(aes(color = cell_type_ref), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("this_study", predicted.cell_types), values = c("black", pal(predicted.cell_types)))+
  labs(title = "Dataset mapped on reference cell types")+NoLegend()

pl[["sample"]] = p1 + geom_point(aes(color = sample), size =2, alpha =0.7)+
  scale_color_manual(limits = c("reference", samples), values = c("grey", pal(samples)))+
  labs(title = "Samples mapped on reference")

pl[["sample_no_legend"]] = p1 + geom_point(aes(color = sample), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("reference", samples), values = c("grey", pal(samples)))+
  labs(title = "Samples mapped on reference")+NoLegend()

pl[["group"]] = p1 + geom_point(aes(color = group), size = 2, alpha =0.7)+
  scale_color_manual(limits = c("reference", gr), values = c("grey", pal(gr)))+
  labs(title = "Group mapped on reference")

pl[["group_no_legend"]] = p1 + geom_point(aes(color = group), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("reference", gr), values = c("grey", pal(gr)))+
  labs(title = "Group mapped on reference")+NoLegend()

pl[["area"]] = p1 + geom_point(aes(color = area), size = 2, alpha =0.7)+
  scale_color_manual(limits = c("this_study", predicted.areas), values = c("black", pal(predicted.areas)))+
  labs(title = "Dataset mapped on reference areas")

pl[["area_no_legend"]] = p1 + geom_point(aes(color = area), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("this_study", predicted.areas), values = c("black", pal(predicted.areas)))+
  labs(title = "Dataset mapped on reference areas")+NoLegend()

pdf(file = paste0(out_dir,script_ind, "UMAP_samples_covars_mapped_on_reference.pdf"), 
    width = 6, height = 5)
lapply(pl, function(x){x} )
dev.off()




##################################################################################################
# plot cells by split by ref vs query dataset, colored by predicted cell type and area 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_bw()

pl[["cell_types"]] = p1 + geom_point(aes(color = predicted.cell_type), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = predicted.cell_types, values = pal(predicted.cell_types))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted cell type split by dataset")
pl[["area"]] = p1 + geom_point(aes(color = predicted.area), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = predicted.areas, values = pal(predicted.areas))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted area split by dataset")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_dataset.pdf"), 
    width = 7.5, height = 2.5)
lapply(pl, function(x){x} )
dev.off()


pl[["cell_types"]] = p1 + geom_point(aes(color = predicted.cell_type), size = 5, alpha = 0.7)+
  scale_color_manual(limits = predicted.cell_types, values = pal(predicted.cell_types))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted cell type split by dataset")
pl[["area"]] = p1 + geom_point(aes(color = predicted.area), size = 5, alpha = 0.7)+
  scale_color_manual(limits = predicted.areas, values = pal(predicted.areas))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted area split by dataset")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_dataset_larger_points.pdf"), 
    width = 12, height = 4)
lapply(pl, function(x){x} )
dev.off()





##################################################################################################
# plot cells by split by sample, colored by predicted cell type and area 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_bw()

pl[["cell_types"]] = p1 + geom_point(aes(color = predicted.cell_type), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = predicted.cell_types, values = pal(predicted.cell_types))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("reference", samples)), nrow = 1)+
  labs(title = "predicted cell type split by sample")
pl[["area"]] = p1 + geom_point(aes(color = predicted.area), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = predicted.areas, values = pal(predicted.areas))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("reference", samples)), nrow = 1)+
  labs(title = "precited area split by sample")


pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_sample.pdf"), 
    width = 2.5*length(samples)+4, height = 2.5)
lapply(pl, function(x){x} )
dev.off()




#################################################
#generate and plot stats by sample vs predicted area
#################################################

t1 = pl_tab[pl_tab$sample != "reference",]

t2 = t1 %>% group_by(sample, predicted.area) %>% dplyr::summarise(N_cells = n())

t3 = t1 %>% group_by(sample) %>% dplyr::summarise(N_cells = n())

t2$N_cells_sample = t3$N_cells[match(t2$sample, t3$sample)]
t2$fract_sample_predicted_area = t2$N_cells/t2$N_cells_sample

stat_area_pred = t2

write_csv(stat_area_pred , file = paste0(out_dir, script_ind, "stat_area_pred.csv"))


### generate and plot cell cluster vs pred cell type matrix as heatmap

area_pred_mat = matrix(nrow = length(predicted.areas), ncol = length(samples),
                       dimnames = list(predicted.areas, samples))

for (s1 in samples){
  
  t1 = stat_area_pred[stat_area_pred$sample == s1,]
  
  area_pred_mat[t1$predicted.area, s1] = t1$fract_sample_predicted_area
  
}

area_pred_mat[is.na(area_pred_mat)] = 0


### plot heatmap of proportion of cells/cluster predicted as cell type from reference

pdf(file = paste0(out_dir,script_ind, "stat_area_pred_heatmap.pdf"), 
    width = 10, height = 8)
{
  heatmap_vir(area_pred_mat, main = "area predicted from reference vs sample", 
              cluster_rows = FALSE, cluster_cols = FALSE)
  heatmap_vir(area_pred_mat, main = "area predicted from reference vs sample", 
              cluster_rows = TRUE, cluster_cols = FALSE)
  heatmap_vir(area_pred_mat, main = "area predicted from reference vs sample", 
              cluster_rows = TRUE, cluster_cols = TRUE)
}
dev.off()





#get info on version of R, used packages etc
sessionInfo()





message("\n####### End B03: ", Sys.time(),"\n")


