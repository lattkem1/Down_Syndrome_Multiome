### load RNA data, map to reference dataset, add ATAC data

message("\n#############################################################\n",
        "####### Start G02 Map data to reference (including subset reference annotation): ", Sys.time(),
        "\n#############################################################\n")

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
script_ind = "G02_"

#specify output directory
out_dir = paste0(main_dir,"G_basic_analysis_grafts_map_to_tissue_v05_10X_Singl_DS2U_DS1/")

### load normalised Seurat reference dataset

load(file = "C_subsetting_all_cells_non_cx_excl/C03_seur_integr_labelled.rda") 
seur_ref = seur
rm(seur)
gc()

gr_tab_ref = read_csv("B_basic_analysis/B02_gr_tab_filtered_non_cx_excl.csv")

clust_tab_ref = read_csv("C_subsetting_all_cells_non_cx_excl/C02_subset_cluster_assignment.csv")


### load own data

#load group and file info
gr_tab = read_csv("A_input/group_tab_grafts_DS2U_DS1_SCZ.csv")

#load own dataset (list of seurat objects)

load(file = paste0(out_dir,"G01_seur_merged.rda")) 


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

#distinct scale (larger number of colours, colour vector gets shuffled)
pal_dist = function(v){
  v2 = length(unique(v))
  if (v2 < 6){p2 = matlab.like(6)[1:v2]} else {
    p2 = matlab.like(v2)
    set.seed(12)
    p2 = sample(p2)
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
# filter own data
###########################################################

message("\n\n          *** Removing low quality cells... ", Sys.time(),"\n\n")

set.seed(1234)
DefaultAssay(seur) <- "RNA"

filter_stats = seur@meta.data %>% group_by(library) %>% summarise(N_cells_unfiltered = n())

set.seed(123)

seur = subset(x = seur,
              nCount_RNA < 30000 &
                nCount_RNA > 500 &
                percent.mt < 2
)


t1 = seur@meta.data %>% group_by(library) %>% summarise(N_high_quality = n())

#remove samples with <500 cells, identify dominant samples with many cells (>2*median remaining samples)

t1$N_kept_for_analysis = t1$N_high_quality
t1$N_kept_for_analysis[t1$N_high_quality<500] = 0

median_retained = median(t1$N_kept_for_analysis[t1$N_kept_for_analysis>0])

t1$N_kept_for_analysis[t1$N_high_quality>2*median_retained] = 2*median_retained

filter_stats$N_high_quality = t1$N_high_quality[match(filter_stats$library, t1$library)]
filter_stats$N_kept_for_analysis = t1$N_kept_for_analysis[match(filter_stats$library, t1$library)]

filter_stats$fract_low_quality = (filter_stats$N_cells_unfiltered - filter_stats$N_high_quality)/filter_stats$N_cells_unfiltered
filter_stats$fract_kept_for_analysis = filter_stats$N_kept_for_analysis/filter_stats$N_cells_unfiltered

write_csv(filter_stats, file = paste0(out_dir,script_ind, "filter_stats.csv"))

#save gr_tab with only retained libaries
libraries_kept = filter_stats$library[filter_stats$N_kept_for_analysis>0]
gr_tab = gr_tab[gr_tab$library %in% libraries_kept,]

write_csv(gr_tab, file = paste0(out_dir,script_ind,"gr_tab_filtered.csv"))


#subsample seurat to retain only cells/samples determined as above

t1 = seur@meta.data
t2 = t1[t1$library %in% libraries_kept,]

cells_keep = NULL

for (lib in libraries_kept){
  
  t3 = t2[t2$library==lib,]
  
  v1 = rownames(t3)
  
  if (length(v1)>filter_stats$N_kept_for_analysis[filter_stats$library == lib]){
    set.seed(123)
    v2 = sample(v1, size = filter_stats$N_kept_for_analysis[filter_stats$library == lib])
  } else {v2 = v1}
  cells_keep = c(cells_keep,v2)
  
}

seur = subset(seur, cells = cells_keep)




###########################################################
# map own dataset (seur_list) on reference
###########################################################

options(future.globals.maxSize = 30 * 1024^3) #SCTransform exceeds default memory limit for parallelisation with futures


#SCT normalise dataset

set.seed(123)

seur <- SCTransform(seur, ncells = 3000,variable.features.n = 2000, 
                        assay = "RNA",new.assay.name = "SCT", 
                        conserve.memory = TRUE,
                        verbose = TRUE)

#add pca model to reference for mapping query data
#seur_ref <- RunUMAP(seur_ref, dims = 1:30, reduction = "harmony", return.model = TRUE)


###mapping query data and add ref.umap and add percent.mt to metadata 

message("\n\n          *** Finding transfer anchors... ", Sys.time(),"\n\n")

anchors <- FindTransferAnchors(reference = seur_ref, query = seur,
                               dims = 1:30, reference.reduction = "pca")


message("\n\n          *** Mapping query data... ", Sys.time(),"\n\n")


seur <- MapQuery(anchorset = anchors, reference = seur_ref, query = seur,
                 refdata = list(cluster_name = "cluster_name", cell_type = "cell_type", 
                                dev_PCW = "dev_PCW", stage = "stage"), 
                 reference.reduction = "pca", reduction.model = "umap")


message("\n\n          *** Saving dataset with ref mapping... ", Sys.time(),"\n\n")

save(seur, file = paste0(out_dir, script_ind,"seur_with_ref_mapping.rda")) 




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
                                       "cluster_name", "group", "cell_type", "sample_type", "dev_PCW", "stage"))
pl_tab = cbind(t1, sample = "reference", 
               cluster_name_ref = t1$cluster_name, cell_type_ref = t1$cell_type, 
               protocol = "reference", seq_tech = "reference", isogenic_pair = "reference",
               predicted.cluster_name.score = 1, predicted.cluster_name = t1$cluster_name, 
               predicted.cell_type.score = 1, predicted.cell_type = t1$cell_type, 
               predicted.dev_PCW.score = 1, predicted.dev_PCW = t1$dev_PCW, 
               predicted.stage.score = 1, predicted.stage = t1$stage, 
               dataset = "reference")


#add plot data from graft datasets

t1 = seur@meta.data
colnames(t1)

t2 = FetchData(seur, vars = c("refUMAP_1", "refUMAP_2"))
names(t2) = c("umap_1", "umap_2")
t3 = cbind(t2, 
           cluster_name = "NA", group = t1$group, cell_type = "graft", sample_type = t1$sample_type, dev_PCW = "NA", stage = "NA",
           sample = t1$sample, 
           cluster_name_ref = "NEU_graft", cell_type_ref = "NEU_graft", 
           t1[,c("protocol", "seq_tech", "isogenic_pair", "predicted.cluster_name.score", "predicted.cluster_name",
                 "predicted.cell_type.score", "predicted.cell_type",
                 "predicted.dev_PCW.score", "predicted.dev_PCW",
                 "predicted.stage.score", "predicted.stage"
                 )],
           dataset = "grafts")

pl_tab = rbind(pl_tab, t3)


#convert to American English

pl_tab$sample_type[pl_tab$sample_type == "foetal"] = "fetal"



#############################################
# plot graft samples on reference UMAP
#############################################

message("\n\n          *** Plot ref UMAP mapping... ", Sys.time(),"\n\n")

#shuffle cells (reference separately) to avoid unbalanced overplotting 
t1 = pl_tab[pl_tab$sample != "reference", ]
set.seed(123)
t1 = t1[sample(nrow(t1)),]

t2 = pl_tab[pl_tab$sample == "reference", ]
set.seed(123)
t2 = t2[sample(nrow(t2)),]
pl_tab_shuffled = rbind(t1, t2)


#define grouping variables
datasets = unique(pl_tab$dataset)
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
sample_types = unique(pl_tab$sample_type)
cluster_names = unique(clust_tab_ref$cluster_name)
cell_types = unique(clust_tab_ref$cell_type)
stages = unique(seur_ref$stage[order(seur_ref$dev_PCW)])
seq_techs = unique(pl_tab$seq_tech)


##################################################################################################
# plot cells by split by ref vs query dataset, colored by predicted cluster_name and stage 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_minimal()

pl[["cluster_name"]] = p1 + geom_point(aes(color = predicted.cluster_name), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted cluster_name split by dataset")
pl[["stage"]] = p1 + geom_point(aes(color = predicted.stage), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted dev_PCW split by dataset")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_dataset.pdf"), 
    width = 7, height = 3)
lapply(pl, function(x){x} )
dev.off()


pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_minimal()

pl[["cluster_name"]] = p1 + geom_point(aes(color = predicted.cluster_name), size = 5, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted cluster_name split by dataset")
pl[["stage"]] = p1 + geom_point(aes(color = predicted.stage), size = 5, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$dataset, levels = datasets), nrow = 1)+
  labs(title = "predicted dev_PCW split by dataset")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_dataset_larger_points.pdf"), 
    width = 11, height = 5)
lapply(pl, function(x){x} )
dev.off()




##################################################################################################
# plot cells by split by sequencing technology, colored by predicted cluster_name and stage 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_minimal()

pl[["cluster_name"]] = p1 + geom_point(aes(color = predicted.cluster_name), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$seq_tech, levels = seq_techs), nrow = 1)+
  labs(title = "predicted cluster_name split by sequencing technology")
pl[["stage"]] = p1 + geom_point(aes(color = predicted.stage), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$seq_tech, levels = seq_techs), nrow = 1)+
  labs(title = "predicted dev_PCW split by sequencing technology")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_seq_tech.pdf"), 
    width = 16, height = 5)
lapply(pl, function(x){x} )
dev.off()




##################################################################################################
# plot cells by split by group vs sample_type, colored by predicted cluster_name and stage 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_minimal()

pl[["cluster_name"]] = p1 + geom_point(aes(color = predicted.cluster_name), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  facet_grid(rows = vars(factor(pl_tab_shuffled$sample_type, levels = sample_types)),
             cols = vars(factor(pl_tab_shuffled$group, levels = gr)))+
  labs(title = "predicted cluster_name split by group vs sample type")
pl[["stage"]] = p1 + geom_point(aes(color = predicted.stage), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_grid(rows = vars(factor(pl_tab_shuffled$sample_type, levels = sample_types)),
             cols = vars(factor(pl_tab_shuffled$group, levels = gr)))+
  labs(title = "predicted dev_PCW split by group vs sample type")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_group_vs_sample_type.pdf"), 
    width = 8, height = 7)
lapply(pl, function(x){x} )
dev.off()

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_minimal()

pl[["cluster_name"]] = p1 + geom_point(aes(color = predicted.cluster_name), size = 5, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  facet_grid(rows = vars(factor(pl_tab_shuffled$sample_type, levels = sample_types)),
             cols = vars(factor(pl_tab_shuffled$group, levels = gr)))+
  labs(title = "predicted cluster_name split by group vs sample type")
pl[["stage"]] = p1 + geom_point(aes(color = predicted.stage), size = 5, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_grid(rows = vars(factor(pl_tab_shuffled$sample_type, levels = sample_types)),
             cols = vars(factor(pl_tab_shuffled$group, levels = gr)))+
  labs(title = "predicted dev_PCW split by group vs sample type")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_group_vs_sample_type_larger_points.pdf"), 
    width = 11, height = 10)
lapply(pl, function(x){x} )
dev.off()





##################################################################################################
# plot cells by split by ref vs samples, colored by predicted cluster_name and stage 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_minimal()

pl[["cluster_name"]] = p1 + geom_point(aes(color = predicted.cluster_name), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("ref", samples)), nrow = 1)+
  labs(title = "predicted cluster_name split by sample")
pl[["stage"]] = p1 + geom_point(aes(color = predicted.stage), size = 0.5, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("ref", samples)), nrow = 1)+
  labs(title = "predicted dev_PCW split by sample")

pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_sample.pdf"), 
    width = length(samples)*5+1, height = 5)
lapply(pl, function(x){x} )
dev.off()



#################################################
#generate and plot stats by sample vs predicted cluster_name
#################################################

message("\n\n          *** Plot predicted cluster by sample... ", Sys.time(),"\n\n")

t1 = pl_tab[pl_tab$sample != "reference",]

t2 = t1 %>% group_by(sample, predicted.cluster_name) %>% dplyr::summarise(N_cells = n())

t3 = t1 %>% group_by(sample) %>% dplyr::summarise(N_cells = n())

t2$N_cells_sample = t3$N_cells[match(t2$sample, t3$sample)]
t2$fract_sample = t2$N_cells/t2$N_cells_sample

stat_pred_cluster_name = t2

write_csv(stat_pred_cluster_name , file = paste0(out_dir, script_ind, "stat_pred_cluster_name.csv"))


### generate and plot cell cluster vs pred cell type matrix as heatmap

t1 = stat_pred_cluster_name

predicted.cluster_names = intersect(cluster_names, t1$predicted.cluster_name)

pred_mat = matrix(nrow = length(samples), ncol = length(predicted.cluster_names),
                       dimnames = list(samples, predicted.cluster_names))

for (s1 in samples){
  
  t2 = t1[t1$sample == s1,]
  pred_mat[s1,] = t2$fract_sample[match(colnames(pred_mat),t2$predicted.cluster_name)]
  
}

pred_mat[is.na(pred_mat)] = 0


### plot heatmap of proportion of cells/cluster predicted as cell type from reference

pdf(file = paste0(out_dir,script_ind, "stat_pred_cluster_name_heatmap.pdf"), 
    width = 6, height = 8)
{
  heatmap_vir(pred_mat, main = "cluster_name predicted from reference vs sample", 
              cluster_rows = FALSE, cluster_cols = FALSE)
  heatmap_vir(pred_mat, main = "cluster_name predicted from reference vs sample", 
              cluster_rows = FALSE, cluster_cols = TRUE)
  heatmap_vir(pred_mat, main = "cluster_name predicted from reference vs sample", 
              cluster_rows = TRUE, cluster_cols = TRUE)
}
dev.off()








#get info on version of R, used packages etc
sessionInfo()





message("\n####### End G03a: ", Sys.time(),"\n")




