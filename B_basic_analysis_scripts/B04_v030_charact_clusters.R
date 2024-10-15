message("\n\n##########################################################################\n",
        "# Start B05: Clustering tests and basic cell population analysis ", Sys.time(),
        "\n##########################################################################\n",
        "\n   clustering with selected cluster resolution ",
        "\n   differential abundance analysis on cluster level", 
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
library(GenomicRanges)
library(colorRamps)
library(viridis)
library(lmerTest)
library(pheatmap)
library(sccomp)

#specify script/output index as prefix for file names
script_ind = "B04_"

#specify output directory
out_dir = paste0(main_dir,"B_basic_analysis/")

#load group and file info
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")

#load dataset
load(file = paste0(out_dir, "B03_seur_integr.rda")) 

#get marker gene panels
GOI = list()
t1 = read_csv("A_input/cell_type_markers_240806_consolidated.csv")
GOI$cell_type_markers = t1$gene[t1$level %in% c("cell_types", "neuronal_lineage")]
GOI$subtype_markers = t1$gene[t1$level %in% c("NPC_patterning","cortical_layers", "Interneuron_subtypes")]

#define clustering resolution for downstream analyses (for full dataset)
clust_resol = 0.5

#load cluster labels
clust_tab = read_csv(paste0(out_dir, "B03_cluster_assignment_all.csv"))


####################################
#Functions
####################################

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


###########################################################
# save cluster/cell class labels to dataset, plot with cluster/cell_class names
###########################################################

DefaultAssay(seur) = "SCT"

#add clusters with default resolution
seur$seurat_clusters = seur@meta.data[,paste0("SCT_snn_res.", clust_resol)]

seur$cluster_name = clust_tab$cluster_name[match(seur$seurat_clusters, clust_tab$cluster)]
seur$cell_type = clust_tab$cell_type[match(seur$seurat_clusters, clust_tab$cluster)]
seur$cell_class = clust_tab$cell_class[match(seur$seurat_clusters, clust_tab$cluster)]

save(seur, file = paste0(out_dir,script_ind, "seur_integr_labelled.rda")) 


### plot whole dataset as UMAP with named clusters/cell classes

#define grouping variables
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
clusters = clust_tab$cluster
cluster_names = clust_tab$cluster_name
cell_types = unique(clust_tab$cell_type)
cell_classes = unique(clust_tab$cell_class)
stages = unique(seur$stage[order(seur$stage)])


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


#umap plots for cell cluster and cell class (labelled)
pl = list()

pl[["group"]] = DimPlot(seur, group.by = "group", shuffle = TRUE,
                             label = FALSE, reduction = "umap", 
                             pt.size = 0.01)+
  scale_color_manual(limits = gr, values = pal(gr))+
  labs(title = "UMAP by group")

pl[["stage"]] = DimPlot(seur, group.by = "stage", shuffle = TRUE,
                             label = FALSE, reduction = "umap", 
                             pt.size = 0.01)+
  scale_color_manual(limits = stages, values = pal(stages))+
  labs(title = "UMAP by stage")

pl[["cluster"]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = clusters, values = pal(clusters))+
  NoLegend()+labs(title = "clusters")

pl[["cluster_dist"]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE, repel = TRUE, 
                          label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = clusters, values = pal_dist(clusters))+
  NoLegend()+labs(title = "clusters")

pl[["cluster_name_umap"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pl[["cluster_name_umap_dist"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pl[["cell_type_umap"]] = DimPlot(seur, group.by = "cell_type", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_types, values = pal(cell_types))+
  NoLegend()+labs(title = "cell_types")

pl[["cell_type_umap_dist"]] = DimPlot(seur, group.by = "cell_type", shuffle = TRUE, repel = TRUE, 
                                  label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_types, values = pal_dist(cell_types))+
  NoLegend()+labs(title = "cell_types")

pl[["cell_class_umap"]] = DimPlot(seur, group.by = "cell_class", shuffle = TRUE, repel = TRUE, 
                                  label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_classes, values = pal(cell_classes))+
  NoLegend()+labs(title = "cell_classes")

pl[["cell_class_umap_dist"]] = DimPlot(seur, group.by = "cell_class", shuffle = TRUE, repel = TRUE, 
                                label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_classes, values = pal_dist(cell_classes))+
  NoLegend()+labs(title = "cell_classes")

pdf(file = paste0(out_dir,script_ind, "UMAP_clusters_labelled.pdf"), width = 6, height = 5)
lapply(pl, function(x){x})
dev.off()


#labelled dotplot for marker gene expression by cluster (split by cell type markers vs neuron subtype markers)

seur = seur0

DefaultAssay(seur) = "SCT"

p1 = DotPlot(seur, features = intersect(GOI$cell_type_markers, rownames(seur)), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

p2 = DotPlot(seur, features = intersect(GOI$subtype_markers, rownames(seur)), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir,script_ind, "Cell_markers_dotplot_clusters_labelled.pdf"), 
    width = 10, height = 8)
plot(p1)
plot(p2)
dev.off()



#######################################################################################
### add plots split by group vs repl/batch_seq/dev_PCW/stage colored by cluster_name
#######################################################################################

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


t1 = FetchData(seur, vars = c("umap_1", "umap_2", "cluster_name","group", "sample", "stage", "dev_PCW", "batch_seq"))

#plot group vs stage

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2, color = cluster_name))+geom_point(size = 0.1, alpha = 0.5)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  ggplot2::facet_grid(cols = vars(stage), rows = vars(group))+
  theme_bw()+
  labs(title = "UMAP by group vs stage")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_stage.pdf"), width = 6, height = 3)
plot(p1)
dev.off()

#plot group vs batch

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2, color = cluster_name))+geom_point(size = 0.1, alpha = 0.5)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  ggplot2::facet_grid(cols = vars(batch_seq), rows = vars(group))+
  theme_bw()+
  labs(title = "UMAP by group vs batch")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_batch.pdf"), width = 15, height = 4)
plot(p1)
dev.off()

#plot group vs PCW

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2, color = cluster_name))+geom_point(size = 0.1, alpha = 0.5)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  ggplot2::facet_grid(cols = vars(dev_PCW), rows = vars(group))+
  theme_bw()+
  labs(title = "UMAP by group vs PCW")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_PCW.pdf"), width = 15, height = 4)
plot(p1)
dev.off()


#plot group vs sample

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2, color = cluster_name))+geom_point(size = 0.05, alpha = 0.5)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
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


#############################################################################
# cluster quantification and differential abundance analysis
#############################################################################

seur = seur0

meta = seur@meta.data

# create table with one row for each cluster for each sample

stat_tab = tibble(cluster = unlist(lapply(cluster_names, rep, length.out = length(samples))), 
                  sample = rep(samples, length(cluster_names)))
stat_tab$group = gr_tab$group[match(stat_tab$sample, gr_tab$sample)] #add group assignment

#count cells per cluster per sample, add to stat_tab
t1 = meta %>% group_by(cluster_name, sample) %>% summarize(N_cells = n())
stat_tab$N_cells = t1$N_cells[match( paste0(stat_tab$cluster,stat_tab$sample), paste0(t1$cluster_name,t1$sample) )]
stat_tab$N_cells[is.na(stat_tab$N_cells)] = 0

#add total cells per sample and per cluster, fraction of cluster, fraction of sample
N_sample = stat_tab %>% group_by(sample) %>% summarize(N = sum(N_cells))
stat_tab$N_sample = N_sample$N[match(stat_tab$sample, N_sample$sample)]
stat_tab$fract_sample = stat_tab$N_cells / stat_tab$N_sample
N_cluster = stat_tab %>% group_by(cluster) %>% summarize(N = sum(N_cells))
stat_tab$N_cluster = N_cluster$N[match(stat_tab$cluster, N_cluster$cluster)]
stat_tab$fract_cluster = stat_tab$N_cells / stat_tab$N_cluster

write_csv(stat_tab, file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster.csv"))


### crossbar-dotplot quantification of cluster contribution fraction of sample

t1 = stat_tab
t2 = t1 %>% group_by(cluster, group) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster, y = mean_fract, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.3)+
  geom_errorbar(data = t2, aes(x = cluster,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               color = group),
                position = position_dodge(width = 0.5), width = 0.3, lwd = 0.2)+
  geom_point(data = t1, aes(x = cluster, y = fract_sample, color = group, shape = sample), 
             position = position_dodge(width = 0.5), size = 0.5, stroke = 0.3)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names)+
  scale_color_manual(limits = gr, values = pal(gr))+
  scale_fill_manual(limits = gr, values = pal(gr))+
  scale_shape_manual(limits = unique(stat_tab$sample), 
                     values = rep_len(c(1:4), length.out = length(unique(stat_tab$sample))))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster.pdf"), 
    width = 7, height = 4)
plot(p1)
dev.off()


### heatmap quantification of cluster contribution fraction of sample

heatmap_vir = function(m1, main = "", cluster_rows = FALSE, cluster_cols = FALSE){
  p1 = pheatmap(m1, show_rownames=TRUE, cluster_rows = cluster_rows,
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

t1 = stat_tab

m1 = matrix(nrow = length(cluster_names), ncol = length(samples), dimnames = list(cluster_names, samples))

for (i in samples){
  t2 = t1[t1$sample==i,]
  m1[,i] = t2$fract_sample[match(rownames(m1),t2$cluster)]
}
m1[is.na(m1)] = 0

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster_heatmap.pdf"), 
    width = 10, height = 8)
heatmap_vir(m1, main = "cell_abundance_by_area_and_sample")
heatmap_vir(m1, main = "cell_abundance_by_area_and_sample", 
            cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()



### statistical analysis t-tests by cluster

t1 = stat_tab

t2 = tibble(cluster = cluster_names, p = NA, padj = NA)

for (i in cluster_names){
  t4 = t1[t1$cluster == i,]
  t5 = pairwise.t.test(t4$fract_sample, t4$group, p.adjust.method = "none")
  t2$p[t2$cluster == i] = t5$p.value[1,1]
}

t2$padj = p.adjust(t2$p, method = "BH")

write_csv(t2, paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster_t-Test.csv"))


#############################################################################
# cluster quantification and differential abundance analysis accounting for developmental stage
#############################################################################

#define additional grouping variables

dev_PCW_stages = unique(gr_tab$dev_PCW[order(gr_tab$dev_PCW)])
stages = unique(seur$stage[order(seur$stage)])

meta = seur@meta.data

# create table with one row for each cluster for each sample

t1 = expand_grid(cluster_names, samples)
names(t1) = c("cluster_name", "sample")
t2 = gr_tab[match(t1$sample, gr_tab$sample),]
t2$stage = seur$stage[match(t1$sample, seur$sample)]
stat_tab = cbind(t1, t2[,-c(1:3)])

#count cells per cluster per sample, add to stat_tab
t1 = meta %>% group_by(cluster_name, sample) %>% summarize(N_cells = n())
stat_tab$N_cells = t1$N_cells[match( paste0(stat_tab$cluster_name,stat_tab$sample), 
                                     paste0(t1$cluster_name,t1$sample) )]
stat_tab$N_cells[is.na(stat_tab$N_cells)] = 0

#add total cells per sample and per cluster, fraction of cluster, fraction of sample
N_sample = stat_tab %>% group_by(sample) %>% summarize(N = sum(N_cells))
stat_tab$N_sample = N_sample$N[match(stat_tab$sample, N_sample$sample)]
stat_tab$fract_sample = stat_tab$N_cells / stat_tab$N_sample
N_cluster = stat_tab %>% group_by(cluster_name) %>% summarize(N = sum(N_cells))
stat_tab$N_cluster = N_cluster$N[match(stat_tab$cluster_name, N_cluster$cluster_name)]
stat_tab$fract_cluster = stat_tab$N_cells / stat_tab$N_cluster

write_csv(stat_tab, file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster_w_metadata.csv"))



#############################################################################
# plot cluster abundance vs metadata (PCW, stage, sex, for each cluster)
#############################################################################

pl_PCW = list()
pl_stage = list()
pl_sex = list()

t1 = stat_tab

cl = cluster_names[1]

for (cl in cluster_names){
  
  t2 = t1[t1$cluster_name == cl,]
  
  p1 = ggplot(t2, aes(y = fract_sample, color = group))+
    geom_hline(yintercept = 0)+
    scale_color_manual(limits = gr, values = pal(gr))+
    labs(title = cl)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  pl_PCW[[cl]] = p1+
    geom_smooth(aes(x = dev_PCW), method = "lm")+
    geom_point(aes(x = dev_PCW))
  
  
  pl_stage[[cl]] = p1+
    stat_summary(aes(x = stage, group = group, color = group), 
                 fun.data = "mean_se", size = 0.5,  
                 position = position_dodge(width = 0.2), geom = "errorbar")+
    stat_summary(aes(x = stage, group = group, color = group),
                 fun.data = "mean_se", size = 0.5, 
                 position = position_dodge(width = 0.2), geom = "line")+
    geom_point(aes(x = stage), position = position_dodge(width = 0.4))+
    scale_x_discrete(limits = stages)
  
  pl_sex[[cl]] = p1+
    stat_summary(aes(x = sex, group = group, color = group), 
                 fun.data = "mean_se", size = 0.5,  
                 position = position_dodge(width = 0.2), geom = "errorbar")+
    stat_summary(aes(x = sex, group = group, color = group),
                 fun.data = "mean_se", size = 0.5, 
                 position = position_dodge(width = 0.2), geom = "line")+
    geom_point(aes(x = sex), position = position_dodge(width = 0.2))
  
}


pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_cluster_by_sample_vs_PCW.pdf"), 
    width = 4, height = 3)
lapply(pl_PCW, function(x){x})
dev.off()

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_cluster_by_sample_vs_stage.pdf"), 
    width = 3, height = 3)
lapply(pl_stage, function(x){x})
dev.off()

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_cluster_by_sample_vs_sex.pdf"), 
    width = 2.5, height = 3)
lapply(pl_sex, function(x){x})
dev.off()



##############################################
#fit lmertest model to assess overall impact of group on fract of cells per cluster 
##############################################

#function for LRT test for comparing 2 lmer models

lmer_anova = function(form, form_red, dat, red = "form vs form_red"){
  lmm1 = lmer(form, dat)
  lmm1_red = lmer(form_red, dat)
  sum_tab = cbind(form = as.character(form)[3],
                  form_red = as.character(form_red)[3], 
                  anova(lmm1, lmm1_red, test = "LRT"), 
                  red = red)
}


#function for LRT test for comparing 2 lm models

lm_anova = function(form, form_red, dat, red = "form vs form_red"){
  lmm1 = lm(form, dat)
  lmm1_red = lm(form_red, dat)
  sum_tab = cbind(form = as.character(form)[3],
                  form_red = as.character(form_red)[3], 
                  anova(lmm1, lmm1_red, test = "LRT"), 
                  red = red)
}




#compare lmer models for group (with interactions) and cluster_name contribution to fraction of cells in sample

t1 = stat_tab
t1$dev_PCW_scaled = scale(t1$dev_PCW)

form1 = fract_sample ~ (1|cluster_name) + dev_PCW_scaled + (1|sample) +
  (1|group) + (1|group:cluster_name) + (1|group:dev_PCW_scaled) 
form1_no_group = fract_sample ~ (1|cluster_name) + dev_PCW_scaled + (1|sample)
form1_no_cluster = fract_sample ~ dev_PCW_scaled + (1|sample) +
  (1|group) + (1|group:dev_PCW_scaled) 

t2 = lmer_anova(form1, form1_no_cluster, dat = t1, red = "cluster_name_removed")
t3 = lmer_anova(form1, form1_no_group, dat = t1, red = "group_w_int_removed")
t4 = rbind(t2,t3)
t4 = t4[!is.na(t4$`Pr(>Chisq)`),]

write_csv(t4, file = paste0(out_dir,script_ind,"cell_abundance_clusters_combined_lmer_LRT_stats.csv"))


#compare lm models for group (with interactions) contribution to fraction of cells in sample (each cluster individually)

t1 = stat_tab
t1$dev_PCW_scaled = scale(t1$dev_PCW)

form1 = fract_sample ~ dev_PCW_scaled + 
  group + group:dev_PCW_scaled
form1_no_group = fract_sample ~ dev_PCW_scaled

lm_stats_by_cluster = NULL

for (cl in cluster_names){
  
  t2 = t1[t1$cluster_name == cl,]
  t3 = lm_anova(form1, form1_no_group, dat = t2, red = "group_w_int_removed")
  t4 = cbind(t3, cluster_name = cl)
  
  lm_stats_by_cluster = rbind(lm_stats_by_cluster, t4)
  
}
names(lm_stats_by_cluster)
lm_stats_by_cluster = lm_stats_by_cluster[!is.na(lm_stats_by_cluster$`Pr(>Chi)`),]
lm_stats_by_cluster$padj = p.adjust(lm_stats_by_cluster$`Pr(>Chi)`, method = "BH")

write_csv(lm_stats_by_cluster, 
          file = paste0(out_dir,script_ind,"cell_abundance_by_cluster_lm_LRT_stats.csv"))




###########################################################
# sccomp differential cell cluster abundance anlaysis
###########################################################

sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ group, 
    .sample =  sample, 
    .cell_group = cluster_name, 
    bimodal_mean_variability_association = TRUE,
    cores = 8 
  ) |> 
  #sccomp_remove_outliers(cores = 8) |> # Optional
  sccomp_test()

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group.csv"))


pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates.pdf"), 
    width = 4, height = 5)
sccomp_result |> 
  plot_1D_intervals()
dev.off()

pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot.pdf"), 
    width = 8, height = 8)
sccomp_result |> 
  sccomp_boxplot(factor = "group")
dev.off()


### analysis by group*dev_PCW

sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ group + dev_PCW + group:dev_PCW, 
    .sample =  sample, 
    .cell_group = cluster_name, 
    bimodal_mean_variability_association = TRUE,
    cores = 8 
  ) |> 
  #sccomp_remove_outliers(cores = 8) |> # Optional
  sccomp_test()

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_PCW.csv"))


pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_PCW_estimates.pdf"), 
    width = 8, height = 10)
sccomp_result |> 
  plot_1D_intervals()
dev.off()





#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed B05 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


