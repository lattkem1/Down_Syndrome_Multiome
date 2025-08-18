message("\n\n##########################################################################\n",
        "# Start G05a: Clustering tests and basic cell population analysis ", Sys.time(),
        "\n##########################################################################\n",
        "\n   clustering with selected cluster resolution ",
        "\n   differential abundance analysis on cluster level", 
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
library(GenomicRanges)
library(colorRamps)
library(viridis)
library(lmerTest)
library(pheatmap)
library(sccomp)

#specify script/output index as prefix for file names
script_ind = "G04_"

#specify output directory
out_dir = paste0(main_dir,"G_basic_analysis_grafts_map_to_tissue_v05_10X_Singl_DS2U_DS1/")


###load group and file info
gr_tab = read_csv(paste0(out_dir,"G02_gr_tab_filtered.csv"))

###load own dataset

load(file = paste0(out_dir,"G03_seur_with_ref_mapping.rda")) 



### get clusters of reference dataset

clust_tab_ref = read_csv("C_subsetting_exc_lin_from_all_non_cx_excl/C02_subset_cluster_assignment.csv")



### get marker gene panels
GOI = list()
t1 = read_csv("A_input/cell_type_markers_240806_consolidated.csv")
GOI$cell_type_markers = t1$gene[t1$level %in% c("cell_types", "neuronal_lineage")]
GOI$subtype_markers = t1$gene[t1$level %in% c("NPC_patterning","cortical_layers", "Interneuron_subtypes")]


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
# save cluster/cell class labels to dataset
###########################################################

DefaultAssay(seur) = "SCT"

seur$cluster_name = seur$predicted.cluster_name
seur$cell_type = clust_tab_ref$cell_type[match(seur$cluster_name, clust_tab_ref$cluster_name)]
seur$cell_class = clust_tab_ref$cell_class[match(seur$cluster_name, clust_tab_ref$cluster_name)]

save(seur, file = paste0(out_dir,script_ind, "seur_integr_labelled.rda")) 


##############################################################################################
#labelled dotplot for marker gene expression by predicted cluster (split by cell type markers vs neuron subtype markers)
##############################################################################################

DefaultAssay(seur) = "SCT"

cluster_names_ref = unique(clust_tab_ref$cluster_name)

p1 = DotPlot(seur, features = intersect(GOI$cell_type_markers, rownames(seur)), 
             group.by = "predicted.cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names_ref)

p2 = DotPlot(seur, features = intersect(GOI$subtype_markers, rownames(seur)), 
             group.by = "predicted.cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names_ref)

pdf(file = paste0(out_dir,script_ind, "Cell_markers_dotplot_by_predicted.clusters.pdf"), 
    width = 10, height = 8)
plot(p1)
plot(p2)
dev.off()



#############################################################################
# cluster quantification and differential abundance analysis
#############################################################################

meta = seur@meta.data
samples = unique(gr_tab$sample)
gr = unique(gr_tab$group)

# create table with one row for each cluster for each sample

stat_tab = tibble(cluster = unlist(lapply(cluster_names_ref, rep, length.out = length(samples))), 
                  sample = rep(samples, length(cluster_names_ref)))
stat_tab$group = gr_tab$group[match(stat_tab$sample, gr_tab$sample)] #add group assignment
stat_tab$isogenic_pair = gr_tab$isogenic_pair[match(stat_tab$sample, gr_tab$sample)]

#count cells per cluster per sample, add to stat_tab
t1 = meta %>% group_by(predicted.cluster_name, sample) %>% summarize(N_cells = n())
stat_tab$N_cells = t1$N_cells[match( paste0(stat_tab$cluster,stat_tab$sample), 
                                     paste0(t1$predicted.cluster_name,t1$sample) )]
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
  geom_point(data = t1, aes(x = cluster, y = fract_sample, color = group, shape = isogenic_pair), 
             position = position_dodge(width = 0.5), size = 0.5, stroke = 0.3)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names_ref)+
  scale_color_manual(limits = gr, values = pal(gr))+
  scale_fill_manual(limits = gr, values = pal(gr))+
  scale_shape_manual(limits = unique(stat_tab$isogenic_pair), 
                     values = rep_len(c(1:4), length.out = length(unique(stat_tab$isogenic_pair))))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_sample_pred_cluster.pdf"), 
    width = 7, height = 4)
plot(p1)
dev.off()


### statistical analysis t-tests by cluster

t1 = stat_tab

t2 = tibble(cluster = cluster_names_ref, p = NA, padj = NA)

for (i in cluster_names_ref){
  t4 = t1[t1$cluster == i,]
  t5 = pairwise.t.test(t4$fract_sample, t4$group, p.adjust.method = "none")
  t2$p[t2$cluster == i] = t5$p.value[1,1]
}

t2$padj = p.adjust(t2$p, method = "BH")

write_csv(t2, paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster_t-Test.csv"))




###########################################################
# sccomp differential cell cluster abundance anlaysis
###########################################################

sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ group, 
    .sample =  sample, 
    .cell_group = predicted.cluster_name, 
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




###########################################################
# sccomp differential cell cluster abundance anlaysis accounting for isogenic pair
###########################################################

sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ group + (1|seq_tech), 
    .sample =  sample, 
    .cell_group = predicted.cluster_name, 
    bimodal_mean_variability_association = TRUE,
    cores = 8 
  ) |> 
  #sccomp_remove_outliers(cores = 8) |> # Optional
  sccomp_test()

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_vs_covar.csv"))


pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_vs_covar.pdf"), 
    width = 8, height = 7)
sccomp_result |> 
  plot_1D_intervals()
dev.off()

pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_vs_covar.pdf"), 
    width = 8, height = 8)
sccomp_result |> 
  sccomp_boxplot(factor = "group")
dev.off()





#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed B05 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


