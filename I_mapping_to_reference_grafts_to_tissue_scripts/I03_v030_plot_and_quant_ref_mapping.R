message("\n\n##########################################################################\n",
        "# Start I03 plot and quantify mapping of grafts to reference:  ", Sys.time(),
        "\n##########################################################################\n",
        "\n   \n",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
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
library(sccomp)

#specify script/output index as prefix for file names
script_ind = "I03_"

#specify output directory
out_dir = paste0(main_dir,"I_mapping_to_reference_grafts_to_tissue/")


#load reference dataset
load(file = paste0(main_dir,"B_basic_analysis/B04_seur_integr_labelled.rda") )
seur_ref = seur
rm(seur)
gc()

gr_tab_ref = read_csv(file =paste0(main_dir, "B_basic_analysis/B02_gr_tab_filtered.csv"))

clust_tab_ref = read_csv(paste0(main_dir, "B_basic_analysis/B03_cluster_assignment_all.csv"))


#load group and file info for data to map on reference

gr_tab_data = read_csv(paste0(main_dir, "A_input/group_tab_grafts.csv"))


#load own dataset

load(file = paste0(out_dir,"I02_seur_filtered_mapped.rda")) 



#############################
# functions
#############################

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




#############################################
# create table to plot reference vs query data
#############################################

#define grouping variables for plotting
gr = unique(gr_tab_ref$group)
samples = unique(gr_tab_ref$sample)
cluster_names = clust_tab_ref$cluster_name
cell_classes = unique(clust_tab_ref$cell_class)
cell_types = unique(clust_tab_ref$cell_type)
stages = unique(seur_ref$stage[order(seur_ref$stage)])

#subsample reference dataset for smaller plot sizes
if (nrow(seur_ref@meta.data)>50000){
  set.seed(123)
  seur_ref_plot = seur_ref[, sample(colnames(seur_ref), size =50000, replace=F)]
} else {seur_ref_plot = seur_ref}


#extract plot data from reference dataset + add graft-specific metadata columns (set sample as "ref")

t1 = FetchData(seur_ref_plot, vars = c("umap_1", "umap_2", "sample", "group", 
                                               "stage", "dev_PCW", "batch_seq", "cluster_name", 
                                               "cell_type", "cell_class"))
pl_tab = cbind(t1, cell_line = "tissue", protocol = "NA", weeks_post_injection = "NA")
pl_tab$sample_orig = pl_tab$sample
pl_tab$sample = "ref" 


#add plot data from graft datasets

t2 = FetchData(seur, vars = c("refUMAP_1", "refUMAP_2", "sample", "group", 
                              "predicted.stage", "predicted.dev_PCW","batch_seq", 
                              "predicted.cluster_name", 
                              "predicted.cell_type", "predicted.cell_class"))

t3 = FetchData(seur, vars = c("cell_line", "protocol", "weeks_post_injection"))
t4 = cbind(t2, t3, sample_orig = t2$sample)

names(t4) = names(pl_tab)

pl_tab = rbind(pl_tab, t4)

pl_tab$predicted_cluster_name = "ref"
pl_tab$predicted_cluster_name[pl_tab$cell_line != "tissue"] = pl_tab$cluster_name[pl_tab$cell_line != "tissue"]

pl_tab$group_grafts = "ref"
pl_tab$group_grafts[pl_tab$cell_line != "tissue"] = pl_tab$group[pl_tab$cell_line != "tissue"]


#############################################
# plot graft samples on reference UMAP
#############################################

#shuffle cells (except reference) to avoid unbalanced overplotting 
t1 = pl_tab[pl_tab$sample != "ref", ]
set.seed(123)
t1 = t1[sample(nrow(t1)),]
pl_tab_shuffled = rbind(pl_tab[pl_tab$sample == "ref",], t1)

#define grouping variables
samples = gr_tab_data$sample
cell_lines = unique(gr_tab_data$cell_line)
protocols = unique(gr_tab_data$protocol)
WPIs = unique(gr_tab_data$weeks_post_injection[order(gr_tab_data$weeks_post_injection)])


### plot graft cells by samples and covariates on reference UMAP

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_classic()

pl[["sample"]] = p1 + geom_point(aes(color = sample), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("ref", samples), values = c("grey", pal(samples)))+
  labs(title = "Samples mapped on reference")

pl[["group_grafts"]] = p1 + geom_point(aes(color = group_grafts), size = 2, alpha =0.7)+
  scale_color_manual(limits = c("ref", gr), values = c("grey", pal(gr)))+
  labs(title = "Group mapped on reference")

pl[["group_grafts_no_legend"]] = p1 + geom_point(aes(color = group_grafts), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("ref", gr), values = c("grey", pal(gr)))+
  labs(title = "Group mapped on reference")+NoLegend()

pl[["predicted_cluster_name"]] = p1 + geom_point(aes(color = predicted_cluster_name), size = 2, alpha =0.7)+
  scale_color_manual(limits = c("ref", cluster_names), values = c("grey", pal(cluster_names)))+
  labs(title = "Cluster_name mapped on reference")

pl[["predicted_cluster_name_no_legend"]] = p1 + geom_point(aes(color = predicted_cluster_name), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("ref", cluster_names), values = c("grey", pal(cluster_names)))+
  labs(title = "Cluster_name mapped on reference")+NoLegend()

pl[["predicted_cluster_name_dist"]] = p1 + geom_point(aes(color = predicted_cluster_name), size = 2, alpha =0.7)+
  scale_color_manual(limits = c("ref", cluster_names), values = c("grey", pal_dist(cluster_names)))+
  labs(title = "Cluster_name mapped on reference")

pl[["predicted_cluster_name_dist_no_legend"]] = p1 + geom_point(aes(color = predicted_cluster_name), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("ref", cluster_names), values = c("grey", pal_dist(cluster_names)))+
  labs(title = "Cluster_name mapped on reference")+NoLegend()

pl[["cell_line"]] = p1 + geom_point(aes(color = cell_line), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("tissue", cell_lines), values = c("grey", pal(cell_lines)))+
  labs(title = "Cell lines mapped on reference")

pl[["protocol"]] = p1 + geom_point(aes(color = protocol), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("NA", protocols), values = c("grey", pal(protocols)))+
  labs(title = "Protcols mapped on reference")

pl[["wpi"]] = p1 + geom_point(aes(color = weeks_post_injection), size = 0.5, alpha =0.7)+
  scale_color_manual(limits = c("NA", WPIs), values = c("grey", pal(WPIs)))+
  labs(title = "WPI mapped on reference")

pdf(file = paste0(out_dir,script_ind, "UMAP_by_sample_graft_covars.pdf"), 
    width = 4, height = 3.5)
lapply(pl, function(x){x} )
dev.off()



##################################################################################################
# plot cells by split by graft group, colored by predicted cluster_name, cell_class and stage 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_classic()

pl[["cluster_name"]] = p1 + geom_point(aes(color = cluster_name), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$group_grafts, levels = c("ref", gr)), nrow = 1)+
  labs(title = "cluster_name split by sample")
pl[["cell_class"]] = p1 + geom_point(aes(color = cell_class), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_classes, values = pal_dist(cell_classes))+
  facet_wrap(facets = factor(pl_tab_shuffled$group_grafts, levels = c("ref", gr)), nrow = 1)+
  labs(title = "cell_class split by sample")
pl[["cell_type"]] = p1 + geom_point(aes(color = cell_type), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_types, values = pal_dist(cell_types))+
  facet_wrap(facets = factor(pl_tab_shuffled$group_grafts, levels = c("ref", gr)), nrow = 1)+
  labs(title = "cell_type split by sample")
pl[["stage"]] = p1 + geom_point(aes(color = stage), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$group_grafts, levels = c("ref", gr)), nrow = 1)+
  labs(title = "stage split by sample")


pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_group.pdf"), 
    width = 3*length(samples)+4, height = 3)
lapply(pl, function(x){x} )
dev.off()



##################################################################################################
# plot cells by split by graft sample, colored by predicted cluster_name, cell_class and stage 
##################################################################################################

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_classic()

pl[["cluster_name"]] = p1 + geom_point(aes(color = cluster_name), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("ref", samples)), nrow = 1)+
  labs(title = "cluster_name split by sample")
pl[["cell_class"]] = p1 + geom_point(aes(color = cell_class), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_classes, values = pal_dist(cell_classes))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("ref", samples)), nrow = 1)+
  labs(title = "cell_class split by sample")
pl[["cell_type"]] = p1 + geom_point(aes(color = cell_type), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_types, values = pal_dist(cell_types))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("ref", samples)), nrow = 1)+
  labs(title = "cell_type split by sample")
pl[["stage"]] = p1 + geom_point(aes(color = stage), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$sample, levels = c("ref", samples)), nrow = 1)+
  labs(title = "stage split by sample")


pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_sample.pdf"), 
    width = 3*length(gr)+5, height = 3)
lapply(pl, function(x){x} )
dev.off()




### plot cells by split by cell line, colored by predicted cluster_name, cell_class and stage 

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_classic()

pl[["cluster_name"]] = p1 + geom_point(aes(color = cluster_name), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$cell_line, levels = c("tissue", cell_lines)), nrow = 1)+
  labs(title = "cluster_name split by cell line")
pl[["cell_class"]] = p1 + geom_point(aes(color = cell_class), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_classes, values = pal_dist(cell_classes))+
  facet_wrap(facets = factor(pl_tab_shuffled$cell_line, levels = c("tissue", cell_lines)), nrow = 1)+
  labs(title = "cell_class split by cell line")
pl[["cell_type"]] = p1 + geom_point(aes(color = cell_type), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_types, values = pal_dist(cell_types))+
  facet_wrap(facets = factor(pl_tab_shuffled$cell_line, levels = c("tissue", cell_lines)), nrow = 1)+
  labs(title = "cell_type split by cell line")
pl[["stage"]] = p1 + geom_point(aes(color = stage), size = 0.1)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$cell_line, levels = c("tissue", cell_lines)), nrow = 1)+
  labs(title = "stage split by cell line")


pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_cell_line.pdf"), 
    width = 3*length(cell_lines)+4, height = 3)
lapply(pl, function(x){x} )
dev.off()


### plot cells by split by protocol, colored by predicted cluster_name, cell_type and stage 

pl = list()

p1 = ggplot(pl_tab_shuffled, aes(x = umap_1, y = umap_2))+theme_classic()

pl[["cluster_name"]] = p1 + geom_point(aes(color = cluster_name), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  facet_wrap(facets = factor(pl_tab_shuffled$protocol, levels = c("NA", protocols)), nrow = 1)+
  labs(title = "cluster_name split by protocol")
pl[["cell_class"]] = p1 + geom_point(aes(color = cell_class), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_classes, values = pal_dist(cell_classes))+
  facet_wrap(facets = factor(pl_tab_shuffled$protocol, levels = c("NA", protocols)), nrow = 1)+
  labs(title = "cell_class split by protocol")
pl[["cell_type"]] = p1 + geom_point(aes(color = cell_type), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = cell_types, values = pal_dist(cell_types))+
  facet_wrap(facets = factor(pl_tab_shuffled$protocol, levels = c("NA", protocols)), nrow = 1)+
  labs(title = "cell_type split by protocol")
pl[["stage"]] = p1 + geom_point(aes(color = stage), size = 0.2, alpha = 0.7)+
  scale_color_manual(limits = stages, values = pal(stages))+
  facet_wrap(facets = factor(pl_tab_shuffled$protocol, levels = c("NA", protocols)), nrow = 1)+
  labs(title = "stage split by protocol")


pdf(file = paste0(out_dir,script_ind, "UMAP_split_by_protocol.pdf"), 
    width = 3*length(protocols)+4, height = 3)
lapply(pl, function(x){x} )
dev.off()



#################################################
#generate stats by cell cluster and sample
#################################################


#stats by cell abundance by cluster_name and sample

t1 = pl_tab[pl_tab$sample != "ref",]
t1$cluster_sample = paste0(t1$cluster_name, "_", t1$sample)

t2 = expand_grid(cluster_name = cluster_names, sample = samples)
t2$cluster_sample = paste0(t2$cluster_name, "_", t2$sample)
stat_by_cluster_sample = as_tibble(cbind(t2, t1[match(t2$sample, t1$sample),c("group", "weeks_post_injection", "cell_line", "protocol")]))

t3 = t1 %>% group_by(cluster_sample) %>% summarise(N_cells = n())

stat_by_cluster_sample$N_cells = t3$N_cells[match(stat_by_cluster_sample$cluster_sample, t3$cluster_sample)]

stat_by_cluster_sample$N_cells[is.na(stat_by_cluster_sample$N_cells)] = 0

N_sample = t1 %>% group_by(sample) %>% summarize(N_cells = n())
stat_by_cluster_sample$N_sample = N_sample$N_cells[match(stat_by_cluster_sample$sample, N_sample$sample)]
stat_by_cluster_sample$fract_sample = stat_by_cluster_sample$N_cells/stat_by_cluster_sample$N_sample


write_csv(stat_by_cluster_sample, file = paste0(out_dir, script_ind, "stats_cell_abundance_by_cluster_name_and_sample.csv"))





#barplot quantification of cell type by sample

t1 = stat_by_cluster_sample

p1 = ggplot()+
  geom_col(data = t1, aes(x = cluster_name, y = fract_sample, fill = sample), 
           position = position_dodge(preserve = "single"), width = 0.5, lwd = 0.4)+
  scale_color_manual(breaks = samples, values = pal(samples))+
  scale_fill_manual(breaks = samples, values = pal(samples))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind, "stats_cell_abundance_by_cluster_name_and_sample_barplot.pdf"), 
    width = 7, height = 3.5)
plot(p1)
dev.off()


### heatmap quantification of cell type and cluster_name by sample

heatmap_vir = function(m1, main = "", cluster_rows = FALSE){
  p1 = pheatmap(m1, show_rownames=TRUE, cluster_rows = cluster_rows,
                cluster_cols = FALSE, show_colnames = TRUE, 
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

t1 = stat_by_cluster_sample

m1 = matrix(nrow = length(cluster_names), ncol = length(samples), dimnames = list(cluster_names, samples))

for (i in samples){
  t2 = t1[t1$sample==i,]
  m1[,i] = t2$fract_sample[match(rownames(m1),t2$cluster_name)]
}
m1[is.na(m1)] = 0

pdf(file = paste0(out_dir,script_ind, "stats_cell_abundance_by_cluster_name_and_sample_heatmap.pdf"), 
    width = 6, height = 6)
heatmap_vir(m1, main = "Cell abundance by cluster_name and sample")
dev.off()


##################################
# plot stats by protocol/cell line
##################################


### crossbar-dotplot quantification of cluster contribution fraction of sample by protocol

t1 = stat_by_cluster_sample
t1$protocol = factor(t1$protocol, levels = protocols)
t2 = t1 %>% group_by(cluster_name, protocol) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster_name, y = mean_fract, color = protocol, group = protocol), fill = "grey90", 
           position = position_dodge(width = 0.7, preserve = "single"), width = 0.7, lwd = 0.3)+
  geom_errorbar(data = t2, aes(x = cluster_name,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               group = protocol),
                position = position_dodge(width = 0.7, preserve = "single"), width = 0.7, lwd = 0.2)+
  geom_point(data = t1, aes(x = cluster_name, y = fract_sample, 
                            fill = sample, color = protocol, group = protocol), 
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05, seed = 123), 
             size = 1, stroke = 0.01, shape = 21)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names)+
  scale_color_manual(limits = protocols, values = pal(protocols))+
  scale_fill_manual(limits = samples, values = pal(samples))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"stats_cell_abundance_by_cluster_name_and_protocol.pdf"), 
    width = 9, height = 3.5)
plot(p1)
dev.off()




### crossbar-dotplot quantification of cluster contribution fraction of sample by cell line

t1 = stat_by_cluster_sample
t1$cell_line = factor(t1$cell_line, levels = cell_lines)
t2 = t1 %>% group_by(cluster_name, cell_line) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster_name, y = mean_fract, color = cell_line, group = cell_line), fill = "grey90", 
           position = position_dodge(width = 0.7, preserve = "single"), width = 0.7, lwd = 0.3)+
  geom_errorbar(data = t2, aes(x = cluster_name,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               group = cell_line),
                position = position_dodge(width = 0.7, preserve = "single"), width = 0.7, lwd = 0.2)+
  geom_point(data = t1, aes(x = cluster_name, y = fract_sample, 
                            fill = sample, color = protocol, group = cell_line), 
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05, seed = 123), 
             size = 1, stroke = 0.01, shape = 21)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names)+
  scale_color_manual(limits = cell_lines, values = pal(cell_lines))+
  scale_fill_manual(limits = samples, values = pal(samples))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"stats_cell_abundance_by_cluster_name_and_cell_line.pdf"), 
    width = 9, height = 3.5)
plot(p1)
dev.off()




##################################
# plot cell number stats by group, t-test for comparing groups by cluster
##################################


### crossbar-dotplot quantification of cluster contribution fraction of sample by protocol

t1 = stat_by_cluster_sample
t1$group = factor(t1$group, levels = gr)
t2 = t1 %>% group_by(cluster_name, group) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster_name, y = mean_fract, color = group, group = group), fill = "grey90", 
           position = position_dodge(width = 0.7, preserve = "single"), width = 0.7, lwd = 0.3)+
  geom_errorbar(data = t2, aes(x = cluster_name,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               group = group),
                position = position_dodge(width = 0.7, preserve = "single"), width = 0.7, lwd = 0.2)+
  geom_point(data = t1, aes(x = cluster_name, y = fract_sample, 
                            fill = sample, color = group, group = group), 
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05, seed = 123), 
             size = 1, stroke = 0.01, shape = 21)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names)+
  scale_color_manual(limits = gr, values = pal(gr))+
  scale_fill_manual(limits = samples, values = pal(samples))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"stats_cell_abundance_by_cluster_name_and_group.pdf"), 
    width = 6, height = 3.5)
plot(p1)
dev.off()




### statistical analysis t-tests by cluster

t1 = stat_by_cluster_sample

t2 = tibble(cluster_name = unique(t1$cluster_name), p = NA, padj = NA)

for (cl in unique(t1$cluster_name)){
  t4 = t1[t1$cluster_name == cl,]
  t5 = pairwise.t.test(t4$fract_sample, t4$group, p.adjust.method = "none")
  t2$p[t2$cluster_name == cl] = t5$p.value[1,1]
}

t2$padj = p.adjust(t2$p, method = "BH")

write_csv(t2, paste0(out_dir,script_ind,"cell_abundance_by_group_vs_cluster_t-Test.csv"))





###########################################################
# sccomp differential cell cluster abundance anlaysis
###########################################################

seur$cluster_name = seur$predicted.cluster_name

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




#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed I03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
