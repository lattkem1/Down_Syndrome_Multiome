message("\n\n##########################################################################\n",
        "# Start B05: MiloR differential abundance analysis ", Sys.time(),
        "\n##########################################################################\n",
        "\n    ",
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
library(DESeq2)
library(miloR)

#specify script/output index as prefix for file names
script_ind = "C04_"

#specify output directory
out_dir = paste0(main_dir,"C_subsetting_all_cells_non_cx_excl/")

#load group and file info; keep only samples
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered_non_cx_excl.csv")

#load subset dataset
clust_tab = read_csv(paste0(out_dir,"C02_subset_cluster_assignment.csv"))
load(file = paste0(out_dir, "C03_seur_integr_labelled.rda"))


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
# prepare MiloR analysis
###########################################################


#convert to milo object

sce <- as.SingleCellExperiment(seur)
milo <- Milo(sce)

rm(sce)
gc()


### build neighbourhood graph, check neighbourhood sizes (ideal 50-100)

message("\n\n          *** Build neighbourhood graph... ", Sys.time(),"\n\n")

milo <- buildGraph(milo, k = 50, d = 30)
milo <- makeNhoods(milo, prop = 0.05, k = 50, d=30, refined = TRUE)

save(milo, file = paste0(out_dir,script_ind, "milo.rda"))

     
pdf(file = paste0(out_dir,script_ind, "milo_neighbourhood_sizes.pdf"), 
    width = 5, height = 5)
plotNhoodSizeHist(milo)
dev.off()


### count neighbourhood cells

message("\n\n          *** Count neighbourhood cells... ", Sys.time(),"\n\n")

#load(file = paste0(out_dir,script_ind, "milo.rda"))

milo <- countCells(milo, meta.data = data.frame(colData(milo)), samples="sample")


message("\n\n          *** Save milo object... ", Sys.time(),"\n\n")

save(milo, file = paste0(out_dir,script_ind, "milo.rda"))


###########################################################
# MiloR differential abundance analysis
###########################################################

message("\n\n          *** Differential Abundance analysis... ", Sys.time(),"\n\n")

milo_design <- data.frame(colData(milo))[,c("sample", "group")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$sample

## Reorder rownames to match columns of nhoodCounts(milo)
milo_design <- milo_design[colnames(nhoodCounts(milo)), ]


#store the distances between nearest neighbors in the Milo object

milo <- calcNhoodDistance(milo, d=30)


message("\n\n          *** Save milo object with abundance analysis... ", Sys.time(),"\n\n")


save(milo, file = paste0(out_dir,script_ind, "milo.rda"))



### differential abundance test


message("\n\n          *** Differential neighbourhoods test... ", Sys.time(),"\n\n")


da_results <- testNhoods(milo, design = ~ group, design.df = milo_design)


#check results
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

t1 = da_results[da_results$SpatialFDR <=0.05,]


message("\n\n          *** Save milo object with abundance analysis... ", Sys.time(),"\n\n")


save(milo, da_results, file = paste0(out_dir,script_ind, "milo.rda"))



###########################################################
# Extract membership of cells in differential neighbourhoods
###########################################################


message("\n\n          *** Extract Milo differential neighbourhood data... ", Sys.time(),"\n\n")

t1 = milo@nhoods

# #neighbourhoods of individual cell
# t3 = t1[7,]
# t4 = t3[t3>0]

t2 = da_results[da_results$SpatialFDR <=0.05,]

neighb_sign = t1[,as.numeric(rownames(t2))]

v1 = rowSums(neighb_sign)
v1[1:100]
cells_in_neighb_sign = names(v1[v1>0])

meta = as.data.frame(milo@colData)
meta$in_diff_neighbourhood = rownames(meta) %in% cells_in_neighb_sign

t3 = meta %>% group_by(cluster_name) %>% 
  summarise(N_cells_total = n())

t4 = meta[meta$in_diff_neighbourhood,] %>% group_by(cluster_name) %>% 
  summarise(N_cells = n())

t3$N_in_diff_neighbourhood = t4$N_cells[match(t4$cluster_name, t3$cluster_name)]
t3$fraction_in_diff_neighbourhood = t3$N_in_diff_neighbourhood/t3$N_cells_total

write_csv(t3, file = paste0(out_dir, script_ind, "fraction_cluster_in_diff_neighbourhood.csv"))



###########################################################
# Extract max(-logp10) and max(abs(log2FC)) for neigbourhood of each cell 
###########################################################

meta = as.data.frame(milo@colData)

m1 = milo@nhoods

# #neighbourhoods of individual cell
# t3 = t1[7,]
# t4 = t3[t3>0]

t2 = da_results
t2$neg_log10_padj = -log10(t2$SpatialFDR)

neg_log10_padj_mat = m1

all_cells = rownames(m1)
N_cells = length(all_cells)

time_start = Sys.time()

for (i in 1:N_cells){
  
  c1 = all_cells[i]
  
  #extract neg_log10_padj for all neighbourhoods of cell, keep highest sign for as metadata for plotting
  v1 = m1[c1,]*t2$neg_log10_padj
  
  meta$neg_log10_padj[rownames(meta) == c1] = max(v1)
  
  #extract log2FC for all differential neighbourhoods of cell, keep highest sign for as metadata for plotting
  v2 = m1[c1,]*t2$logFC
  v2[v1 < -log10(0.05)] = 0
  
  meta$log2FC[rownames(meta) == c1] = v2[abs(v2) == max(abs(v2))]
  
  #monitor progress
  if (i == 100){
    time100 = Sys.time()
    delta_time = round((difftime(time100, time_start, units = "mins")*(N_cells-100)/100),1)
    message("                 *Estimated minutes remaining ", delta_time )
  }
  if (i%%1000 == 0){
    timeX = Sys.time()
    delta_time = round((difftime(timeX, time_start, units = "mins")*(N_cells-i)/i),1)
    message("                 *Processed cell ", i, " of ", N_cells, 
            " ; estimated remaining minutes ", delta_time)
  }
}


t1 = milo@colData
t1$milo_neg_log10_padj = meta$neg_log10_padj
t1$milo_log2FC = meta$log2FC
milo@colData = t1


save(milo, file = paste0(out_dir,script_ind, "milo.rda"))




###########################################################
# add Milo neg_log10_padj, log2FC to Seurat
###########################################################

#load(file = paste0(out_dir,script_ind, "milo.rda"))

t1 = as.data.frame(milo@colData)

seur$milo_neg_log10_padj = t1$milo_neg_log10_padj

seur$milo_sign = seur$milo_neg_log10_padj >= -log10(0.05)
seur$milo_neg_log10_padj_sign = seur$milo_neg_log10_padj
seur$milo_neg_log10_padj_sign[seur$milo_neg_log10_padj_sign< -log10(0.05)] = 0

seur$milo_log2FC = t1$milo_log2FC

#rm(milo)
#gc()

###########################################################
# Visualise differential abundance (manual)
###########################################################


message("\n\n          *** Visualise Milo differential neighbourhood data... ", Sys.time(),"\n\n")


### plot whole dataset as UMAP with named clusters/cell classes

#define grouping variables
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
libraries = unique(gr_tab$library)
clusters = clust_tab$cluster
cluster_names = clust_tab$cluster_name
cell_types = unique(clust_tab$cell_type)
cell_classes = unique(clust_tab$cell_class)
stages = unique(seur$stage[order(seur$stage)])


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(42)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 

t1 = seur@meta.data
t2 = t1[t1$cell_type == "NEU_RORB",]


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

pl[["milo_neg_log10_padj"]] = FeaturePlot(seur, cells = sample(Cells(seur)), 
                                          features = "milo_neg_log10_padj", order = FALSE,
                                          label = FALSE, reduction = "umap", pt.size = 0.01)+
  scale_color_viridis(option = "A")+
  labs(title = "UMAP by MiloR -log10(padj)")

pl[["milo_neg_log10_padj_sign"]] = FeaturePlot(seur, features = "milo_neg_log10_padj_sign", 
                                               cells = sample(Cells(seur)), 
                                               order = FALSE, label = FALSE, reduction = "umap", 
                                               pt.size = 0.01)+
  scale_color_viridis(option = "A")+
  labs(title = "UMAP by MiloR -log10(padj) only sign diff")

pl[["milo_log2FC"]] = FeaturePlot(seur, features = "milo_log2FC", cells = sample(Cells(seur)), 
                                  order = FALSE, label = FALSE, reduction = "umap", pt.size = 0.01)+
  scale_colour_gradient2(limits = c(-max(abs(seur$milo_log2FC)), max(abs(seur$milo_log2FC))), 
                         low = "blue", mid = "grey90", high = "red",
                         midpoint = 0, na.value = "grey90")+
  labs(title = "UMAP by MiloR log2FC")

pl[["cluster_name_umap"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pl[["cluster_name_umap_dist"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                         label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pdf(file = paste0(out_dir,script_ind, "UMAP_clusters_Milo_abundance.pdf"), width = 6, height = 5)
lapply(pl, function(x){x})
dev.off()




#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed B05 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


