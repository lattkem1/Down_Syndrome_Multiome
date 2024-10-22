message("\n\n##########################################################################\n",
        "# Start J03b: DEG characterisation (reference (sc) vs bulk data)", Sys.time(),
        "\n##########################################################################\n",
        "\n   plot DEG stats, overlap bulk data vs reference",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(DESeq2)
library(colorRamps)
library(viridis)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(msigdbr)

#specify script/output index as prefix for file names
script_ind = "J03b_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_mapped_data_vs_reference_with_neurons_in_vitro_bulk/")


#load tissue bulk_data ( to extract GO terms and DEGs from tissue analysis)

load(file = paste0(main_dir,"E_DESeq_pseudobulk_all_by_cluster/E03_bulk_data_w_expr_z_scores.rda")) 
bulk_data_ref = bulk_data


#load combined foetal tissue/graft DEseq2 dataset

load(file = paste0(out_dir,"J02_bulk_data_with_DESeq_results.rda")) 
bulk_data_merged = bulk_data


#load in vitro neuron only DEseq2 dataset

load(file = paste0(out_dir,"J02b_bulk_data_with_DESeq_results.rda")) 
bulk_data_mapped = bulk_data

#get marker gene panels

GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol
t1 = read_csv(paste0(main_dir,"H_Pando_GRN_analysis_exc_lin_PCW11_12/H02_network_nodes.csv"))
GOI$network_TFs = t1$gene[t1$N_targets>0]
GOI$network_TFs_top20_Chr21 = na.omit(unique(c(GOI$network_TFs[1:20], 
                                               t1$gene[t1$N_targets>0 & t1$Chr21 == "Chr21"])))



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


### function: calculate and plot mean expr Z-score for selected metadata variable (mean_z_scores_by) 
#         for each CON and DS samples, and group difference

calculate_mean_Z_by = function(mat_full, meta, mean_z_scores_by){
  
  meta[["mean_var"]]= meta[[mean_z_scores_by]]
  meta2 = meta
  
  mean_var_vals = unique(meta$mean_var)
  
  
  mean_mat_CON = matrix(nrow= nrow(mat_full), ncol = length(mean_var_vals),
                        dimnames = list(rownames(mat_full), mean_var_vals))
  
  mean_mat_DS = mean_mat_CON 
  
  for (mean_var_val in mean_var_vals){
    
    samples_CON = meta$cluster_sample[meta$mean_var == mean_var_val & meta$group == "CON"]
    m1 = mat_full[,samples_CON]
    mean_mat_CON[,mean_var_val] = apply(m1, 1, mean)
    
    samples_DS = meta$cluster_sample[meta$mean_var == mean_var_val & meta$group == "DS"]
    m1 = mat_full[,samples_DS]
    mean_mat_DS[,mean_var_val] = apply(m1, 1, mean)
    
    #keep only meta columns with unique values for each value of the mean_var
    v1 = unlist(lapply(meta2[meta2$mean_var == mean_var_val,], function(x){length(unique(x))}))
    meta2 = meta2[, intersect(colnames(meta2), names(v1)[v1 == 1])]
    
  }
  meta2 = meta2[!duplicated(meta2$mean_var),]
  rownames(meta2) = meta2$mean_var
  
  mean_mat_list = list()
  mean_mat_list$meta = meta2
  mean_mat_list$CON = mean_mat_CON
  mean_mat_list$DS = mean_mat_DS
  mean_mat_list$DELTA = mean_mat_DS - mean_mat_CON
  
  return(mean_mat_list)
  
}



###function: plot bulk gene expression heatmap with annotations

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = NULL,  cellwidth = 15, cellheight = 10, 
                            fontsize = 10, title = "Z-score vst-norm gene expression"){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}
  
  pl_mat = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(pl_meta[[x_col]], colnames(pl_mat), nomatch = 0),]
  
  #define annotation bars
  
  if (!is.null(meta_annot_cols)){
    
    annot_col = data.frame(row.names = pl_meta[[x_col]]) #if not defined cbind converts factor values to factor levels
    
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      v1 = factor(v1, levels = unique(pl_meta[[col1]]))
      annot_col = as.data.frame(cbind(annot_col, v1))
    }
    colnames(annot_col) = meta_annot_cols
    rownames(annot_col) = colnames(pl_mat)
    
    annot_colors = lapply(meta_annot_cols, function(x){
      v1 = pal(unique(annot_col[[x]]) )
      names(v1) = levels(annot_col[[x]])
      return(v1)
    })
    names(annot_colors) = meta_annot_cols
    
  }else{
    annot_col = NULL
    annot_colors = NULL
  }
  
  # define plot limits
  
  if (is.null(lims)){lims = c(-0.7*max(abs(na.omit(pl_mat))), 0.7*max(abs(na.omit(pl_mat))))}
  
  #create plot
  
  p1 = pheatmap(pl_mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                color = color,
                breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                annotation_col = annot_col, annotation_colors = annot_colors,
                border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                fontsize = fontsize, main = title
  )
  
  return(p1)
  
}



#################################################
# Plot number of DEG (incl Chr21 genes and number overlapping with tissue) (only neurons in vitro)
#################################################

bulk_data = bulk_data_mapped

l1 = bulk_data$DEGs
l2 = lapply(l1, intersect, GOI$Chr21)
meta = bulk_data$meta

t1 = tibble(comp = "DS_vs_CON", DS_up_down = c("up","down"), 
               N_genes = lengths(l1), N_genes_Chr21 = lengths(l2))

for (i in 1:nrow(t1)){
  
  #genes up/down respectively in mapped cluster
  DEGs_cluster = bulk_data$DEGs[[paste0(t1$comp[i],"_", t1$DS_up_down[i])]]
  
  #all genes up or down respectively in any reference population 
  DEGs_ref = unlist(bulk_data_ref$DEGs[grepl(t1$DS_up_down[i], names(bulk_data_ref$DEGs))])
  
  t1$N_genes_overlap_any_foetal[i] = length(intersect(DEGs_cluster, DEGs_ref))
}



p1 = ggplot(t1, aes(x = comp, group = DS_up_down))+
  scale_x_discrete(limits = unique(t1$comp))+
  scale_fill_manual(limits = c("down", "up"), values = c("blue", "red"))+
  scale_color_manual(limits = c("down", "up"), values = c("blue", "red"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2 = p1 + geom_col(aes(y = N_genes, fill = DS_up_down), width = 0.7, position = position_dodge(width = 0.7))+
  labs(title = "All DEGs")

p3 = p1 + geom_col(aes(y = N_genes_Chr21, fill = DS_up_down), width = 0.7, position = position_dodge(width = 0.7))+
  labs(title = "Chr21 DEGs")

p4 = p1 + geom_col(aes(y = N_genes, color = DS_up_down), width = 0.7, position = position_dodge(width = 0.9),
                   fill = "grey80")+
  geom_col(aes(y = N_genes_Chr21, fill = DS_up_down, color = DS_up_down), 
           width = 0.7, position = position_dodge(width = 0.9))+
  labs(title = "DEG Chr21 vs others")

p5 = p1 + geom_col(aes(y = N_genes, color = DS_up_down), width = 0.7, position = position_dodge(width = 0.9),
                   fill = "grey80")+
  geom_col(aes(y = N_genes_overlap_any_foetal, fill = DS_up_down, color = DS_up_down), 
           width = 0.7, position = position_dodge(width = 0.9))+
  labs(title = "DEG with concordant reg in any foetal population vs others")


pdf(file = paste0(out_dir,script_ind,"DEGs_numbers_neurons_in_vitro.pdf"), 
    width = 3, height = 4)
{
  plot(p2)
  plot(p3)
  plot(p4)
  plot(p5)
}
dev.off()







#################################################
# Identify overlap all DEGs up/down in grafts vs tissue
#################################################

DEG_overlap_list = list()

l1 = bulk_data_ref$DEGs
l2 = bulk_data_mapped$DEGs

DEG_overlap_list[["foetal_up"]] = unique(unlist(l1[grepl("_up", names(l1))]))
DEG_overlap_list[["in_vitro_up"]] = unique(unlist(l2[grepl("_up", names(l2))]))
DEG_overlap_list[["foetal_up_in_vitro_up"]] = intersect(DEG_overlap_list[["foetal_up"]], 
                                                   DEG_overlap_list[["in_vitro_up"]])

DEG_overlap_list[["foetal_down"]] = unique(unlist(l1[grepl("_down", names(l1))]))
DEG_overlap_list[["in_vitro_down"]] = unique(unlist(l2[grepl("_down", names(l2))]))
DEG_overlap_list[["foetal_down_in_vitro_down"]] = intersect(DEG_overlap_list[["foetal_down"]], 
                                                       DEG_overlap_list[["in_vitro_down"]])

DEG_overlap_list[["foetal_up_in_vitro_down"]] = intersect(DEG_overlap_list[["foetal_up"]], 
                                                     DEG_overlap_list[["in_vitro_down"]])
DEG_overlap_list[["foetal_down_in_vitro_up"]] = intersect(DEG_overlap_list[["foetal_down"]], 
                                                     DEG_overlap_list[["in_vitro_up"]])


# save table with DEGs by comparison including TFs and HSA21 genes

l1 = DEG_overlap_list
l2 = lapply(l1, function(x){x = x[x%in% GOI$TF]})
names(l2) = paste0(names(l1), "_TF")
l3 = lapply(l1, function(x){x = x[x%in% GOI$Chr21]})
names(l3) = paste0(names(l1), "_Chr21")

l3 = c(l1, l2, l3)

m1 = matrix(nrow = max(lengths(l3)), ncol = length(l3))
colnames(m1) = names(l3)

for (i in names(l3)){
  v1 = l3[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "DEGs_overlaps_foetal_in_vitro_genes.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "DEGs_overlaps_foetal_in_vitro_N_genes.csv"))


# plot DEG numbers ref vs mapped

t2 = t1[1:8,]
t2$color = c("red", "orange", "red4", "blue", "dodgerblue", "blue4", "grey20", "grey30")

p1 = ggplot(t2, aes(x = gene_set, y = N_genes, fill = gene_set))+
  geom_col()+
  scale_x_discrete(limits = unique(t2$gene_set))+
  scale_fill_manual(limits = t2$gene_set, values = t2$color)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"DEGs_overlaps_foetal_in_vitro.pdf"), 
    width = 4, height = 3.5)
{
  plot(p1)
}
dev.off()





#######################################
# calculate for merged dataset per gene Z-scores, mean Z-scores per cluster_name per group and group difference (DELTA), add to bulk_data
#######################################

message("\n\n          *** Calculate gene expression Z-scores... ", Sys.time(),"\n\n")

bulk_data = bulk_data_merged

meta = bulk_data$meta

vst_mat = assay(vst(bulk_data$deseq_dataset))

bulk_data$vst_mat = vst_mat


#calculate Z-score per gene by pseudobulk (cluster_sample)
cluster_sample_mat = t(apply(vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat)

bulk_data$gene_Z_scores$by_cluster_sample = cluster_sample_mat


#calculate mean Z-scores by cluster_name/comp

l1 = calculate_mean_Z_by(mat_full = bulk_data$gene_Z_scores$by_cluster_sample, 
                         meta = meta, 
                         mean_z_scores_by = "cluster_name")
bulk_data$gene_Z_scores$mean_by_cluster_name = l1

l1 = calculate_mean_Z_by(mat_full = bulk_data$gene_Z_scores$by_cluster_sample, 
                         meta = meta, 
                         mean_z_scores_by = "comp")
bulk_data$gene_Z_scores$mean_by_comp = l1

bulk_data_merged = bulk_data

save(bulk_data_merged, file = paste0(out_dir, script_ind, "bulk_data_merged_w_expr_z_scores.rda"))



#########################################
#plot gene overlaps
#########################################

bulk_data = bulk_data_merged

v1 = max(abs(na.omit(bulk_data$gene_Z_scores$mean_by_comp$DELTA)))
lims_DELTA = c(-0.1*v1, 0.1*v1)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_z_score_CON_DELTA_foetal_in_vitro_overlap.pdf"), 
    width = 15, height = 30)
{
  
  p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$DELTA, 
                        pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                        pl_genes = DEG_overlap_list$foetal_up_in_vitro_up,
                        x_col = "comp", 
                        meta_annot_cols = c("cluster_name", "sample_type"),
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, 
                        title = paste0("Genes up tissue + in vitro expression - Z-score mean DELTA"))
  
  pl_genes_clust = DEG_overlap_list$foetal_up_in_vitro_up[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$CON, 
                   pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "comp", 
                   meta_annot_cols = c("cluster_name", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = viridis(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, 
                   title = paste0("Genes up tissue + in vitro expression - Z-score mean CON"))
  
  
  p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$DELTA, 
                        pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                        pl_genes = DEG_overlap_list$foetal_down_in_vitro_down,
                        x_col = "comp", 
                        meta_annot_cols = c("cluster_name", "sample_type"),
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, 
                        title = paste0("Genes down tissue + in vitro expression - Z-score mean DELTA"))
  
  pl_genes_clust =  DEG_overlap_list$foetal_down_in_vitro_down[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$CON, 
                   pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "comp", 
                   meta_annot_cols = c("cluster_name", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = viridis(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, 
                   title = paste0("Genes down tissue + in vitro expression - Z-score mean CON"))
  
  
  
}

dev.off()  




message("\n\n##########################################################################\n",
        "# Completed J03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


