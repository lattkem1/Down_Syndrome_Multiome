message("\n\n##########################################################################\n",
        "# Start J03: DEG characterisation (reference mapped data)", Sys.time(),
        "\n##########################################################################\n",
        "\n   plot DEG stats per cluster, GO/MSigDB analysis per cluster up vs downreg genes",
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
script_ind = "J03_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_mapped_data_vs_reference_with_neurons_in_vitro_bulk/")


#load tissue bulk_data ( to extract GO terms and DEGs from tissue analysis)

load(file = paste0(main_dir,"E_DESeq_pseudobulk_all_by_cluster/E03_bulk_data_w_expr_z_scores.rda")) 
bulk_data_ref = bulk_data


#load combined foetal tissue/graft DEseq2 dataset

load(file = paste0(out_dir,"J02_bulk_data_with_DESeq_results.rda")) 




#get marker gene panels and network TFs

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
# Plot number of DEG (incl Chr21 genes) by cluster
#################################################

l1 = bulk_data$DEGs
l2 = lapply(l1, intersect, GOI$Chr21)
meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)
comps = unique(meta$comp)

t1 = NULL

t1 = expand_grid(comp = comps, DS_up_down = c("up","down"))
t1$cluster_name = meta$cluster_name[match(t1$comp, meta$comp)]
t1$sample_type = meta$sample_type[match(t1$comp, meta$comp)]
t1$N_genes = lengths(l1[paste0(t1$comp, "_DS_", t1$DS_up_down)])
t1$N_genes_Chr21 = lengths(l2[paste0(t1$comp, "_DS_", t1$DS_up_down)])

t1 = t1[order(match(t1$cluster_name, cluster_names)),]



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

pdf(file = paste0(out_dir,script_ind,"DEGs_numbers_by_cluster_sample_type.pdf"), 
    width = 4.5, height = 4)
{
  plot(p2)
  plot(p3)
  plot(p4)
}
dev.off()





#######################################
# calculate per gene Z-scores, mean Z-scores per cluster_name per group and group difference (DELTA), add to bulk_data
#######################################

message("\n\n          *** Calculate gene expression Z-scores... ", Sys.time(),"\n\n")

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



save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_w_expr_z_scores.rda"))





#########################################
#plot Chr21 DEGs
#########################################

message("\n\n          *** Plot heatmaps Chr21 genes... ", Sys.time(),"\n\n")

v1 = unlist(bulk_data_ref$DEGs)
Chr21_DEGs = intersect(GOI$Chr21, v1)

v1 = max(abs(na.omit(bulk_data$gene_Z_scores$mean_by_comp$DELTA)))
lims_DELTA = c(-0.15*v1, 0.15*v1)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_Chr21_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 10, height = 20)
{
  
  p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$CON, 
                        pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                        pl_genes = Chr21_DEGs,
                        x_col = "comp", 
                        meta_annot_cols = c("cluster_name", "sample_type"),
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = NULL,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, 
                        title = paste0("Chr21 DEG expression - Z-score mean CON"))
  
  pl_genes_clust = Chr21_DEGs[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$DELTA, 
                   pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "comp", 
                   meta_annot_cols = c("cluster_name", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                   lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, 
                   title = paste0("Chr21 DEG expression - Z-score mean DELTA"))
  
}

dev.off()  




#########################################
#plot network TFs
#########################################

message("\n\n          *** Plot heatmaps Netowrk TFs... ", Sys.time(),"\n\n")

v1 = rownames(bulk_data$gene_Z_scores$mean_by_comp$DELTA)
pl_genes = intersect(GOI$network_TFs, v1)

v1 = max(abs(na.omit(bulk_data$gene_Z_scores$mean_by_comp$DELTA)))
lims_DELTA = c(-0.15*v1, 0.15*v1)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_network_TFs_z_score_mean_CON_DELTA.pdf"), 
    width = 10, height = 20)
{
  
  p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$DELTA, 
                        pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                        pl_genes = GOI$network_TFs,
                        x_col = "comp", 
                        meta_annot_cols = c("cluster_name", "sample_type"),
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, 
                        title = paste0("Network TFs expression - Z-score mean DELTA"))
  
  pl_genes_clust = GOI$network_TFs[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$CON, 
                   pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "comp", 
                   meta_annot_cols = c("cluster_name", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = viridis(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, 
                   title = paste0("Network TFs expression - Z-score mean CON"))
  
  
  
  p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$DELTA, 
                        pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                        pl_genes = GOI$network_TFs_top20_Chr21,
                        x_col = "comp", 
                        meta_annot_cols = c("cluster_name", "sample_type"),
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, 
                        title = paste0("Network TFs expression - Z-score mean DELTA"))
  
  pl_genes_clust = GOI$network_TFs_top20_Chr21[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$CON, 
                   pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "comp", 
                   meta_annot_cols = c("cluster_name", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = viridis(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, 
                   title = paste0("Network TFs expression - Z-score mean CON"))
  
}

dev.off()  



#########################################
#plot genes by GO_term (max top20 simplified GO terms)
#########################################

message("\n\n          *** Plot heatmaps GO genes... ", Sys.time(),"\n\n")

go_res = bulk_data_ref$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = go_res$ID

v1 = max(abs(na.omit(bulk_data$gene_Z_scores$mean_by_comp$DELTA)))
lims_DELTA = c(-0.15*v1, 0.15*v1)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 10, height = 25)
{
  for (go in names(go_genes_list)){
    
    go_genes = go_genes_list[[go]]
    
    #plot clustered by DELTA expression Z-score (change in DS vs CON)
    
    p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$DELTA, 
                          pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                          pl_genes = go_genes,
                          x_col = "comp", 
                          meta_annot_cols = c("cluster_name", "sample_type"),
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                          lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, paste0(go_res$Description[go_res$ID == go]," (", go,
                                                ") - DEG expression - Z-score mean DELTA"))
    
    go_genes_clust = go_genes[p1$tree_row$order]
    
    bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$CON, 
                     pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                     pl_genes = go_genes_clust,
                     x_col = "comp", 
                     meta_annot_cols = c("cluster_name", "sample_type"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     color = viridis(250),
                     lims = NULL,  cellwidth = 10, cellheight = 10, 
                     fontsize = 10, title = paste0(go_res$Description[go_res$ID == go]," (", go,
                                                   ") - DEG expression - Z-score mean CON"))
    
  }
  
}
dev.off()





#########################################
#plot genes by HPO term (max top20 HPO terms)
#########################################

message("\n\n          *** Plot heatmaps HPO genes... ", Sys.time(),"\n\n")

t1 = bulk_data_ref$MSigDB_results
go_res = t1[t1$gs_subcat == "HPO",]
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = go_res$ID

v1 = max(abs(na.omit(bulk_data$gene_Z_scores$mean_by_comp$DELTA)))
lims_DELTA = c(-0.15*v1, 0.15*v1)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_HPO_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 15, height = 30)
{
  for (go in names(go_genes_list)){
    
    go_genes = go_genes_list[[go]]
    
    #plot clustered by DELTA expression Z-score (change in DS vs CON)
    
    p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$DELTA, 
                          pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                          pl_genes = go_genes,
                          x_col = "comp", 
                          meta_annot_cols = c("cluster_name", "sample_type"),
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                          lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, paste0(go_res$Description[go_res$ID == go]," (", go,
                                                ") - DEG expression - Z-score mean DELTA"))
    
    go_genes_clust = go_genes[p1$tree_row$order]
    
    bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_comp$CON, 
                     pl_meta = bulk_data$gene_Z_scores$mean_by_comp$meta,
                     pl_genes = go_genes_clust,
                     x_col = "comp", 
                     meta_annot_cols = c("cluster_name", "sample_type"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     color = viridis(250),
                     lims = NULL,  cellwidth = 10, cellheight = 10, 
                     fontsize = 10, title = paste0(go_res$Description[go_res$ID == go]," (", go,
                                                   ") - DEG expression - Z-score mean CON"))
    
  }
  
}
dev.off()





#get info on version of R, used packages etc
sessionInfo()





message("\n\n##########################################################################\n",
        "# Completed J03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


