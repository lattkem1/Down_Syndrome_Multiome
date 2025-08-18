message("\n\n##########################################################################\n",
        "# Start H03: DEG characterisation ", Sys.time(),
        "\n##########################################################################\n",
        "\n   plot DEG stats per cluster, GO/MSigDB analysis per cluster up vs downreg genes",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E14_241219_DS_foetal_brain_grafts_for_man_v02_low_string/")
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
library(ggrepel)

#specify script/output index as prefix for file names
script_ind = "H03_"

#specify output directory
out_dir = paste0(main_dir,"H_expression_analysis_grafts_vs_tissue_v05_10X_Singl_DS2U_DS1/")


#load DEseq2 dataset
load(file = paste0(out_dir,"H02_bulk_data_with_DESeq_results.rda")) 


#select clusters for plotting

cluster_names = unique(bulk_data$meta$cluster_name)
cluster_names
sel_clusters = c("RG_s1", "NEU_RORB_s4", "NEU_TLE4_s3" )

#get marker gene panels

GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol



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



###function: plot bulk gene expression heatmap with annotations

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, 
                            p_mat = FALSE,
                            pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = NULL,  cellwidth = 15, cellheight = 10, 
                            fontsize = 10, number_color = "grey70",
                            title = "Z-score vst-norm gene expression"){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}
  
  pl_mat = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(pl_meta[[x_col]], colnames(pl_mat), nomatch = 0),]
  
  #define annotation bars
  
  if (!is.null(meta_annot_cols)){
    
    annot_col = data.frame(row.names = pl_meta[[x_col]]) #if not defined cbind converts factor values to factor levels
    
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      if (is.numeric(v1)){v1 = factor(v1, levels = unique(pl_meta[[col1]][order(pl_meta[[col1]])]))} else {
        v1 = factor(v1, levels = unique(pl_meta[[col1]]))
      }
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
  
  p1 = pheatmap::pheatmap(pl_mat, 
                          display_numbers = p_mat,
                          cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                          show_rownames = show_rownames, show_colnames = show_colnames,
                          color = color,
                          breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                          annotation_col = annot_col, annotation_colors = annot_colors,
                          border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                          fontsize = fontsize, number_color = number_color, fontsize_number = fontsize, 
                          main = title
  )
  
  return(p1)
  
}




#######################################
# plot sample distance heatmap and PCA plot
#######################################

meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)
samples = unique(bulk_data$gr_tab$sample)
gr = unique(bulk_data$gr_tab$group)
cell_lines = unique(bulk_data$gr_tab$cell_line)
isogenic_pairs = unique(bulk_data$gr_tab$isogenic_pair)
seq_techs = unique(bulk_data$gr_tab$seq_tech)


### plot cluster_sample distance/correlation heatmap (by cluster)

pdf(file = paste0(out_dir,script_ind, "Heatmap_distance_matrix_cluster_samples.pdf"), 
    width = 12, height = 12)
{
  for (cl in cluster_names){
    
    #subset expression and metadata for cluster
    cl_meta = meta[meta$cluster_name == cl,]
    cl_samples = cl_meta$cluster_sample
    cl_vst_mat = bulk_data$vst_mat[,cl_samples]
    colnames(cl_vst_mat) = paste0(cl_meta$cluster_sample, "(", cl_meta$cell_line, ", ", cl_meta$seq_tech )
    
    #identify top 3000 variable genes
    v1 = rowVars(cl_vst_mat)
    var_genes = names(v1[order(-v1)][1:3000])
    
    #calculate z-scores
    z_mat = t(apply(cl_vst_mat[var_genes,], 1, scale))
    colnames(z_mat) = colnames(cl_vst_mat)
    
    #plot sample distance matrix 
    sampleDists <- dist(t(cl_vst_mat[var_genes,]))
    sampleDistMatrix <- as.matrix(sampleDists)
    
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=viridis(250),
             cellwidth = 10, cellheight = 10,
             main = paste0(cl, " - Cluster_sample Distance heatmap"))
    
    #plot sample distance matrix 
    
    cor_mat = cor(z_mat)
    
    pheatmap(cor_mat,
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(-1, 1, length.out = 251),
             cellwidth = 10, cellheight = 10,
             main = paste0(cl, " - Cluster_sample Correlation heatmap"))
    
  }
}
dev.off()


### plot PCA

pl = NULL # list of PCA plots

for (cl in cluster_names){
  
  #subset expression and metadata for cluster
  cl_meta = meta[meta$cluster_name == cl,]
  cl_samples = cl_meta$cluster_sample
  cl_vst_mat = bulk_data$vst_mat[,cl_samples]
  
  #identify top 3000 variable genes
  v1 = rowVars(cl_vst_mat)
  var_genes = names(v1[order(-v1)][1:3000])
  
  #calculate z-scores
  z_mat = t(apply(cl_vst_mat[var_genes,], 1, scale))
  colnames(z_mat) = colnames(cl_vst_mat)
  
  #calculate PC
  pc_analysis = prcomp(z_mat)
  
  meta_pca = cbind(cl_meta, pc_analysis$rotation)
  meta_pca$sample_label = paste0(meta_pca$sample, "(", meta_pca$cell_line, ")")
  
  pl[[cl]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group, shape = seq_tech))+
    geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
    scale_color_manual(limits = gr, values = pal(gr))+
    theme_minimal()+labs(title = paste0(cl, " - PCA"))
  
}
  

pdf(file = paste0(out_dir,script_ind, "PCA_plot_cluster_samples.pdf"), 
    width = 14, height = 12)
{
  lapply(pl, function(x){x})
}
dev.off()



#######################################
# calculate per gene Z-scores for combined dataset, selected clusters combined, and individual clusters
#######################################

message("\n\n          *** Calculate gene expression Z-scores for selected and individual clusters... ", Sys.time(),"\n\n")

#calculate Z-score per gene by pseudobulk (cluster_sample) for selected clusters

v1 = bulk_data$meta$cluster_sample[bulk_data$meta$cluster_name %in% sel_clusters]
m1 = bulk_data$vst_mat[,v1]

cluster_sample_mat = t(apply(m1, 1, scale))
colnames(cluster_sample_mat) = colnames(m1)

bulk_data$gene_Z_scores[["sel_clusters"]] = cluster_sample_mat


#calculate Z-score per gene by pseudobulk (cluster_sample) for all clusters individually

for (cl in cluster_names){
  
  v1 = bulk_data$meta$cluster_sample[bulk_data$meta$cluster_name == cl]
  m1 = bulk_data$vst_mat[, v1]
  
  cluster_sample_mat = t(apply(m1, 1, scale))
  colnames(cluster_sample_mat) = colnames(m1)
  
  bulk_data$gene_Z_scores[[cl]] = cluster_sample_mat
  
}



#######################################
# calculate mean gene Z-scores by cluster CON and DS vs CON, extract padj for each gene and comparison
#######################################

message("\n\n          *** Calculate gene expression Z-scores... ", Sys.time(),"\n\n")

meta = bulk_data$meta

m1 = bulk_data$gene_Z_scores$clusters_combined

m_CON = matrix(nrow = nrow(m1), ncol = length(cluster_names), 
               dimnames = list(rownames(m1), cluster_names))
m_DS = m_CON
m_padj = m_CON

for (cl in cluster_names){
  m_CON[,cl] = apply(m1[,meta$cluster_sample[meta$cluster_name == cl & meta$group == "CON"] ], 1, mean)
  m_DS[,cl] = apply(m1[,meta$cluster_sample[meta$cluster_name == cl & meta$group == "DS"] ], 1, mean)
  
  t1 = bulk_data$deseq_results[[cl]]
  m_padj[,cl] = t1[rownames(m_padj),][["padj"]]
  
}

m_padj[is.na(m_padj)] = 1

meta_cl = data.frame(cluster_name = cluster_names, 
                 cell_type = meta$cell_type[match(cluster_names, meta$cluster_name)])
rownames(meta_cl) = meta_cl$cluster_name

bulk_data$gene_Z_scores_cluster$meta = meta_cl
bulk_data$gene_Z_scores_cluster$CON = m_CON
bulk_data$gene_Z_scores_cluster$DS = m_DS
bulk_data$gene_Z_scores_cluster$DELTA = m_DS - m_CON
bulk_data$gene_Z_scores_cluster$padj = m_padj 




#########################################
#plot all DEGs combined by cluster_sample
#########################################

pl_genes = unique(unlist(bulk_data$DEGs))

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = meta,
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 4, cellheight = 0.2, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - Z-score"))
  
}
dev.off()


#########################################
#plot all DEGs combined by cluster_sample (selected clusters)
#########################################

pl_genes = unique(unlist(bulk_data$DEGs))

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 4, cellheight = 0.2, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - Z-score"))
  
}
dev.off()



#########################################
#plot all DEGs combined mean Z CON, DELTA (DS - CON)
#########################################

pl_genes = unique(unlist(bulk_data$DEGs))

pl_mat_CON = bulk_data$gene_Z_scores_cluster$CON[pl_genes,]
lims_CON = 0.7*c(-max(abs(pl_mat_CON)), max(abs(pl_mat_CON)))

pl_mat_DELTA = bulk_data$gene_Z_scores_cluster$DELTA[pl_genes,]
lims_DELTA = 0.4*c(-max(abs(pl_mat_DELTA)), max(abs(pl_mat_DELTA)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_DELTA, 
                        pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                        pl_genes = pl_genes,
                        x_col = "cluster_name", 
                        meta_annot_cols = c("cluster_name"),
                        show_rownames = FALSE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 0.2, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - mean DELTA Z-score"))
  
  pl_genes_clust = pl_genes[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = pl_mat_CON, 
                   pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE, show_colnames = TRUE,
                   color = viridis(250),
                   lims = lims_CON,  cellwidth = 10, cellheight = 0.2, 
                   fontsize = 10, title = paste0("DEGs DS vs CON combined - mean CON Z-score"))
  
}
dev.off()



#########################################
#plot DEGs by cluster_name
#########################################

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_by_cluster.pdf"), 
    width = 25, height = 60)
{
  for (cl in cluster_names){
    
    pl_genes = bulk_data$DEGs[[paste0(cl,"_up")]]
    pl_mat_X = bulk_data$gene_Z_scores[[cl]][pl_genes,]
    
    if(length(pl_genes)>1){
      
      bulkdata_heatmap(pl_mat = pl_mat_X, 
                       pl_meta = bulk_data$meta[bulk_data$meta$cluster_name == cl,],
                       pl_genes = pl_genes,
                       x_col = "cluster_sample", 
                       meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                       cluster_rows = TRUE, cluster_cols = FALSE,
                       color = viridis(250),
                       lims = NULL,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, title = paste0(cl, "_DS_up - DEGs expression - Z-score by cluster_sample"))
      
      
      pl_genes = bulk_data$DEGs[[paste0(cl,"_down")]]
      pl_mat_X = bulk_data$gene_Z_scores[[cl]][pl_genes,]
      
      if(length(pl_genes)>1){
        
        bulkdata_heatmap(pl_mat = pl_mat_X, 
                         pl_meta = bulk_data$meta[bulk_data$meta$cluster_name == cl,],
                         pl_genes = pl_genes,
                         x_col = "cluster_sample", 
                         meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                         cluster_rows = TRUE, cluster_cols = FALSE,
                         color = viridis(250),
                         lims = NULL,  cellwidth = 10, cellheight = 10, 
                         fontsize = 10, title = paste0(cl, "_DS_down - DEGs expression - Z-score by cluster_sample"))
        
      }
    }
  }
}
    
dev.off()


###smaller overview plot

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_by_cluster_small.pdf"), 
    width = 25, height = 20)
{
  for (cl in cluster_names){
    
    pl_genes = bulk_data$DEGs[[paste0(cl,"_up")]]
    pl_mat_X = bulk_data$gene_Z_scores[[cl]][pl_genes,]
    
    if(length(pl_genes)>1){
      
      bulkdata_heatmap(pl_mat = pl_mat_X, 
                       pl_meta = bulk_data$meta[bulk_data$meta$cluster_name == cl,],
                       pl_genes = pl_genes,
                       x_col = "cluster_sample", 
                       meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                       cluster_rows = TRUE, cluster_cols = FALSE,
                       show_rownames = FALSE, show_colnames = FALSE,
                       color = viridis(250),
                       lims = NULL,  cellwidth = 10, cellheight = 1, 
                       fontsize = 10, title = paste0(cl, "_DS_up - DEGs expression - Z-score"))
      
      
      pl_genes = bulk_data$DEGs[[paste0(cl,"_down")]]
      pl_mat_X = bulk_data$gene_Z_scores[[cl]][pl_genes,]
      
      if(length(pl_genes)>1){
        
        bulkdata_heatmap(pl_mat = pl_mat_X, 
                         pl_meta = bulk_data$meta[bulk_data$meta$cluster_name == cl,],
                         pl_genes = pl_genes,
                         x_col = "cluster_sample", 
                         meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                         cluster_rows = TRUE, cluster_cols = FALSE,
                         show_rownames = FALSE, show_colnames = FALSE,
                         color = viridis(250),
                         lims = NULL,  cellwidth = 10, cellheight = 1, 
                         fontsize = 10, title = paste0(cl, "_DS_down - DEGs expression - Z-score by cluster_sample"))
        
      }
    }
  }
}

dev.off()



#######################################
# plot diff expr TFs (sel clusters + DELTA)
#######################################

pl_genes = intersect(unique(unlist(bulk_data$DEGs)), GOI$TF)

### plot selected clusters by cluster_sample

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters_TFs.pdf"), 
    width = 15, height = 70)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                        show_rownames = TRUE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 4, cellheight = 5, 
                        fontsize = 5, title = paste0("DEGs DS vs CON combined - Z-score"))
  
}
dev.off()


### plot mean CON/DELTA Z-scores by cluster

pl_mat_CON = bulk_data$gene_Z_scores_cluster$CON[pl_genes,]
lims_CON = 0.7*c(-max(abs(pl_mat_CON)), max(abs(pl_mat_CON)))

pl_mat_DELTA = bulk_data$gene_Z_scores_cluster$DELTA[pl_genes,]
lims_DELTA = 0.4*c(-max(abs(pl_mat_DELTA)), max(abs(pl_mat_DELTA)))

pl_mat_padj = bulk_data$gene_Z_scores_cluster$padj[pl_genes,]
pl_mat_padj[pl_mat_padj<=0.1] = "*"
pl_mat_padj[pl_mat_padj!= "*"] = ""

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA_TFs.pdf"), 
    width = 15, height = 70)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_DELTA, 
                        pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                        pl_genes = pl_genes,
                        p_mat = pl_mat_padj,
                        x_col = "cluster_name", 
                        meta_annot_cols = c("cluster_name"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 5, 
                        fontsize = 5, title = paste0("DEGs DS vs CON combined - mean DELTA Z-score"))
  
  pl_genes_clust = pl_genes[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = pl_mat_CON, 
                   pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = TRUE, show_colnames = TRUE,
                   color = viridis(250),
                   lims = lims_CON,  cellwidth = 10, cellheight = 5, 
                   fontsize = 5, title = paste0("DEGs DS vs CON combined - mean CON Z-score"))
  
}
dev.off()




#######################################
# plot diff expr Chr21 genes
#######################################

pl_genes = intersect(unique(unlist(bulk_data$DEGs)), GOI$Chr21)

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters_Chr21_genes.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                        pl_genes = pl_genes,
                        p_mat = FALSE,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        viridis(250),
                        lims = lims_X,  cellwidth = 4, cellheight = 1, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - Z-score"))
  
}
dev.off()


### plot mean CON/DELTA Z-scores by cluster

pl_mat_CON = bulk_data$gene_Z_scores_cluster$CON[pl_genes,]
lims_CON = 0.7*c(-max(abs(pl_mat_CON)), max(abs(pl_mat_CON)))

pl_mat_DELTA = bulk_data$gene_Z_scores_cluster$DELTA[pl_genes,]
lims_DELTA = 0.4*c(-max(abs(pl_mat_DELTA)), max(abs(pl_mat_DELTA)))

pl_mat_padj = bulk_data$gene_Z_scores_cluster$padj[pl_genes,]
pl_mat_padj[pl_mat_padj<=0.1] = "*"
pl_mat_padj[pl_mat_padj!= "*"] = ""

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA_Chr21_genes.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_DELTA, 
                        pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                        pl_genes = pl_genes,
                        p_mat = pl_mat_padj,
                        x_col = "cluster_name", 
                        meta_annot_cols = c("cluster_name"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - mean DELTA Z-score"))
  
  pl_genes_clust = pl_genes[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = pl_mat_CON, 
                   pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                   pl_genes = pl_genes_clust,
                   p_mat = FALSE,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = TRUE, show_colnames = TRUE,
                   color = viridis(250),
                   lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, title = paste0("DEGs DS vs CON combined - mean CON Z-score"))
  
}
dev.off()






#################################################
### GO-BP over-representation analysis of DEGs combined => mapping to clusters
#################################################

message("\n\n          *** GO analysis DEGs vs clusters ... ", Sys.time(),"\n\n")

### GO analysis

ego = enrichGO(gene         = unique(unlist(bulk_data$DEGs)),
               OrgDb         = org.Hs.eg.db,
               keyType       = 'SYMBOL',
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)

GO_results_tab = ego@result[ego@result$p.adjust<=0.05,]
bulk_data$GO_results$full = GO_results_tab

# simpify GO results (usually helpful, sometimes removes also many intersting terms) 
ego_simplified =  dropGO(ego, level = c(1,2)) #drops vast majority of terms in some cases
ego_simplified = simplify(ego_simplified) #drops vast majority of terms in some cases

GO_results_simplified_tab = ego_simplified@result
bulk_data$GO_results$simplified = GO_results_simplified_tab



### identify overlap of GO genes with cluster DEGs

GO_results_tab_with_N_cluster_DEGs = GO_results_tab
v1 = GO_results_tab$geneID
l1 = str_split(v1, "/")
names(l1) = GO_results_tab$ID

GO_gene_list = l1

for (comp in names(bulk_data$DEGs)){
  
  l2 = lapply(GO_gene_list, function(go_genes){
    length(intersect(bulk_data$DEGs[[comp]], go_genes))
  })
  GO_results_tab_with_N_cluster_DEGs[[comp]] = unlist(l2)
  
}

bulk_data$GO_results$full = GO_results_tab_with_N_cluster_DEGs

write_csv(GO_results_tab_with_N_cluster_DEGs, paste0(out_dir, script_ind, "GO_results_vs_N_DEGs_by_cluster.csv"))


### identify overlap of genes of simplified GO terms with cluster DEGs

GO_results_simplified_tab_with_N_cluster_DEGs = GO_results_simplified_tab
v1 = GO_results_simplified_tab$geneID
l1 = str_split(v1, "/")
names(l1) = GO_results_simplified_tab$ID

GO_gene_list = l1

for (comp in names(bulk_data$DEGs)){
  
  l2 = lapply(GO_gene_list, function(go_genes){
    length(intersect(bulk_data$DEGs[[comp]], go_genes))
  })
  GO_results_simplified_tab_with_N_cluster_DEGs[[comp]] = unlist(l2)
  
}

bulk_data$GO_results$simplified = GO_results_simplified_tab_with_N_cluster_DEGs

write_csv(GO_results_simplified_tab_with_N_cluster_DEGs, paste0(out_dir, script_ind, "GO_results_simplified_vs_N_DEGs_by_cluster.csv"))


#save bulk dataset with GO analysis

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_w_expr_z_scores.rda"))



### plot number of genes of GO term reg in each cluster

t1 = GO_results_tab_with_N_cluster_DEGs

t2 = t1[,names(bulk_data$DEGs)]
m1 = as.matrix(t2)
rownames(m1) = paste0(t1$Description, " (", t1$ID, ")")

if(nrow(m1)>20){m1 = m1[1:20,]}

top_go_mat = m1


t1 = GO_results_simplified_tab_with_N_cluster_DEGs

t2 = t1[,names(bulk_data$DEGs)]
m1 = as.matrix(t2)
rownames(m1) = paste0(t1$Description, " (", t1$ID, ")")

if(nrow(m1)>20){m1 = m1[1:20,]}

top_go_simpl_mat = m1


pdf(file = paste0(out_dir,script_ind, "GO_terms_top20_vs_clusters.pdf"), 
    width = 18, height = 14)
{
  
  pheatmap::pheatmap(top_go_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(top_go_mat), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Top20 GO terms vs DEGs by cluster"
  )
  
  pheatmap::pheatmap(top_go_simpl_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(top_go_simpl_mat), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Top20 simplified GO terms vs DEGs by cluster"
  )
  
}

dev.off()



#########################################
#plot genes by GO_term (max top20 GO terms) (selected clusters combined)
#########################################

go_res = bulk_data$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_sel_clusters.pdf"), 
    width = 12, height = 40)
{
  for (go in names(go_genes_list)){
    
    pl_genes = go_genes_list[[go]]
    
    pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
    lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))
    
    #plot clustered by DELTA expression Z-score (change in DS vs CON)
    
    p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                          pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                          pl_genes = pl_genes,
                          x_col = "cluster_sample", 
                          meta_annot_cols = c("group", "weeks_post_injection", "seq_tech","cluster_name"),
                          show_rownames = TRUE, show_colnames = FALSE,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = viridis(250),
                          lims = lims_X,  cellwidth = 4, cellheight = 10, 
                          fontsize = 10, title = paste0(go, "- DEGs DS vs CON - Z-score"))
    
  }
  
}
dev.off()


#########################################
#plot genes by GO_term (max top20 GO terms) (mean DELTA CON)
#########################################

go_res = bulk_data$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

#go = "forebrain development (GO:0030900)"     

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 12, height = 40)
{
  for (go in names(go_genes_list)){
    
    pl_genes = go_genes_list[[go]]
    
    pl_mat_CON = bulk_data$gene_Z_scores_cluster$CON[pl_genes,]
    lims_CON = 0.7*c(-max(abs(pl_mat_CON)), max(abs(pl_mat_CON)))
    
    pl_mat_DELTA = bulk_data$gene_Z_scores_cluster$DELTA[pl_genes,]
    lims_DELTA = 0.4*c(-max(abs(pl_mat_DELTA)), max(abs(pl_mat_DELTA)))
    
    pl_mat_padj = bulk_data$gene_Z_scores_cluster$padj[pl_genes,]
    pl_mat_padj[pl_mat_padj<=0.1] = "*"
    pl_mat_padj[pl_mat_padj!= "*"] = ""
    
    p1 = bulkdata_heatmap(pl_mat = pl_mat_DELTA, 
                          pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                          pl_genes = pl_genes,
                          p_mat = pl_mat_padj,
                          x_col = "cluster_name", 
                          meta_annot_cols = c("cluster_name"),
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          show_rownames = TRUE, show_colnames = TRUE,
                          color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                          lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, number_color = "grey70",
                          title = paste0(go, " - DEGs - Expr DS vs CON (cluster Z-score)"))
    
    pl_genes_clust = pl_genes[p1$tree_row$order]
    
    bulkdata_heatmap(pl_mat = pl_mat_CON, 
                     pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                     pl_genes = pl_genes_clust,
                     x_col = "cluster_name", 
                     p_mat = FALSE,
                     meta_annot_cols = c("cluster_name"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE, show_colnames = TRUE,
                     color = viridis(250),
                     lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                     fontsize = 10, title = paste0(go, " - DEGs - Expr CON (cluster Z-score)"))
  }
}
dev.off()




#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed H03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


