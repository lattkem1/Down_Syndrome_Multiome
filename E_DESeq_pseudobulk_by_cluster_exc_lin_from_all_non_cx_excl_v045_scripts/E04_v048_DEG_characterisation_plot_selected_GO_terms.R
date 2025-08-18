message("\n\n##########################################################################\n",
        "# Start E04: plot selected GO terms ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
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
library(msigdbr)

#specify script/output index as prefix for file names
script_ind = "E04_"

#specify output directory
out_dir = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_v045/")


#load DEseq2 dataset
load(file = paste0(out_dir,"E03_bulk_data_w_expr_z_scores.rda")) 


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

# get table of selected GO terms

go_tab = read_csv(paste0(out_dir,"E03_GO_results_vs_N_DEGs_by_cluster_selected.csv"))



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





#################################################
### plot number of genes of GO term reg in each cluster
#################################################

t1 = go_tab

t2 = t1[,names(bulk_data$DEGs)]
m1 = as.matrix(t2)
rownames(m1) = paste0(t1$Description, " (", t1$ID, ")")

top_go_mat = m1


pdf(file = paste0(out_dir,script_ind, "Heatmap_selected_GO_terms_N_genes_vs_clusters.pdf"), 
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
           main = "Selected GO terms vs DEGs by cluster"
  )
  
}

dev.off()



#########################################
#plot genes by GO_term (selected GO terms) (selected clusters combined)
#########################################

go_res = go_tab

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sel_GO_genes_z_score_sel_clusters.pdf"), 
    width = 12, height = 30)
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
                          p_mat = FALSE,
                          meta_annot_cols = c("group", "dev_PCW", "cluster_name"),
                          show_rownames = TRUE, show_colnames = FALSE,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          viridis(250),
                          lims = lims_X,  cellwidth = 2, cellheight = 10, 
                          fontsize = 10, title = paste0(go, "- DEGs DS vs CON - Z-score"))
    
  }
  
}
dev.off()


#########################################
#plot genes by GO_term (selected GO terms) (mean DELTA CON)
#########################################

go_res = go_tab

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

#go = "forebrain development (GO:0030900)"     

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sel_GO_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 12, height = 20)
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
        "# Completed E03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


