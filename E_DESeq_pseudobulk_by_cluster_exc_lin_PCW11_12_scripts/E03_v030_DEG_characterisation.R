message("\n\n##########################################################################\n",
        "# Start E03: DEG characterisation ", Sys.time(),
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
script_ind = "E03_"

#specify output directory
out_dir = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_PCW11_12/")

#load group and file info
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")
gr_tab = gr_tab[gr_tab$dev_PCW <= 12,]

#load DEseq2 dataset, get genes for module clustering analysis

load(file = paste0(out_dir,"E02_bulk_data_with_DESeq_results.rda")) 


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

t1 = tibble(comps = names(l1), cluster_name = NA, DS_up_down = NA,
            N_genes = lengths(l1), N_genes_Chr21 = lengths(l2))

for (cl in cluster_names){
  t1$cluster_name[grepl(cl, t1$comps)] = cl
}

t1$DS_up_down[grepl("up", t1$comps)] = "up"
t1$DS_up_down[grepl("down", t1$comps)] = "down"

p1 = ggplot(t1, aes(x = cluster_name, group = DS_up_down))+
  scale_x_discrete(limits = cluster_names)+
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

pdf(file = paste0(out_dir,script_ind,"DEGs_numbers_by_cluster.pdf"), 
    width = 5, height = 4)
{
  plot(p2)
  plot(p3)
  plot(p4)
}
dev.off()




#################################################
### GO-BP over-representation analysis of DEGs combined => mapping to clusters
#################################################

message("\n\n          *** GO analysis DEGs vs clusters ... ", Sys.time(),"\n\n")

bulk_data$DEGs_combined = unique(unlist(bulk_data$DEGs))

### GO analysis

ego = enrichGO(gene         = bulk_data$DEGs_combined,
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
  
  pheatmap(top_go_mat, show_rownames=TRUE, show_colnames = TRUE,
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
  
  pheatmap(top_go_simpl_mat, show_rownames=TRUE, show_colnames = TRUE,
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


#save bulk dataset with GO analysis

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_w_expr_z_scores.rda"))






#################################################
# MSigDB over-representation analysis of DEGs combined => mapping to clusters
#################################################

message("\n\n          *** MSigDB overrepresentation analysis DEGs vs clusters... ", Sys.time(),"\n\n")

#get MSigDB collections of interest

msigdb_tab = msigdbr(species = "human")

msigdb_tab_sel = msigdb_tab[msigdb_tab$gs_cat %in% c("H", "C8")|
                              msigdb_tab$gs_subcat %in% c("CGP", "CP:REACTOME", "TFT:GTRD",
                                                          "GO:BP", "GO:CC", "GO:MF", "HPO"),]

msigdbr_t2g = msigdb_tab_sel %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()


#get DEG enrichment MSigDB collections

enrich_res = enricher(gene = bulk_data$DEGs_combined, TERM2GENE = msigdbr_t2g)
t2 = as.data.frame(enrich_res)
t2 = t2[t2$p.adjust <= 0.05,]
t2$gs_cat = msigdb_tab$gs_cat[match(t2$ID, msigdb_tab$gs_name)]
t2$gs_subcat = msigdb_tab$gs_subcat[match(t2$ID, msigdb_tab$gs_name)]
t2 = t2[order(t2$gs_cat, t2$gs_subcat),]

msigdb_results_tab = t2


### identify overlap of genes of msigdb terms with cluster DEGs

msigdb_results_tab_with_N_cluster_DEGs = msigdb_results_tab 
v1 = msigdb_results_tab$geneID
l1 = str_split(v1, "/")
names(l1) = msigdb_results_tab$ID

msigdb_gene_list = l1

for (comp in names(bulk_data$DEGs)){
  
  l2 = lapply(msigdb_gene_list, function(msigdb_genes){
    length(intersect(bulk_data$DEGs[[comp]], msigdb_genes))
  })
  msigdb_results_tab_with_N_cluster_DEGs[[comp]] = unlist(l2)
  
}



write_csv(msigdb_results_tab_with_N_cluster_DEGs, file = paste0(out_dir,script_ind, "MSigDB_DEGs_enrichment_results_combined.csv"))

bulk_data$MSigDB_results = msigdb_results_tab_with_N_cluster_DEGs

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_w_expr_z_scores.rda"))


###plot genes in top100/top30 terms/category vs clusters

t1 = msigdb_results_tab_with_N_cluster_DEGs
t1$gs_cat_comb = paste0(t1$gs_cat, "_", t1$gs_subcat)

pdf(file = paste0(out_dir,script_ind, "MSigDB_top100_terms_vs_clusters.pdf"), 
    width = 18, height = 30)
{
  for (cat in unique(t1$gs_cat_comb)){
    
    t3 = t1[t1$gs_cat_comb == cat,]
    if (nrow(t3)>100){t3 = t3[1:100,]}
    
    t4 = t3[,names(bulk_data$DEGs)]
    m1 = as.matrix(t4)
    rownames(m1) = t3$Description
    
    pheatmap(m1, show_rownames=TRUE, show_colnames = TRUE,
             cluster_rows = TRUE, cluster_cols = FALSE,  
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D2",
             treeheight_row = 10, treeheight_col = 10,
             color = colorRampPalette(c("white", "blue"))(250),
             breaks = seq(0, max(m1), length.out = 251),
             border_color = NA, fontsize = 10,
             cellwidth = 10, cellheight = 10,
             main = paste0(cat, " - Top100 MSigDB gene sets vs DEGs by cluster")
    )
    
    
    
  }
}

dev.off()


pdf(file = paste0(out_dir,script_ind, "MSigDB_top20_terms_vs_clusters.pdf"), 
    width = 18, height = 14)
{
  for (cat in unique(t1$gs_cat_comb)){
    
    t3 = t1[t1$gs_cat_comb == cat,]
    if (nrow(t3)>20){t3 = t3[1:20,]}
    
    t4 = t3[,names(bulk_data$DEGs)]
    m1 = as.matrix(t4)
    rownames(m1) = t3$Description
    
    pheatmap(m1, show_rownames=TRUE, show_colnames = TRUE,
             cluster_rows = TRUE, cluster_cols = FALSE,  
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D2",
             treeheight_row = 10, treeheight_col = 10,
             color = colorRampPalette(c("white", "blue"))(250),
             breaks = seq(0, max(m1), length.out = 251),
             border_color = NA, fontsize = 10,
             cellwidth = 10, cellheight = 10,
             main = paste0(cat, " - Top20 MSigDB gene sets vs DEGs by cluster")
    )
    
    
    
  }
}

dev.off()





#######################################
# calculate per gene Z-scores, mean Z-scores per cluster_stage/cluster_name per group and group difference (DELTA), add to bulk_data
#######################################

message("\n\n          *** Calculate gene expression Z-scores... ", Sys.time(),"\n\n")

meta = bulk_data$meta

vst_mat = assay(vst(bulk_data$deseq_dataset))

bulk_data$vst_mat = vst_mat


#calculate Z-score per gene by pseudobulk (cluster_sample)
cluster_sample_mat = t(apply(vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat)

bulk_data$gene_Z_scores$by_cluster_sample = cluster_sample_mat


#calculate mean Z-scores by cluster_name/cluster_stage

l1 = calculate_mean_Z_by(mat_full = bulk_data$gene_Z_scores$by_cluster_sample, 
                         meta = meta, 
                         mean_z_scores_by = "cluster_name")
bulk_data$gene_Z_scores$mean_by_cluster_name = l1

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_w_expr_z_scores.rda"))



#########################################
#plot DEGs by cluster_name
#########################################

pl_genes_list = bulk_data$DEGs

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_by_cluster_name_z_score_by_cluster_sample.pdf"), 
    width = 25, height = 60)
{
  for (cl in names(pl_genes_list)){
    
    pl_genes = pl_genes_list[[cl]]
    
    if(length(pl_genes)>1){
      
      bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$by_cluster_sample, 
                       pl_meta = bulk_data$meta,
                       pl_genes = pl_genes,
                       x_col = "cluster_sample", 
                       meta_annot_cols = c("group","cluster_name", "sample"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = viridis(250),
                       lims = NULL,  cellwidth = 3, cellheight = 10, 
                       fontsize = 10, title = paste0(cl, " - DEGs expression - Z-score by cluster_sample"))
      
      bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$by_cluster_sample, 
                       pl_meta = bulk_data$meta,
                       pl_genes = pl_genes,
                       x_col = "cluster_sample", 
                       meta_annot_cols = c("group","cluster_name", "sample"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                       lims = NULL,  cellwidth = 3, cellheight = 10, 
                       fontsize = 10, title = paste0(cl, " - DEGs expression - Z-score by cluster_sample"))
      
    }
    
  }
  
}
dev.off()





pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_by_cluster_name_z_score_mean_CON_DELTA.pdf"), 
    width = 12, height = 60)
{
  for (cl in names(pl_genes_list)){
    
    pl_genes = pl_genes_list[[cl]]
    
    if(length(pl_genes)>1){
      
      
      bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_cluster_name$CON, 
                       pl_meta = bulk_data$gene_Z_scores$mean_by_cluster_name$meta,
                       pl_genes = pl_genes,
                       x_col = "cluster_name", 
                       meta_annot_cols = c("cluster_name"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = viridis(250),
                       lims = NULL,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, title = paste0(cl, " - DEGs expression - Z-score mean CON"))
      
      bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_cluster_name$DELTA, 
                       pl_meta = bulk_data$gene_Z_scores$mean_by_cluster_name$meta,
                       pl_genes = pl_genes,
                       x_col = "cluster_name", 
                       meta_annot_cols = c("cluster_name"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                       lims = NULL,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, paste0(cl, " - DEGs expression - Z-score mean DELTA"))
    }
    
  }
  
}
dev.off()







#########################################
#plot genes by GO_term (max top20 GO terms)
#########################################

go_res = bulk_data$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = go_res$ID

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 12, height = 30)
{
  for (go in names(go_genes_list)){
    
    go_genes = go_genes_list[[go]]
    
    #plot clustered by DELTA expression Z-score (change in DS vs CON)
    
    p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_cluster_name$DELTA, 
                          pl_meta = bulk_data$gene_Z_scores$mean_by_cluster_name$meta,
                          pl_genes = go_genes,
                          x_col = "cluster_name", 
                          meta_annot_cols = c("cluster_name"),
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                          lims = NULL,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, paste0(go_res$Description[go_res$ID == go]," (", go,
                                                ") - DEG expression - Z-score mean DELTA"))
    
    go_genes_clust = go_genes[p1$tree_row$order]
    
    bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_cluster_name$CON, 
                     pl_meta = bulk_data$gene_Z_scores$mean_by_cluster_name$meta,
                     pl_genes = go_genes_clust,
                     x_col = "cluster_name", 
                     meta_annot_cols = c("cluster_name"),
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

t1 = bulk_data$MSigDB_results
go_res = t1[t1$gs_subcat == "HPO",]
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = go_res$ID

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_HPO_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 12, height = 30)
{
  for (go in names(go_genes_list)){
    
    go_genes = go_genes_list[[go]]
    
    #plot clustered by DELTA expression Z-score (change in DS vs CON)
    
    p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_cluster_name$DELTA, 
                          pl_meta = bulk_data$gene_Z_scores$mean_by_cluster_name$meta,
                          pl_genes = go_genes,
                          x_col = "cluster_name", 
                          meta_annot_cols = c("cluster_name"),
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                          lims = NULL,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, paste0(go_res$Description[go_res$ID == go]," (", go,
                                                ") - DEG expression - Z-score mean DELTA"))
    
    go_genes_clust = go_genes[p1$tree_row$order]
    
    bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$mean_by_cluster_name$CON, 
                     pl_meta = bulk_data$gene_Z_scores$mean_by_cluster_name$meta,
                     pl_genes = go_genes_clust,
                     x_col = "cluster_name", 
                     meta_annot_cols = c("cluster_name"),
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
        "# Completed E03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


