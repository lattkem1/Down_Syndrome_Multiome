message("\n#################################\n",
        "####### Start J02 DESeq2 analysis: ", Sys.time(),
        "\n##################################\n")

#set environment and main directory
# 
# Sys.setenv(RENV_PATHS_ROOT = "/nemo/lab/patanir/home/users/lattkem/R/renv/")
# renv::load(project = "/nemo/lab/patanir/home/users/lattkem/scRNA_P01_241108_R_4.4.0/")

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
script_ind = "J02_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_in_vitro_vs_exc_lin_v050/")

#load group and file info

gr_tab_in_vitro = read_csv("A_input/group_tab_in_vitro_bulk.csv")

gr_tab_tissue = read_csv("B_basic_analysis/B02_gr_tab_filtered_non_cx_excl.csv")

#load tissue cluster assignment

clust_tab = read_csv(paste0(main_dir,"C_subsetting_exc_lin_from_all_non_cx_excl/C02_subset_cluster_assignment.csv"))


###load in vitro dataset

load(file = paste0(out_dir,"J01_bulk_data_with_DESeq_results.rda")) 
bulk_data_in_vitro = bulk_data

###load tissue pseudobulk dataset

load("E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_v045/E03_bulk_data_w_expr_z_scores.rda")
bulk_data_tissue = bulk_data


### get gene-of-interest sets

GOI = list()

t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol

#add all DS up/down genes in tissue

l1 = bulk_data_tissue$DEGs

GOI$DS_up_any_cluster = unique(unlist(l1[grepl("_up",names(l1))]))
GOI$DS_down_any_cluster = unique(unlist(l1[grepl("_down",names(l1))]))

#add TF targets

cand_TFs = c("PKNOX1", "BACH1", "GABPA")

t1 = read_csv(paste0(main_dir,"F_Chromatin_scMEGA_GRN_analysis_exc_lin_from_all_non_cx_excl_v045/F04_network_edges.csv"))
t1 = t1[t1$DS_reg,]
GOI$PKNOX1_targets = c("PNOX1", unique(t1$target[t1$tf == "PKNOX1"]))
GOI$BACH1_targets = c("BACH1", unique(t1$target[t1$tf == "BACH1"]))
GOI$GABPA_targets = c("GABPA", unique(t1$target[t1$tf == "GABPA"]))

net_edges_cand_TFs = t1[t1$tf %in% cand_TFs,]

net_edges_cand_TFs$Chr21_targets = net_edges_cand_TFs$target %in% GOI$Chr21


### merge bulk_data counts and metadata 

t1 = bulk_data_in_vitro$meta
t1$sample_type = "in_vitro"
t1$cluster_name = paste0(t1$comp)
t1$cluster_sample = paste0(t1$cluster_name, "_", t1$sample)

t2 = bulk_data_tissue$meta
t2$sample_type = "fetal"
t2$genotype = t2$disease
t2$isogenic_pair = "fetal"
t2$seq_tech = "Multiome_10X"

t3 = full_join(t1, t2, by = c("cluster_name","cluster_sample", "sample", "group", 
                              "sample_type", "genotype", "cell_type", "isogenic_pair", "seq_tech"))
t3$group = t3$genotype

meta_merged = t3

# merge count matrices

m1 = bulk_data_in_vitro$counts
colnames(m1) = t1$cluster_sample
m2 = bulk_data_tissue$counts

genes_merged = intersect(rownames(m1), rownames(m2))
counts_merged = cbind(m1[genes_merged,], m2[genes_merged,])



#merge deseq results (adjust names to cluster_names)

names(bulk_data_tissue$deseq_results) = str_remove_all(names(bulk_data_tissue$deseq_results), "_DS_vs_CON")

deseq_results_merged = c(bulk_data_in_vitro$deseq_results, bulk_data_tissue$deseq_results)


bulk_data_merged = list(meta = meta_merged,
                        counts = counts_merged,
                        deseq_results = deseq_results_merged
                        )


### define covariates to correct for

covars = c("seq_tech")


### define selected clusters of merged dataset for plotting

unique(meta_merged$cluster_name)

sel_clusters = c("iNPC_C9_C13_batch_LILAC1", "iNPC_DS2U_DS1_batch_LILAC1", "iNPC_DS2U_DS1_batch_LILAC3",
                 "iNEU_C9_C13_batch_LILAC1"  , "iNEU_DS2U_DS1_batch_LILAC3", 
                 "RG_s1", "NEU_RORB_s4", "NEU_TLE4_s3")

bulk_data_merged$sel_clusters = sel_clusters



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
# summarise all DEGs in vitro vs tissue
#################################################

v1 = unique(unlist(bulk_data_tissue$DEGs))

v2 = unique(unlist(bulk_data_in_vitro$DEGs))

v3 = intersect(v1, v2)

t1 = tibble(N_DEGs_all_tissue = length(v1),
            N_DEGs_all_tissue_Chr21 = length(intersect(v1, GOI$Chr21)),
            N_DEGs_all_tissue_TF = length(intersect(v1, GOI$TF)),
            N_DEGs_all_in_vitro = length(v2),
            N_DEGs_all_in_vitro_Chr21 = length(intersect(v2, GOI$Chr21)),
            N_DEGs_all_in_vitro_TF = length(intersect(v2, GOI$TF)),
            N_DEGs_tissue_and_in_vitro = length(v3),
            N_DEGs_tissue_and_in_vitro_Chr21 = length(intersect(v3, GOI$Chr21)),
            N_DEGs_tissue_and_in_vitro_TF = length(intersect(v3, GOI$TF)))

write_csv(t1, file = paste0(out_dir,script_ind, "DEGs_tissue_vs_in_vitro_combined_N_genes.csv"))





###########################################################
# identify overlaps between tissue DEGs and in vitro DEGs
###########################################################

overlap_mat = matrix(nrow = length(bulk_data_tissue$DEGs), ncol = length(bulk_data_in_vitro$DEGs),
                     dimnames = list(names(bulk_data_tissue$DEGs), names(bulk_data_in_vitro$DEGs)))

overlap_mat_N = overlap_mat
overlap_mat_fract = overlap_mat

for (DEG_sets_tissue in names(bulk_data_tissue$DEGs)){
  
  for (DEG_sets_in_vitro in names(bulk_data_in_vitro$DEGs)){
    
    DEGs_tissue = bulk_data_tissue$DEGs[[DEG_sets_tissue]]
    DEGs_in_vitro = bulk_data_in_vitro$DEGs[[DEG_sets_in_vitro]]
    
    overlap_mat_N[DEG_sets_tissue, DEG_sets_in_vitro] = length(intersect(DEGs_tissue, DEGs_in_vitro))
    overlap_mat_fract[DEG_sets_tissue, DEG_sets_in_vitro] = length(intersect(DEGs_tissue, DEGs_in_vitro))/
      length(DEGs_tissue)

  }
}


###plot overlap stats

pdf(file = paste0(out_dir,script_ind, "Overlap_tissue_DEGs_by_cluster_vs_in_vitro.pdf"), 
    width = 10, height = 10)
{
  
  pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(overlap_mat_N), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Number of DEGs by cluster vs reg in vitro models")
  
  pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(overlap_mat_N), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Number of DEGs by cluster vs reg in vitro models")
  
  
  pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(overlap_mat_fract), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Fraction of DEGs by cluster vs reg in vitro models")
  
  pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(overlap_mat_fract), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Fraction of DEGs by cluster vs reg in vitro models")
  
}

dev.off()




#######################################
# normalise merged in vitro/tissue datasets for gene expression plots 
#######################################

bulk_data = bulk_data_merged

meta = bulk_data$meta
counts_merged = bulk_data$counts

dds = DESeqDataSetFromMatrix(counts_merged, colData = meta, 
                             design = ~group, ignoreRank = TRUE)

vst_mat = assay(vst(dds))
bulk_data$vst_mat_uncorr = vst_mat


#batch-correct vst matrix
vst_mat_corr = vst_mat

for (cov1 in covars){
  if (length(unique(meta[[cov1]]))>1){
    vst_mat_corr = limma::removeBatchEffect(vst_mat_corr, batch = meta[[cov1]], group = meta$group)
  } else {vst_mat_corr = vst_mat_corr}
}

bulk_data$vst_mat = vst_mat_corr


#calculate Z-score per gene by pseudobulk (cluster_sample) for all clusters combined (from uncorrected and corrected vst matrix)
cluster_sample_mat = t(apply(vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat)
bulk_data$gene_Z_scores_uncorr[["clusters_combined"]] = cluster_sample_mat

cluster_sample_mat = t(apply(vst_mat_corr, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat_corr)
bulk_data$gene_Z_scores[["clusters_combined"]] = cluster_sample_mat



### calculate Z-score by cluster_sample for selected clusters

meta1 = bulk_data$meta[bulk_data$meta$cluster_name %in% bulk_data$sel_clusters,]

vst1 = bulk_data$vst_mat[,meta1$cluster_sample]

z_mat = t(apply(vst1, 1, scale))
colnames(z_mat) = colnames(vst1)

bulk_data$gene_Z_scores[["sel_clusters"]] = z_mat



#######################################
# calculate mean gene Z-scores by cluster CON and DS vs CON, extract padj for each gene and comparison (combined dataset)
#######################################

message("\n\n          *** Calculate gene expression Z-scores... ", Sys.time(),"\n\n")

meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)

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
                     cell_type = meta$cell_type[match(cluster_names, meta$cluster_name)],
                     sample_type = meta$sample_type[match(cluster_names, meta$cluster_name)]
)
rownames(meta_cl) = meta_cl$cluster_name

bulk_data$gene_Z_scores_cluster$meta = meta_cl
bulk_data$gene_Z_scores_cluster$CON = m_CON
bulk_data$gene_Z_scores_cluster$DS = m_DS
bulk_data$gene_Z_scores_cluster$DELTA = m_DS - m_CON
bulk_data$gene_Z_scores_cluster$padj = m_padj 




###########################################################
# plot correlation gene expression between populations (vst_norm, top2000 variable genes)
###########################################################

meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)


### calculate vst cluster mean for top2000 variable genes, add to matrix

vst_mat = bulk_data$vst_mat

v1 = rowVars(vst_mat)
var_genes = names(v1[order(-v1)][1:2000])

vst1 = vst_mat[var_genes,]

vst_mean_mat = matrix(nrow = nrow(vst1), ncol = length(cluster_names),
                      dimnames = list(rownames(vst1), cluster_names ))

for (cl in cluster_names){
  
  cl_samples = meta$cluster_sample[meta$cluster_name == cl]
  
  vst_cl = apply(vst1[,cl_samples], 1, mean)
  
  vst_mean_mat[,cl] = vst_cl
  
}


### plot expression correlation heatmap each cluster vs each other cluster 

m1 = vst_mean_mat
colnames(m1) = paste(cluster_names, "(", meta$sample_type[match(cluster_names, meta$cluster_name)],")")

cor_mat = cor(m1)

write.csv(cor_mat, file = paste0(out_dir,script_ind, "Expr_cluster1_vs_cluster2_var_genes_vst_corr.csv"))

pdf(file = paste0(out_dir,script_ind, "Expr_cluster1_vs_cluster2_var_genes_vst_corr.pdf"), 
    width = 10, height = 10)
{
  
  pheatmap(cor_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = viridis(250),
           breaks = seq(0, 1, length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Corr expr top 2000 variable genes (vst norm)")
  
  pheatmap(cor_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 30, treeheight_col = 30,
           color = viridis(250),
           breaks = seq(0, 1, length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Corr expr top 2000 variable genes (vst norm, clustered)")
  
}

dev.off()


### plot expression in each cluster vs each other cluster 

expr_cl1_cl2_tab = expand_grid(cluster_name1 = cluster_names,
                               cluster_name2 = cluster_names, gene = var_genes, Expr1 = 0, Expr2 = 0)

for (cl in cluster_names){
  
  vst_cl = vst_mean_mat[,cl]
  
  t1 = expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name1 == cl, ]
  t1$Expr1 = vst_cl[match(t1$gene, names(vst_cl))]
  expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name1 == cl, ] = t1
  
  t1 = expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name2 == cl, ]
  t1$Expr2 = vst_cl[match(t1$gene, names(vst_cl))]
  expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name2 == cl, ] = t1
  
}

expr_cl1_cl2_tab$cluster_name1 = factor(expr_cl1_cl2_tab$cluster_name1, levels = cluster_names)
expr_cl1_cl2_tab$cluster_name2 = factor(expr_cl1_cl2_tab$cluster_name2, levels = cluster_names)


p1 = ggplot(expr_cl1_cl2_tab, aes(x = Expr1, y = Expr2))+
  geom_point(size = 0.1)+
  geom_smooth(method = "lm")+
  ggplot2::facet_grid(rows = vars(cluster_name2), 
                      cols = vars(cluster_name1), 
                      scales = "free")

pdf(file = paste0(out_dir,script_ind, "Expr_cluster1_vs_cluster2_var_genes_vst.pdf"), 
    width = 30, height = 30)
{
  plot(p1)
}
dev.off()




###########################################################
# plot correlation log2FC tissue DEGs between populations
###########################################################


meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)


### extract log2FC by cluster, add to matrix

DEGs_tissue = unique(unlist(bulk_data_tissue$DEGs))

logFC_mat = matrix(nrow = length(DEGs_tissue), ncol = length(cluster_names),
                   dimnames = list(DEGs_tissue, cluster_names ))

for (cl in cluster_names){
  t1 = bulk_data$deseq_results[[cl]]
  t2 = t1[DEGs_tissue,]
  logFC_mat[,cl] = t2$log2FoldChange
}

logFC_mat[is.na(logFC_mat)] = 0


### plot expression correlation heatmap each cluster vs each other cluster 

m1 = logFC_mat
colnames(m1) = paste(cluster_names, "(", meta$sample_type[match(cluster_names, meta$cluster_name)],")")

cor_mat = cor(m1)

write.csv(cor_mat, file = paste0(out_dir,script_ind, "Expr_cluster1_vs_cluster2_tissue_DEGs_log2FC_corr.csv"))

pdf(file = paste0(out_dir,script_ind, "Expr_cluster1_vs_cluster2_tissue_DEGs_log2FC_corr.pdf"), 
    width = 10, height = 10)
{
  
  pheatmap(cor_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = viridis(250),
           breaks = seq(0, 1, length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Corr log2FC (tissue DEGs)")
  
  pheatmap(cor_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 30, treeheight_col = 30,
           color = viridis(250),
           breaks = seq(0, 1, length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Corr log2FC (tissue DEGs, clustered)")
  
}

dev.off()


### plot expression corr in each cluster vs each other cluster 

expr_cl1_cl2_tab = expand_grid(cluster_name1 = cluster_names,
                               cluster_name2 = cluster_names, gene = DEGs_tissue, Expr1 = 0, Expr2 = 0)

for (cl in cluster_names){
  
  vst_cl = logFC_mat[,cl]
  
  t1 = expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name1 == cl, ]
  t1$Expr1 = vst_cl[match(t1$gene, names(vst_cl))]
  expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name1 == cl, ] = t1
  
  t1 = expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name2 == cl, ]
  t1$Expr2 = vst_cl[match(t1$gene, names(vst_cl))]
  expr_cl1_cl2_tab[expr_cl1_cl2_tab$cluster_name2 == cl, ] = t1
}

expr_cl1_cl2_tab$cluster_name1 = factor(expr_cl1_cl2_tab$cluster_name1, levels = cluster_names)
expr_cl1_cl2_tab$cluster_name2 = factor(expr_cl1_cl2_tab$cluster_name2, levels = cluster_names)


p1 = ggplot(expr_cl1_cl2_tab, aes(x = Expr1, y = Expr2))+
  geom_point(size = 0.1)+
  geom_smooth(method = "lm")+
  ggplot2::facet_grid(rows = vars(cluster_name2), 
                      cols = vars(cluster_name1), 
                      scales = "free")

pdf(file = paste0(out_dir,script_ind, "Expr_cluster1_vs_cluster2_tissue_DEGs_log2FC.pdf"), 
    width = 30, height = 30)
{
  plot(p1)
}
dev.off()





#########################################
#plot all tissue DEGs combined by cluster_sample
#########################################

meta = bulk_data$meta

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined

pl_genes = intersect(rownames(pl_mat_X ), 
                     unlist(bulk_data_tissue$DEGs))

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_tissue_DEGs_combined.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = meta,
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 1, cellheight = 0.5, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - Z-score"))
  
}
dev.off()




#########################################
#plot all DEGs combined by cluster_sample (selected clusters)
#########################################

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters

pl_genes = unique(unlist(bulk_data_tissue$DEGs))

pl_genes = intersect(rownames(pl_mat_X), pl_genes)
pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 4, cellheight = 0.5, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - Z-score"))
  
}
dev.off()



#########################################
#plot all DEGs combined mean Z CON, DELTA (DS - CON)
#########################################

pl_mat_X = bulk_data$gene_Z_scores_cluster$CON

pl_genes = unique(unlist(bulk_data_tissue$DEGs))

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_CON = bulk_data$gene_Z_scores_cluster$CON[pl_genes,]
lims_CON = 0.7*c(-max(abs(pl_mat_CON)), max(abs(pl_mat_CON)))

pl_mat_DELTA = bulk_data$gene_Z_scores_cluster$DELTA[pl_genes,]
lims_DELTA = 0.2*c(-max(abs(pl_mat_DELTA)), max(abs(pl_mat_DELTA)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_DELTA, 
                        pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                        pl_genes = pl_genes,
                        x_col = "cluster_name", 
                        meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                        show_rownames = FALSE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 0.5, 
                        fontsize = 10, title = paste0("DEGs DS vs CON combined - mean DELTA Z-score"))
  
  pl_genes_clust = pl_genes[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = pl_mat_CON, 
                   pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE, show_colnames = TRUE,
                   color = viridis(250),
                   lims = lims_CON,  cellwidth = 10, cellheight = 0.5, 
                   fontsize = 10, title = paste0("DEGs DS vs CON combined - mean CON Z-score"))
  
}
dev.off()




#######################################
# plot diff expr TFs (sel clusters + DELTA)
#######################################

pl_genes = intersect(unique(unlist(bulk_data_tissue$DEGs)), GOI$TF)

### plot selected clusters by cluster_sample

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters

pl_genes = intersect(pl_genes, rownames(pl_mat_X))
pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters_TFs.pdf"), 
    width = 15, height = 70)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                        show_rownames = TRUE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 4, cellheight = 5, 
                        fontsize = 5, title = paste0("DEGs TFs - DS vs CON combined - Z-score"))
  
}
dev.off()


### plot mean CON/DELTA Z-scores by cluster

pl_mat_CON = bulk_data$gene_Z_scores_cluster$CON[pl_genes,]
lims_CON = 0.7*c(-max(abs(pl_mat_CON)), max(abs(pl_mat_CON)))

pl_mat_DELTA = bulk_data$gene_Z_scores_cluster$DELTA[pl_genes,]
lims_DELTA = 0.2*c(-max(abs(pl_mat_DELTA)), max(abs(pl_mat_DELTA)))

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
                        meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = lims_DELTA,  cellwidth = 10, cellheight = 5, 
                        fontsize = 5, title = paste0("DEGs TFs - DS vs CON combined - mean DELTA Z-score"))
  
  pl_genes_clust = pl_genes[p1$tree_row$order]
  
  bulkdata_heatmap(pl_mat = pl_mat_CON, 
                   pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                   pl_genes = pl_genes_clust,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = TRUE, show_colnames = TRUE,
                   color = viridis(250),
                   lims = lims_CON,  cellwidth = 10, cellheight = 5, 
                   fontsize = 5, title = paste0("DEGs TFs - DS vs CON combined - mean CON Z-score"))
  
}
dev.off()




#######################################
# plot diff expr Chr21 genes
#######################################

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters

pl_genes = intersect(unique(unlist(bulk_data_tissue$DEGs)), GOI$Chr21)

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

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
                        meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
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
                        meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
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
                   meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = TRUE, show_colnames = TRUE,
                   color = viridis(250),
                   lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, title = paste0("DEGs DS vs CON combined - mean CON Z-score"))
  
}
dev.off()



#######################################
# plot Chr21 TFs
#######################################

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters

pl_genes = cand_TFs

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters_Chr21_TFs.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                        pl_genes = pl_genes,
                        p_mat = FALSE,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
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

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA_Chr21_TFs.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_DELTA, 
                        pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                        pl_genes = pl_genes,
                        p_mat = pl_mat_padj,
                        x_col = "cluster_name", 
                        meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
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
                   meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = TRUE, show_colnames = TRUE,
                   color = viridis(250),
                   lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, title = paste0("DEGs DS vs CON combined - mean CON Z-score"))
  
}
dev.off()



#########################################
#plot genes by GO_term (max top20 GO terms) (selected clusters combined)
#########################################

go_res = bulk_data_tissue$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_sel_clusters.pdf"), 
    width = 12, height = 40)
{
  for (go in names(go_genes_list)){

    pl_mat_X = bulk_data$gene_Z_scores$sel_clusters
    
    pl_genes = go_genes_list[[go]]
    
    pl_genes = intersect(rownames(pl_mat_X), pl_genes)
    
    pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
    lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))
    
    #plot clustered by DELTA expression Z-score (change in DS vs CON)
    
    p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                          pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                          pl_genes = pl_genes,
                          x_col = "cluster_sample", 
                          meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
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

go_res = bulk_data_tissue$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

#go = "forebrain development (GO:0030900)"     

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_mean_CON_DELTA.pdf"), 
    width = 12, height = 40)
{
  for (go in names(go_genes_list)){
    
    pl_mat_X = bulk_data$gene_Z_scores_cluster$CON
    
    pl_genes = go_genes_list[[go]]
    
    pl_genes = intersect(rownames(pl_mat_X), pl_genes)
    
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
                          meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
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
                     meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE, show_colnames = TRUE,
                     color = viridis(250),
                     lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                     fontsize = 10, title = paste0(go, " - DEGs - Expr CON (cluster Z-score)"))
  }
}
dev.off()




#######################################
# plot predicted TF targets
#######################################


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters_Chr21_TF_targets.pdf"), 
    width = 15, height = 30)
{
  
  for (tf in cand_TFs){
    
    pl_mat_X = bulk_data$gene_Z_scores$sel_clusters
    
    pl_genes = net_edges_cand_TFs$target[net_edges_cand_TFs$tf == tf]
    
    pl_genes = intersect(rownames(pl_mat_X), pl_genes)
    
    pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
    lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))
    
    p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                          pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                          pl_genes = pl_genes,
                          p_mat = FALSE,
                          x_col = "cluster_sample", 
                          meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                          show_rownames = TRUE, show_colnames = FALSE,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          viridis(250),
                          lims = lims_X,  cellwidth = 4, cellheight = 10, 
                          fontsize = 10, title = paste0(tf ," - targets - Expression Z-score"))
    
  }
  
}
dev.off()


### plot mean CON/DELTA Z-scores by cluster

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA_Chr21_TF_targets.pdf"), 
    width = 15, height = 30)
{
  
  for (tf in cand_TFs){
    
    pl_mat_X = bulk_data$gene_Z_scores_cluster$CON
    
    pl_genes = net_edges_cand_TFs$target[net_edges_cand_TFs$tf == tf]
    
    pl_genes = intersect(rownames(pl_mat_X), pl_genes)
    
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
                          meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                          show_rownames = TRUE, show_colnames = TRUE,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                          lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, title = paste0(tf ," - targets - mean DELTA Z-score"))
    
    pl_genes_clust = pl_genes[p1$tree_row$order]
    
    bulkdata_heatmap(pl_mat = pl_mat_CON, 
                     pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                     pl_genes = pl_genes_clust,
                     p_mat = FALSE,
                     x_col = "cluster_name", 
                     meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE, show_colnames = TRUE,
                     color = viridis(250),
                     lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                     fontsize = 10, title = paste0(tf ," - targets - mean CON Z-score"))
    
  }
  
}
dev.off()



#######################################
# filter for genes concordantly regulated in models and tissue in at least one comparison
#######################################

### extract concordantly regulated DEGs

degs_tissue_up = unique(unlist(bulk_data_tissue$DEGs[grepl("_up", names(bulk_data_tissue$DEGs))]))
degs_tissue_down = unique(unlist(bulk_data_tissue$DEGs[grepl("_down", names(bulk_data_tissue$DEGs))]))

degs_in_vitro_up = unique(unlist(bulk_data_in_vitro$DEGs[grepl("_up", names(bulk_data_in_vitro$DEGs))]))
degs_in_vitro_down = unique(unlist(bulk_data_in_vitro$DEGs[grepl("_down", names(bulk_data_in_vitro$DEGs))]))

degs_conc_up = intersect(degs_tissue_up, degs_in_vitro_up)
degs_conc_down = intersect(degs_tissue_down, degs_in_vitro_down)


### extract concordantly regulated GOIs

l1 = list(degs_in_vitro_any_up = degs_in_vitro_up, degs_in_vitro_any_down = degs_in_vitro_down,
          degs_conc_up = degs_conc_up, degs_conc_down = degs_conc_down)

l2 = lapply(GOI, intersect, degs_conc_up)
names(l2) = paste0(names(GOI), "_conc_up")

l3 = lapply(GOI, intersect, degs_conc_down)
names(l3) = paste0(names(GOI), "_conc_down")

l4 = c(l1, GOI, l2, l3)

m1 = matrix(nrow = max(lengths(l4)), ncol = length(l4))
colnames(m1) = names(l4)

for (i in names(l4)){
  v1 = l4[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "DEGs_conc_up_down_genes.csv"))

t1 = tibble(gene_set = names(l4), N_genes = lengths(l4))

write_csv(t1, file = paste0(out_dir,script_ind, "DEGs_conc_up_down_N_genes.csv"))



#########################################
#plot concordantly regulated genes by GO_term (max top20 GO terms) (selected clusters combined)
#########################################

go_res = bulk_data_tissue$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_sel_clusters_conc_reg.pdf"), 
    width = 12, height = 40)
{
  for (go in names(go_genes_list)){
    
    pl_genes = intersect(go_genes_list[[go]], c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
      pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
      lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))
      
      #plot clustered by DELTA expression Z-score (change in DS vs CON)
      
      p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                            pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                            pl_genes = pl_genes,
                            x_col = "cluster_sample", 
                            meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                            show_rownames = TRUE, show_colnames = FALSE,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = viridis(250),
                            lims = lims_X,  cellwidth = 4, cellheight = 10, 
                            fontsize = 10, title = paste0(go, "- DEGs DS vs CON - Z-score"))
      
    }
  }
}
dev.off()


#########################################
#plot concordantly regulated genes by GO_term (max top20 GO terms) (mean DELTA CON)
#########################################

go_res = bulk_data_tissue$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

#go = "forebrain development (GO:0030900)"     

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_mean_CON_DELTA_conc_reg.pdf"), 
    width = 12, height = 40)
{
  for (go in names(go_genes_list)){
    
    pl_genes = intersect(go_genes_list[[go]], c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
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
                            meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
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
                       meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = TRUE, show_colnames = TRUE,
                       color = viridis(250),
                       lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, title = paste0(go, " - DEGs - Expr CON (cluster Z-score)"))
    }
  }
}
dev.off()




#######################################
# plot concordantly regulated predicted TF targets
#######################################


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters_Chr21_TF_targets_conc_reg.pdf"), 
    width = 15, height = 30)
{
  
  for (tf in cand_TFs){
    
    pl_genes = intersect(net_edges_cand_TFs$target[net_edges_cand_TFs$tf == tf], 
                         c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
      pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
      lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))
      
      p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                            pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                            pl_genes = pl_genes,
                            p_mat = FALSE,
                            x_col = "cluster_sample", 
                            meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                            show_rownames = TRUE, show_colnames = FALSE,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            viridis(250),
                            lims = lims_X,  cellwidth = 4, cellheight = 10, 
                            fontsize = 10, title = paste0(tf ," - targets - Expression Z-score"))
      
    }
  }
}
dev.off()


### plot mean CON/DELTA Z-scores by cluster

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA_Chr21_TF_targets_conc_reg.pdf"), 
    width = 15, height = 30)
{
  
  for (tf in cand_TFs){
    
    pl_genes = intersect(net_edges_cand_TFs$target[net_edges_cand_TFs$tf == tf], 
                         c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
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
                            meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                            show_rownames = TRUE, show_colnames = TRUE,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                            fontsize = 10, title = paste0(tf ," - targets - mean DELTA Z-score"))
      
      pl_genes_clust = pl_genes[p1$tree_row$order]
      
      bulkdata_heatmap(pl_mat = pl_mat_CON, 
                       pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                       pl_genes = pl_genes_clust,
                       p_mat = FALSE,
                       x_col = "cluster_name", 
                       meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = TRUE, show_colnames = TRUE,
                       color = viridis(250),
                       lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, title = paste0(tf ," - targets - mean CON Z-score"))
      
      
    }
  }
}
dev.off()





#######################################
# filter for genes high confidence concordantly regulated in models vs tissue 
# same direction in at least one comparison, not in opposite direction in any comparison
#######################################

### extract concordantly regulated DEGs

degs_tissue_up = unique(unlist(bulk_data_tissue$DEGs[grepl("_up", names(bulk_data_tissue$DEGs))]))
degs_tissue_down = unique(unlist(bulk_data_tissue$DEGs[grepl("_down", names(bulk_data_tissue$DEGs))]))

degs_in_vitro_up = unique(unlist(bulk_data_in_vitro$DEGs[grepl("_up", names(bulk_data_in_vitro$DEGs))]))
degs_in_vitro_down = unique(unlist(bulk_data_in_vitro$DEGs[grepl("_down", names(bulk_data_in_vitro$DEGs))]))

degs_conc_up = intersect(degs_tissue_up[!(degs_tissue_up %in% degs_tissue_down)], 
                         degs_in_vitro_up[!(degs_in_vitro_up %in% degs_in_vitro_down)])
degs_conc_down = intersect(degs_tissue_down[!(degs_tissue_down %in% degs_tissue_up)], 
                         degs_in_vitro_down[!(degs_in_vitro_down %in% degs_in_vitro_up)])


### extract concordantly regulated GOIs

l1 = list(degs_in_vitro_any_up = degs_in_vitro_up, degs_in_vitro_any_down = degs_in_vitro_down,
          degs_conc_up = degs_conc_up, degs_conc_down = degs_conc_down)

l2 = lapply(GOI, intersect, degs_conc_up)
names(l2) = paste0(names(GOI), "_conc_up")

l3 = lapply(GOI, intersect, degs_conc_down)
names(l3) = paste0(names(GOI), "_conc_down")

l4 = c(l1, GOI, l2, l3)

m1 = matrix(nrow = max(lengths(l4)), ncol = length(l4))
colnames(m1) = names(l4)

for (i in names(l4)){
  v1 = l4[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "DEGs_conc_up_down_genes_high_conf.csv"))

t1 = tibble(gene_set = names(l4), N_genes = lengths(l4))

write_csv(t1, file = paste0(out_dir,script_ind, "DEGs_conc_up_down_N_genes_high_conf.csv"))



#########################################
#plot high confidence concordantly regulated genes by GO_term (max top20 GO terms) (selected clusters combined)
#########################################

go_res = bulk_data_tissue$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_sel_clusters_conc_reg_high_conf.pdf"), 
    width = 12, height = 30)
{
  for (go in names(go_genes_list)){
    
    pl_genes = intersect(go_genes_list[[go]], c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
      pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
      lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))
      
      #plot clustered by DELTA expression Z-score (change in DS vs CON)
      
      p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                            pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                            pl_genes = pl_genes,
                            x_col = "cluster_sample", 
                            meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                            show_rownames = TRUE, show_colnames = FALSE,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = viridis(250),
                            lims = lims_X,  cellwidth = 4, cellheight = 10, 
                            fontsize = 10, title = paste0(go, "- DEGs DS vs CON - Z-score"))
      
    }
  }
}
dev.off()


#########################################
#plot high confidence concordantly regulated genes by GO_term (max top20 GO terms) (mean DELTA CON)
#########################################

go_res = bulk_data_tissue$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")

#go = "forebrain development (GO:0030900)"     

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_genes_z_score_mean_CON_DELTA_conc_reg_high_conf.pdf"), 
    width = 12, height = 30)
{
  for (go in names(go_genes_list)){
    
    pl_genes = intersect(go_genes_list[[go]], c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
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
                            meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
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
                       meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = TRUE, show_colnames = TRUE,
                       color = viridis(250),
                       lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, title = paste0(go, " - DEGs - Expr CON (cluster Z-score)"))
    }
  }
}
dev.off()




#######################################
# plot high confidence concordantly regulated predicted TF targets (excluding Chr21 genes => concordant regulation most likely due to increased gene dosage)
#######################################


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_sel_clusters_Chr21_TF_targets_conc_reg_high_conf.pdf"), 
    width = 15, height = 20)
{
  
  for (tf in cand_TFs){
    
    pl_genes = intersect(net_edges_cand_TFs$target[net_edges_cand_TFs$tf == tf & 
                                                     !(net_edges_cand_TFs$target %in% GOI$Chr21) ], 
                         c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
      pl_mat_X = bulk_data$gene_Z_scores$sel_clusters[pl_genes,]
      lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))
      
      p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                            pl_meta = bulk_data$meta[bulk_data$meta$cluster_name %in% sel_clusters,],
                            pl_genes = pl_genes,
                            p_mat = FALSE,
                            x_col = "cluster_sample", 
                            meta_annot_cols = c("group", "cluster_name", "cell_type", "seq_tech", "sample_type"),
                            show_rownames = TRUE, show_colnames = FALSE,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            viridis(250),
                            lims = lims_X,  cellwidth = 4, cellheight = 10, 
                            fontsize = 10, title = paste0(tf ," - targets - Expression Z-score"))
      
    }
  }
}
dev.off()


### plot mean CON/DELTA Z-scores by cluster

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined_CON_DELTA_Chr21_TF_targets_conc_reg_high_conf.pdf"), 
    width = 15, height = 30)
{
  
  for (tf in cand_TFs){
    
    pl_genes = intersect(net_edges_cand_TFs$target[net_edges_cand_TFs$tf == tf & 
                                                     !(net_edges_cand_TFs$target %in% GOI$Chr21) ], 
                         c(degs_conc_up, degs_conc_down))
    
    if (length(pl_genes)>1){
      
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
                            meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                            show_rownames = TRUE, show_colnames = TRUE,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                            fontsize = 10, title = paste0(tf ," - targets - mean DELTA Z-score"))
      
      pl_genes_clust = pl_genes[p1$tree_row$order]
      
      bulkdata_heatmap(pl_mat = pl_mat_CON, 
                       pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                       pl_genes = pl_genes_clust,
                       p_mat = FALSE,
                       x_col = "cluster_name", 
                       meta_annot_cols = c("cluster_name", "cell_type", "sample_type"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = TRUE, show_colnames = TRUE,
                       color = viridis(250),
                       lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, title = paste0(tf ," - targets - mean CON Z-score"))
      
      
    }
  }
}
dev.off()







#############################################
# specific analyses cand TFs and targets
#############################################

#load(file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 
#deseq_res_list = bulk_data$deseq_results

#### extract DESEq2 results for cand TFs 

t1 = NULL

for (comp in names(bulk_data_in_vitro$deseq_results)){
  t2 = bulk_data_in_vitro$deseq_results[[comp]]
  t3 = t2[c("PKNOX1", "BACH1", "GABPA"),]
  t4 = as_tibble(cbind(comp = comp, gene = rownames(t3), t3))
  t1 = rbind(t1, t4)
}

write_csv(t1, file = paste0(out_dir,script_ind, "DESeq2_results_cand_TFs.csv"))




###########################################################
# identify overlaps between tissue TF targets and in vitro DEGs
###########################################################

### extract concordantly regulated TF targets

degs_tissue_up = unique(unlist(bulk_data_tissue$DEGs[grepl("_up", names(bulk_data_tissue$DEGs))]))
degs_tissue_down = unique(unlist(bulk_data_tissue$DEGs[grepl("_down", names(bulk_data_tissue$DEGs))]))

degs_in_vitro_up = unique(unlist(bulk_data_in_vitro$DEGs[grepl("_up", names(bulk_data_in_vitro$DEGs))]))
degs_in_vitro_down = unique(unlist(bulk_data_in_vitro$DEGs[grepl("_down", names(bulk_data_in_vitro$DEGs))]))

degs_conc_up = intersect(degs_tissue_up, degs_in_vitro_up)
degs_conc_down = intersect(degs_tissue_down, degs_in_vitro_down)

degs_conc_up_high_conf = intersect(degs_tissue_up[!(degs_tissue_up %in% degs_tissue_down)], 
                         degs_in_vitro_up[!(degs_in_vitro_up %in% degs_in_vitro_down)])
degs_conc_down_high_conf = intersect(degs_tissue_down[!(degs_tissue_down %in% degs_tissue_up)], 
                           degs_in_vitro_down[!(degs_in_vitro_down %in% degs_in_vitro_up)])


l1 = NULL

for (tf in cand_TFs){
  
  targets_fetal = GOI[[paste0(tf, "_targets")]]
  
  l1[[paste0(tf, "_targets_fetal")]] = targets_fetal
  l1[[paste0(tf, "_conc")]] = intersect(targets_fetal, c(degs_conc_up, degs_conc_down))
  l1[[paste0(tf, "_conc_high_conf")]] = intersect(targets_fetal[!(targets_fetal %in% GOI$Chr21)], 
                                                  c(degs_conc_up_high_conf, degs_conc_down_high_conf))
  
}

# save table with targets, including concordantly reg up/down genes 

v1 = unique(unlist(bulk_data_tissue$DEGs[grepl("_up", names(bulk_data_tissue$DEGs))]))
l2 = lapply(l1, intersect, v1)
names(l2) = paste0(names(l1), "_up_any_cluster")

v1 = unique(unlist(bulk_data_tissue$DEGs[grepl("_down", names(bulk_data_tissue$DEGs))]))
l3 = lapply(l1, intersect, v1)
names(l3) = paste0(names(l3), "_down_any_cluster")


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

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "TF_targets_genes_up_down_concord_fetal_in_vitro.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "TF_targets_N_genes_up_down_concord_fetal_in_vitro.csv"))










message("\n####### End J02: ", Sys.time(),"\n")


