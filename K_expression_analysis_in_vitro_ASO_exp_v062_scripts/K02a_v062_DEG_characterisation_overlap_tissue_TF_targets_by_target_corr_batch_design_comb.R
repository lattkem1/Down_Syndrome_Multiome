message("\n#################################\n",
        "####### Start K02 DESeq2 analysis: ", Sys.time(),
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


#specify script/output index as prefix for file names
script_ind = "K02a_"

#specify output directory
out_dir = paste0(main_dir,"K_expression_analysis_in_vitro_ASO_exp_v062/")


###load basal and ASO in vitro datasets (both with combined batches/designs), combine

load("J_expression_analysis_in_vitro_vs_exc_lin_v061_for_rev3/J01a_bulk_data_with_DESeq_results.rda")
bulk_data_bas = bulk_data

load(file = paste0(out_dir,"K01a_bulk_data_with_DESeq_results.rda")) 
bulk_data_ASO = bulk_data


#combine C9/C13 metadata
t1 = bulk_data_bas$meta
t2 = bulk_data_ASO$meta
t2 = t2[!t2$sample %in% t1$sample,]

t3 = full_join(t1, t2, by = intersect(names(t1), names(t2)))
t4 = t3[t3$isogenic_pair == "C9_C13" & t3$cell_type == "iNPC",]
rownames(t4) = t4$sample

meta_comb = t4

#combine count matrices
m1 = bulk_data_bas$counts
m2 = m1
colnames(m2) = t1$sample[match(colnames(m2), t1$cluster_sample)]

m3 = bulk_data_ASO$counts
m3 = m3[, !colnames(m3) %in% colnames(m2)]

keep_genes = intersect(rownames(m2), rownames(m3))

counts_comb = cbind(m2[keep_genes,], m3[keep_genes,])
counts_comb = counts_comb[,meta_comb$sample]

#combine DESeq2 results
l1 = bulk_data_bas$deseq_results
names(l1)
l1 = l1["iNPC_C9_C13"]
names(l1) = "DS_vs_CON"

l2 = bulk_data_ASO$deseq_results
names(l2)

deseq_results_comb = c(l1, l2)

#combine DEGs

l1 = bulk_data_bas$DEGs
names(l1)
l1 = l1[c("iNPC_C9_C13_up"  ,   "iNPC_C9_C13_down")]
names(l1) = c("DS_vs_CON_up", "DS_vs_CON_down")

l2 = bulk_data_ASO$DEGs

DEGs_comb = c(l1, l2)

bulk_data_comb = list(meta = meta_comb, gr_tab = meta_comb, counts = counts_comb,
                      deseq_results = deseq_results_comb, DEGs = DEGs_comb)


###load tissue pseudobulk dataset

load("E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_v045/E03_bulk_data_w_expr_z_scores.rda")
bulk_data_tissue = bulk_data


### get gene-of-interest sets

GOI = list()

t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol

t1 = read_tsv(paste0(main_dir,"A_input/Genomics_EnglandPanelApp_Intellectual_disability_v8.243.tsv"))
GOI$ID_genes = unique(t1$`Gene Symbol`)


#add TF targets

cand_TFs = c("PKNOX1", "BACH1", "GABPA")

t1 = read_csv(paste0(main_dir,"F_Chromatin_scMEGA_GRN_analysis_exc_lin_from_all_non_cx_excl_v045/F04_network_edges.csv"))
t1 = t1[t1$DS_reg,]
GOI$PKNOX1_targets = c("PNOX1", unique(t1$target[t1$tf == "PKNOX1"]))
GOI$BACH1_targets = c("BACH1", unique(t1$target[t1$tf == "BACH1"]))
GOI$GABPA_targets = c("GABPA", unique(t1$target[t1$tf == "GABPA"]))

net_edges_cand_TFs = t1[t1$tf %in% cand_TFs,]

net_edges_cand_TFs$Chr21_targets = net_edges_cand_TFs$target %in% GOI$Chr21


###define covariates to correct for

covars = c("batch")

###define comparisons
names(bulk_data_comb$deseq_results)
gr = unique(meta_comb$group)
gr

comps = list(DS_vs_CON = c("DS", "CON"),
             PKNOX1 = c("PKNOX1", "DS"),
             BACH1 = c("BACH1", "DS"),
             GABPA = c("GABPA", "DS")
             )


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
  } else if (v2 ==4){
    p2 = c("dodgerblue", "green4","grey20", "orange")
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




#########################################
#plot all DEGs combined
#########################################

bulk_data = bulk_data_comb

meta = bulk_data$meta


### extract vst norm expression matrix (based on all clusters combined), correct for covariates, calculate gene z-scores 

dds = DESeqDataSetFromMatrix(bulk_data$counts, colData = meta, 
                             design = ~group)

bulk_data$deseq_dataset_combined = dds

vst_mat = assay(vst(dds))

bulk_data$vst_mat_uncorr = vst_mat

#batch-correct vst matrix
vst_mat_corr = vst_mat

for (cov1 in covars){
  if (length(unique(meta[[cov1]]))>1){
    vst_mat_corr = limma::removeBatchEffect(vst_mat_corr, batch = meta[[cov1]])
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


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 




#######################################
# calculate mean gene Z-scores by by comparison, extract padj for each gene and comparison (combined dataset)
#######################################

message("\n\n          *** Calculate gene expression Z-scores... ", Sys.time(),"\n\n")

meta = bulk_data$meta

m1 = bulk_data$gene_Z_scores$clusters_combined

m_CON = matrix(nrow = nrow(m1), ncol = length(comps), 
               dimnames = list(rownames(m1), names(comps)))
m_X = m_CON
m_padj = m_CON
m_pval = m_CON

for (comp_name in names(comps)){
  
  comp_X = comps[[comp_name]][1]
  comp_CON = comps[[comp_name]][2]
  
  m_CON[,comp_name] = apply(m1[,meta$sample[meta$group == comp_CON] ], 1, mean)
  m_X[,comp_name] = apply(m1[,meta$sample[meta$group == comp_X] ], 1, mean)
  
  t1 = bulk_data$deseq_results[[comp_name]]
  m_padj[,comp_name] = t1[rownames(m_padj),][["padj"]]
  m_pval[,comp_name] = t1[rownames(m_pval),][["pvalue"]]
  
}

m_padj[is.na(m_padj)] = 1

meta_cl = tibble(comparison = names(comps))
rownames(meta_cl) = names(comps)

bulk_data$gene_Z_scores_cluster$meta = meta_cl
bulk_data$gene_Z_scores_cluster$CON = m_CON
bulk_data$gene_Z_scores_cluster$X = m_X
bulk_data$gene_Z_scores_cluster$DELTA = m_X - m_CON
bulk_data$gene_Z_scores_cluster$padj = m_padj 
bulk_data$gene_Z_scores_cluster$pval = m_pval


#########################################
#plot all DEGs combined
#########################################

meta = bulk_data$meta

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined

pl_genes = intersect(rownames(pl_mat_X ), 
                     unlist(bulk_data$DEGs))

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_combined.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "batch", "ASO_design", "target", "genotype"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 5, cellheight = 0.1, 
                        fontsize = 10, title = paste0("DEGs combined - Z-score"))
  
}
dev.off()




#########################################
#plot candidate TFs and TF targets
#########################################

pl_sets = list(cand_TFs = cand_TFs)

for (tf in cand_TFs){
  pl_sets[[tf]] = unique(c(tf, GOI[[paste0(tf, "_targets")]]))
}

meta = bulk_data$meta

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined

pl_genes = intersect(rownames(pl_mat_X ), 
                     unlist(pl_sets))

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_cand_TFs_pred_targets.pdf"), 
    width = 15, height = 30)
{
  
  for (pl_set in names(pl_sets)){
    
    pl_mat_X = bulk_data$gene_Z_scores$clusters_combined
    
    pl_genes = intersect(rownames(pl_mat_X ), 
                         pl_sets[[pl_set]])
    
    pl_mat_X = pl_mat_X[pl_genes,]
    
    p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                          pl_meta = meta,
                          pl_genes = pl_genes,
                          x_col = "sample", 
                          meta_annot_cols = c("group", "batch", "ASO_design", "target", "genotype"),
                          show_rownames = TRUE, show_colnames = TRUE,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = viridis(250),
                          lims = lims_X,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, title = paste0(pl_set, " - Expression Z-score"))
    
  }
}
dev.off()






#########################################
#analysis focussed on predicted TF targets (multiple testing correction for predicted targets only, use uncorrected p for trend): 
#   identify TF targets, targets regulated in DS vs CON, in TF KD vs DS, 
#   targets reverted by TF KD 
#########################################

t1 = bulk_data$deseq_results$PKNOX1
t2 = t1[t1$gene == "PKNOX1",]


pred_targets_deseq_list = NULL

tf = "PKNOX1"

###extract DESeq2 results for targets, adjust p-values

for (tf in cand_TFs){
  
  TF_targets = unique(c(tf, GOI[[paste0(tf, "_targets")]]))
  
  t1 = bulk_data$deseq_results$DS_vs_CON
  t2 = t1[t1$gene %in% TF_targets,]
  t2$padj = p.adjust(t2$pvalue, method = "BH")
  
  pred_targets_deseq_list[[paste0(tf, "_targets_DS_vs_CON")]] = t2
  
  t1 = bulk_data$deseq_results[[tf]]
  t2 = t1[t1$gene %in% TF_targets,]
  t2$padj = p.adjust(t2$pvalue, method = "BH")
  
  pred_targets_deseq_list[[paste0(tf, "_targets_DS_ASO_vs_DS")]] = t2
  
}


### extract genes with trend for opposing regulation in DS vs CON and ASO vs no ASO

pred_targets_reg_list = NULL

for (tf in cand_TFs){
  
  t1 = pred_targets_deseq_list[[paste0(tf, "_targets_DS_vs_CON")]]
  
  pred_targets_reg_list[[paste0(tf, "_targets_all")]] = t1$gene
  pred_targets_reg_list[[paste0(tf, "_targets_trend_DS_reg")]] = t1$gene[t1$pvalue<0.1]
  pred_targets_reg_list[[paste0(tf, "_targets_DS_reg")]] = t1$gene[t1$padj<0.1]
  
  #extract genes trend down in DS, up by ASO and vice versa
  
  t2 = pred_targets_deseq_list[[paste0(tf, "_targets_DS_ASO_vs_DS")]]
  
  v1 = t1$gene[t1$pval<0.1 & t1$log2FoldChange <0]
  v2 = t2$gene[t2$pval<0.1 & t2$log2FoldChange >0]
  pred_targets_reg_list[[paste0(tf, "_targets_trend_DS_down_ASO_up")]] = intersect(v1, v2)
  
  v1 = t1$gene[t1$padj<0.1 & t1$log2FoldChange <0]
  v2 = t2$gene[t2$padj<0.1 & t2$log2FoldChange >0]
  pred_targets_reg_list[[paste0(tf, "_targets_DS_down_ASO_up")]] = intersect(v1, v2)
  
  v1 = t1$gene[t1$pval<0.1 & t1$log2FoldChange >0]
  v2 = t2$gene[t2$pval<0.1 & t2$log2FoldChange <0]
  pred_targets_reg_list[[paste0(tf, "_targets_trend_DS_up_ASO_down")]] = intersect(v1, v2)
  
  v1 = t1$gene[t1$padj<0.1 & t1$log2FoldChange >0]
  v2 = t2$gene[t2$padj<0.1 & t2$log2FoldChange <0]
  pred_targets_reg_list[[paste0(tf, "_targets_DS_up_ASO_down")]] = intersect(v1, v2)
  
}


### save table with DEGs by comparison   (incl ID genes only)

l3 = pred_targets_reg_list

m1 = matrix(nrow = max(lengths(l3)), ncol = length(l3))
colnames(m1) = names(l3)

for (i in names(l3)){
  v1 = l3[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "TF_targets_up_down_DS_vs_CON_TF_KD_genes.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "TF_targets_up_down_DS_vs_CON_TF_KD_N_genes.csv"))


#save ID genes

l3 = lapply(pred_targets_reg_list, intersect, GOI$ID_genes)

m1 = matrix(nrow = max(lengths(l3)), ncol = length(l3))
colnames(m1) = names(l3)

for (i in names(l3)){
  v1 = l3[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "TF_targets_up_down_DS_vs_CON_TF_KD_ID_genes.csv"))



############################################################
# plot heatmaps candidate TFs and targets by sample
############################################################

pl_sets = c(list(cand_TFs = cand_TFs), pred_targets_reg_list)

meta = bulk_data$meta

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined

pl_genes = intersect(rownames(pl_mat_X ), 
                     unlist(pl_sets))

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_pred_targets_reg.pdf"), 
    width = 15, height = 30)
{
  
  for (pl_set in names(pl_sets)){
    
    pl_mat_X = bulk_data$gene_Z_scores$clusters_combined
    
    pl_genes = intersect(rownames(pl_mat_X ), 
                         pl_sets[[pl_set]])
    
    if (length(pl_genes)>1){
      
      pl_mat_X = pl_mat_X[pl_genes,]
      
      p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                            pl_meta = meta,
                            pl_genes = pl_genes,
                            x_col = "sample", 
                            meta_annot_cols = c("group", "batch", "ASO_design", "target", "genotype"),
                            show_rownames = TRUE, show_colnames = FALSE,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = viridis(250),
                            lims = lims_X,  cellwidth = 10, cellheight = 10, 
                            fontsize = 10, title = paste0(pl_set, " - Expression Z-score"))
    }
  }
}
dev.off()





#########################################
#plot group mean DELTA z-scores and p-value for each comparison for TF targets (genes on x axis)
#########################################

#define common z-score axis range
m1 = bulk_data$gene_Z_scores_cluster$DELTA

v1 = unique(unlist(pred_targets_reg_list))
v2 = max(abs(m1[intersect(rownames(m1), v1),]))

lims_Z = 0.5*c(-v2, v2)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_pred_targets_reg_DELTA.pdf"), 
    width = 30, height = 5)
{
  
  for (tf in cand_TFs){
    
    pl_sets = pred_targets_reg_list[grepl(tf, names(pred_targets_reg_list))]
    
    # adjust p values for all TF targets, create matrix with significance levels for plotting
    
    m1 = bulk_data$gene_Z_scores_cluster$pval
    TF_targets_all = pl_sets[[paste0(tf, "_targets_all")]]
    
    m2 = m1[TF_targets_all,]
    p_mat = m2
    m3 = apply(m2, 2, p.adjust, method = "BH")
    padj_mat = m3
    m4 = m2
    m4[m4>0.1] = ""
    m4[m2<=0.1] = "#"
    m4[m3<=0.1] = "*"
    
    sign_level_mat = m4
    
    #plot individual plot sets
    
    for (pl_set in names(pl_sets)){
      
      pl_genes = intersect(pl_sets[[pl_set]], rownames(sign_level_mat))
      
      if (length(pl_genes)>0){
        
        m1 = bulk_data$gene_Z_scores_cluster$DELTA
        
        m_DELTA = t(m1[pl_genes, ])
        
        m_sign = t(sign_level_mat[pl_genes,])
        
        p1 = pheatmap::pheatmap(m_DELTA, 
                                display_numbers = m_sign,
                                cluster_rows = FALSE, cluster_cols = TRUE,
                                show_rownames = TRUE, show_colnames = TRUE,
                                color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                                breaks = seq(lims_Z[1], lims_Z[2], length.out = 251),
                                border_color = NA, cellwidth = 10, cellheight = 10, 
                                fontsize = 10, number_color = "grey70", fontsize_number = 10, 
                                main = paste0(pl_set, " - DELTA Z-score"))
      }
    }
  }
}
dev.off()



#########################################
#plot group mean DELTA z-scores and p-value for each comparison for TF targets (genes on x axis)
#   only genes padj <0.1 highlighted 
#########################################

#define common z-score axis range
m1 = bulk_data$gene_Z_scores_cluster$DELTA

v1 = unique(unlist(pred_targets_reg_list))
v2 = max(abs(m1[intersect(rownames(m1), v1),]))

lims_Z = 0.5*c(-v2, v2)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_pred_targets_reg_DELTA_only_padj.pdf"), 
    width = 30, height = 5)
{
  
  for (tf in cand_TFs){
    
    pl_sets = pred_targets_reg_list[grepl(tf, names(pred_targets_reg_list))]
    
    # adjust p values for all TF targets, create matrix with significance levels for plotting
    
    m1 = bulk_data$gene_Z_scores_cluster$pval
    TF_targets_all = pl_sets[[paste0(tf, "_targets_all")]]
    
    m2 = m1[TF_targets_all,]
    p_mat = m2
    m3 = apply(m2, 2, p.adjust, method = "BH")
    padj_mat = m3
    m4 = m3
    m4[m3>0.1] = ""
    m4[m3<=0.1] = "*"
    
    sign_level_mat = m4
    
    #plot individual plot sets
    
    for (pl_set in names(pl_sets)){
      
      pl_genes = intersect(pl_sets[[pl_set]], rownames(sign_level_mat))
      
      if (length(pl_genes)>2){
        
        m1 = bulk_data$gene_Z_scores_cluster$DELTA
        
        m_DELTA = t(m1[pl_genes, ])
        
        m_sign = t(sign_level_mat[pl_genes,])
        
        p1 = pheatmap::pheatmap(m_DELTA, 
                                display_numbers = m_sign,
                                cluster_rows = FALSE, cluster_cols = TRUE,
                                show_rownames = TRUE, show_colnames = TRUE,
                                color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                                breaks = seq(lims_Z[1], lims_Z[2], length.out = 251),
                                border_color = NA, cellwidth = 10, cellheight = 10, 
                                fontsize = 10, number_color = "grey70", fontsize_number = 10, 
                                main = paste0(pl_set, " - DELTA Z-score"))
      }
    }
  }
}
dev.off()


#########################################
#plot group mean DELTA z-scores and p-value for each comparison for TF targets (genes on x axis)
#   only ASO for target gene set, padj <0.1 highlighted 
#########################################

#define common z-score axis range
m1 = bulk_data$gene_Z_scores_cluster$DELTA

v1 = unique(unlist(pred_targets_reg_list))
v2 = max(abs(m1[intersect(rownames(m1), v1),]))

lims_Z = 0.5*c(-v2, v2)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_pred_targets_reg_DELTA_only_target_ASO_comps.pdf"), 
    width = 30, height = 5)
{
  
  for (tf in cand_TFs){
    
    pl_sets = pred_targets_reg_list[grepl(tf, names(pred_targets_reg_list))]
    
    # adjust p values for all TF targets, create matrix with significance levels for plotting
    
    m1 = bulk_data$gene_Z_scores_cluster$pval
    TF_targets_all = pl_sets[[paste0(tf, "_targets_all")]]
    
    m2 = m1[TF_targets_all,]
    p_mat = m2
    m3 = apply(m2, 2, p.adjust, method = "BH")
    padj_mat = m3
    m4 = m2
    m4[m4>0.1] = ""
    m4[m2<=0.1] = "#"
    m4[m3<=0.1] = "*"
    
    sign_level_mat = m4
    
    #plot individual plot sets
    
    for (pl_set in names(pl_sets)){
      
      pl_genes = intersect(pl_sets[[pl_set]], rownames(sign_level_mat))
      
      if (length(pl_genes)>2){
        
        m1 = bulk_data$gene_Z_scores_cluster$DELTA
        
        m_DELTA = t(m1[pl_genes, c("DS_vs_CON", tf)])
        
        m_sign = t(sign_level_mat[pl_genes, c("DS_vs_CON", tf)])
        
        p1 = pheatmap::pheatmap(m_DELTA, 
                                display_numbers = m_sign,
                                cluster_rows = FALSE, cluster_cols = TRUE,
                                show_rownames = TRUE, show_colnames = TRUE,
                                color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                                breaks = seq(lims_Z[1], lims_Z[2], length.out = 251),
                                border_color = NA, cellwidth = 10, cellheight = 10, 
                                fontsize = 10, number_color = "grey70", fontsize_number = 10, 
                                main = paste0(pl_set, " - DELTA Z-score"))
      }
    }
  }
}
dev.off()





message("\n####### End K02: ", Sys.time(),"\n")


