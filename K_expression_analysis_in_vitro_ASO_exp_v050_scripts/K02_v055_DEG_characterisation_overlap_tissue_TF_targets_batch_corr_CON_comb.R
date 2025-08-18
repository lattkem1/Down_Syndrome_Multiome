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
script_ind = "K02_"

#specify output directory
out_dir = paste0(main_dir,"K_expression_analysis_in_vitro_ASO_exp_v050/")

#load group and file info

gr_tab = read_csv("A_input/group_tab_in_vitro_bulk_ASO_exp.csv")


###load in vitro dataset

load(file = paste0(out_dir,"K01_bulk_data_with_DESeq_results.rda")) 
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





#########################################
#plot all DEGs combined
#########################################

bulk_data = bulk_data_in_vitro

meta = bulk_data$meta

pl_mat_X = bulk_data$gene_Z_scores$all_samples

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
                        meta_annot_cols = c("group", "batch", "genotype"),
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

pl_mat_X = bulk_data$gene_Z_scores$all_samples

pl_genes = intersect(rownames(pl_mat_X ), 
                     unlist(pl_sets))

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_cand_TFs_pred_targets.pdf"), 
    width = 15, height = 30)
{
  
  for (pl_set in names(pl_sets)){
    
    pl_mat_X = bulk_data$gene_Z_scores$all_samples
    
    pl_genes = intersect(rownames(pl_mat_X ), 
                         pl_sets[[pl_set]])
    
    pl_mat_X = pl_mat_X[pl_genes,]
    
    p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                          pl_meta = meta,
                          pl_genes = pl_genes,
                          x_col = "sample", 
                          meta_annot_cols = c("group", "batch", "genotype"),
                          show_rownames = TRUE, show_colnames = TRUE,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          color = viridis(250),
                          lims = lims_X,  cellwidth = 10, cellheight = 10, 
                          fontsize = 10, title = paste0(pl_set, " - Expression Z-score"))
    
  }
}
dev.off()





#########################################
#analysis focussed on predicted TF targets (use uncorrected p for trend): 
#   identify TF targets, targets regulated in DS vs CON, in TF KD vs DS, 
#   targets reverted by TF KD 
#########################################

pred_targets_reg_list = NULL

tf = "PKNOX1"

for (tf in cand_TFs){
  
  TF_targets = GOI[[paste0(tf, "_targets")]]
  pred_targets_reg_list[[paste0(tf, "_pred_targets")]] = TF_targets
  
  
  ### identify targets reg in DS vs CON (correct padj for considering targets only)
  
  deseq_res_DS_vs_CON = bulk_data$deseq_results$DS_vs_CON[TF_targets,]
  deseq_res_DS_vs_CON_trend = deseq_res_DS_vs_CON[!is.na(deseq_res_DS_vs_CON$pvalue)&
                                                   deseq_res_DS_vs_CON$pvalue <0.1,]
  deseq_res_DS_vs_CON_diff = deseq_res_DS_vs_CON[!is.na(deseq_res_DS_vs_CON$padj)&
                                                   deseq_res_DS_vs_CON$padj <0.1,]
  
  TF_targets_trend_DS_up = intersect(TF_targets, 
                               deseq_res_DS_vs_CON_trend$gene[deseq_res_DS_vs_CON_trend$log2FoldChange>0])
  pred_targets_reg_list[[paste0(tf, "_targets_trend_DS_up")]] = TF_targets_trend_DS_up
  TF_targets_DS_up = intersect(TF_targets, 
                               deseq_res_DS_vs_CON_diff$gene[deseq_res_DS_vs_CON_diff$log2FoldChange>0])
  pred_targets_reg_list[[paste0(tf, "_targets_DS_up")]] = TF_targets_DS_up
  
  TF_targets_trend_DS_down = intersect(TF_targets, 
                                     deseq_res_DS_vs_CON_trend$gene[deseq_res_DS_vs_CON_trend$log2FoldChange<0])
  pred_targets_reg_list[[paste0(tf, "_targets_trend_DS_down")]] = TF_targets_trend_DS_down
  TF_targets_DS_down = intersect(TF_targets, 
                                 deseq_res_DS_vs_CON_diff$gene[deseq_res_DS_vs_CON_diff$log2FoldChange<0])
  pred_targets_reg_list[[paste0(tf, "_targets_DS_down")]] = TF_targets_DS_down

  
  ### identify targets reg in TF KD vs DS
  
  KD_comps = bulk_data$deseq_results[grepl(tf, names(bulk_data$deseq_results))]
  
  for (KD_comp in names(KD_comps)){
    
    deseq_res_KD_vs_DS = KD_comps[[KD_comp]][TF_targets,]
    deseq_res_KD_vs_DS_trend = deseq_res_KD_vs_DS[!is.na(deseq_res_KD_vs_DS$pvalue)&
                                                    deseq_res_KD_vs_DS$pvalue <0.1,]
    deseq_res_KD_vs_DS_diff = deseq_res_KD_vs_DS[!is.na(deseq_res_KD_vs_DS$padj)&
                                                   deseq_res_KD_vs_DS$padj <0.1,]
    
    TF_targets_trend_KD_up = intersect(TF_targets, 
                                       deseq_res_KD_vs_DS_trend$gene[deseq_res_KD_vs_DS_trend$log2FoldChange>0])
    pred_targets_reg_list[[paste0(KD_comp, "_", tf,  "_targets_trend_KD_up")]] = TF_targets_trend_KD_up
    TF_targets_KD_up = intersect(TF_targets, 
                                 deseq_res_KD_vs_DS_diff$gene[deseq_res_KD_vs_DS_diff$log2FoldChange>0])
    pred_targets_reg_list[[paste0(KD_comp, "_", tf,  "_targets_KD_up")]] = TF_targets_KD_up
    
    TF_targets_trend_KD_down = intersect(TF_targets, 
                                         deseq_res_KD_vs_DS_trend$gene[deseq_res_KD_vs_DS_trend$log2FoldChange<0])
    pred_targets_reg_list[[paste0(KD_comp, "_", tf,  "_targets_trend_KD_down")]] = TF_targets_trend_KD_down
    TF_targets_KD_down = intersect(TF_targets, 
                                   deseq_res_KD_vs_DS_diff$gene[deseq_res_KD_vs_DS_diff$log2FoldChange<0])
    pred_targets_reg_list[[paste0(KD_comp, "_", tf,  "_targets_KD_down")]] = TF_targets_KD_down
    
    
    ### identify targets with inverse change in DS vs CON and TF KD vs DS 
    
    pred_targets_reg_list[[paste0(KD_comp, "_", tf, "_targets_trend_DS_up_KD_down")]] = 
      intersect(TF_targets_trend_DS_up, TF_targets_trend_KD_down)
    pred_targets_reg_list[[paste0(KD_comp, "_", tf, "_targets_trend_DS_down_KD_up")]] = 
      intersect(TF_targets_trend_DS_down, TF_targets_trend_KD_up)
    
    pred_targets_reg_list[[paste0(KD_comp, "_", tf,  "_targets_DS_up_KD_down")]] = 
      intersect(TF_targets_DS_up, TF_targets_KD_down)
    pred_targets_reg_list[[paste0(KD_comp, "_", tf,  "_targets_DS_down_KD_up")]] = 
      intersect(TF_targets_DS_down, TF_targets_KD_up)
    
  }
  
}


### save table with DEGs by comparison

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


### save DESeq2 results for cand TFs by comparison

t1 = NULL

for (comp in names(bulk_data$comp_groups)){
  
  t2 = bulk_data$deseq_results[[comp]]
  t3 = t2[t2$gene %in% cand_TFs,]
  t4 = cbind(comp = comp, t3)
  
  t1 = rbind(t1, t4)
  
}

write_csv(t1, file = paste0(out_dir,script_ind, "DESeq_results_cand_TFs.csv"))



### plot heatmaps of all comparisons + cand TFs

pl_sets = c(list(cand_TFs = cand_TFs), pred_targets_reg_list)

meta = bulk_data$meta

pl_mat_X = bulk_data$gene_Z_scores$all_samples

pl_genes = intersect(rownames(pl_mat_X ), 
                     unlist(pl_sets))

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_pred_targets_reg.pdf"), 
    width = 15, height = 30)
{
  
  for (pl_set in names(pl_sets)){
    
    pl_mat_X = bulk_data$gene_Z_scores$all_samples
    
    pl_genes = intersect(rownames(pl_mat_X ), 
                         pl_sets[[pl_set]])
    
    if (length(pl_genes)>1){
      
      pl_mat_X = pl_mat_X[pl_genes,]
      
      p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                            pl_meta = meta,
                            pl_genes = pl_genes,
                            x_col = "sample", 
                            meta_annot_cols = c("group", "batch", "ASO_design", "genotype"),
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
#calculate group mean DELTA z-scores and p-value for each comparison
#########################################

meta = bulk_data$meta

sample_z_mat = bulk_data$gene_Z_scores$all_samples

comp_groups = bulk_data$comp_groups

for (comp in names(comp_groups)){
  
  comp_gr1 = comp_groups[[comp]][1]
  comp_gr2 = comp_groups[[comp]][2]
  
  m1 = matrix(nrow = nrow(sample_z_mat), ncol = 1, 
              dimnames = list(rownames(sample_z_mat), "Z"))
  
  m_comp_gr1 = m1
  m2 = sample_z_mat[,meta$group == comp_gr1]
  m_comp_gr1[,1] = apply(m2, 1, mean)
  
  m_comp_gr2 = m1
  m2 = sample_z_mat[,meta$group == comp_gr2]
  m_comp_gr2[,1] = apply(m2, 1, mean)
  
  m_DELTA = m_comp_gr1 - m_comp_gr2
  
  m_p = m1
  colnames(m_p) = c("p")
  t2 = bulk_data$deseq_results[[comp]]
  m_p[t2$gene,1] = t2$padj 
  
  l1 = list(comp_gr1 = comp_gr1,
            comp_gr2 = comp_gr2,
            m_comp_gr1 = m_comp_gr1,
            m_comp_gr2 = m_comp_gr2,
            m_DELTA = m_DELTA,
            m_p = m_p)
  
  bulk_data$gene_Z_scores_comp[[comp]] = l1
} 


#########################################
#plot group mean DELTA z-scores and p-value for each comparison for TF targets (genes on x axis)
#########################################

pl_sets = c(list(cand_TFs = cand_TFs), pred_targets_reg_list)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_pred_targets_reg_DELTA.pdf"), 
    width = 30, height = 5)
{
  
  for (pl_set in names(pl_sets)){
    
    for (comp in names(bulk_data$gene_Z_scores_comp)){
      
      l1 = bulk_data$gene_Z_scores_comp[[comp]]
      
      pl_mat_Z = t(l1$m_DELTA)
      
      #scale limits
      v1 = max(abs(pl_mat_Z))
      lims_Z = 0.5*c(-v1, v1)
      
      #select genes to plot
      pl_genes = intersect(colnames(pl_mat_Z), pl_sets[[pl_set]])
      
      
      if (length(pl_genes)>0){
        
        pl_mat_Z = t(as.matrix(pl_mat_Z[,pl_genes]))
        colnames(pl_mat_Z) = pl_genes
        
        # label cells with padj <0.1
        pl_mat_p = t(l1$m_p)
        pl_mat_p[] = ""
        pl_mat_p[t(l1$m_p)<0.1] = "*"
        pl_mat_p = t(as.matrix(pl_mat_p[,colnames(pl_mat_Z)]))
        
        #create plot
        
        p1 = pheatmap::pheatmap(pl_mat_Z, 
                                display_numbers = pl_mat_p,
                                cluster_rows = FALSE, cluster_cols = FALSE,
                                show_rownames = FALSE, show_colnames = TRUE,
                                color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                                breaks = seq(lims_Z[1], lims_Z[2], length.out = 251),
                                border_color = NA, cellwidth = 10, cellheight = 10, 
                                fontsize = 10, number_color = "grey70", fontsize_number = 10, 
                                main = paste0(pl_set, " - ", comp, " - DELTA Z-score"))
        
      }
    }
  }
}
dev.off()










message("\n####### End K02: ", Sys.time(),"\n")


