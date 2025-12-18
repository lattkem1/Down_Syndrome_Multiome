message("\n#################################\n",
        "####### Start J01 DESeq2 analysis: ", Sys.time(),
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
library(ggrepel)
library(colorRamps)
library(viridis)
library(pheatmap)


#specify script/output index as prefix for file names
script_ind = "K01b_"

#specify output directory
out_dir = paste0(main_dir,"K_expression_analysis_in_vitro_ASO_exp_v062/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#load group and file info 

gr_tab = read_csv("A_input/group_tab_in_vitro_bulk_ASO_exp_v060.csv")

gr_tab = gr_tab[gr_tab$isogenic_pair == "C9_C13" & gr_tab$batch == "WL037",]

gr_tab$target_design = paste0(gr_tab$target, "_", gr_tab$ASO_design)

gr_tab$target_design[gr_tab$target == "non_targeting"] = "non_targeting"
gr_tab$target_design[gr_tab$target == "no_ASO"] = "no_ASO"

gr_tab$cell_line_target_design = paste0(gr_tab$cell_line, "_",gr_tab$target_design)

gr_tab$group = gr_tab$cell_line_target_design

rownames(gr_tab) = gr_tab$sample


###load dataset, convert to matrix 

t1 = read_tsv(paste0("/rds/general/user/mlattke/projects/dsseq23/live/E21_250716_ASO_exp3_BGI/",
  "A_nfcore_rnaseq_output/star_salmon/salmon.merged.gene_counts.tsv"))
m1 = as.matrix(t1[,-c(1,2)])
rownames(m1) = t1[[2]]
m2 = m1
count_mat = round(m2)

t1 = read_tsv(paste0("/rds/general/user/mlattke/projects/dsseq23/live/E22_251112_ASO_exp4_BGI/",
                     "A_nfcore_rnaseq_output/star_salmon/salmon.merged.gene_counts.tsv"))
m1 = as.matrix(t1[,-c(1,2)])
rownames(m1) = t1[[2]]
m2 = m1
comb_genes = intersect(rownames(m2), rownames(count_mat))
count_mat = cbind(count_mat[comb_genes,], round(m2[comb_genes,]))

gr_tab$sample[!gr_tab$sample %in% colnames(count_mat)]


count_mat = count_mat[,gr_tab$sample]

bulk_data_in_vitro = list(meta = gr_tab, gr_tab = gr_tab, counts = count_mat)


### get marker gene panels

GOI = list()

t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol

#add all DS up/down genes in tissue

load("E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_v045/E03_bulk_data_w_expr_z_scores.rda")
bulk_data_tissue = bulk_data
l1 = bulk_data$DEGs

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


### define covariates to correct for
names(gr_tab)
covars = NULL


### define groups to compare
unique(gr_tab$group)

comps = list(C9_non_targeting_vs_no_ASO = c("C9_non_targeting", "C9_no_ASO"),
             C13_non_targeting_vs_no_ASO = c("C13_non_targeting", "C13_no_ASO"),
             C13_PKNOX1_D1_vs_no_ASO = c("C13_PKNOX1_D1", "C13_no_ASO"),
             C13_BACH1_D1_vs_no_ASO = c("C13_BACH1_D1", "C13_no_ASO"),
             C13_GABPA_D1_vs_no_ASO = c("C13_GABPA_D1", "C13_no_ASO"),
             C13_GABPA_D3_vs_no_ASO = c("C13_GABPA_D3", "C13_no_ASO")
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




###########################################################
# run DESeq2 analysis for group comps[i] vs comp_CON
###########################################################

message("\n\n          *** Preparing dataset for DESeq2 analysis... ", Sys.time(),"\n\n")

bulk_data = bulk_data_in_vitro


for (comp_name in names(comps)){
  
  message("\n\n          *** DESeq2 analysis comp ", comp_name, " - ", Sys.time(),"\n\n")
  
  comp1 = comps[[comp_name]][1]
  comp_CON = comps[[comp_name]][2]
  
  meta_comp = bulk_data$meta[bulk_data$meta$group %in% c(comp_CON, comp1),]
  meta_comp$group = factor(meta_comp$group, levels = c(comp_CON, comp1))
  
  counts_comp = bulk_data$counts[,meta_comp$sample]
  
  
  form_red = "~"
  
  for (cov1 in covars){
    
    if (length(unique(meta_comp[[cov1]]))>1){
      if (form_red == "~"){form_red = paste0(form_red, cov1)} else {
        form_red = paste0(form_red, " + ", cov1)}
    }
  }
  
  if (form_red == "~"){
    form_full = "~group"
    form_red = "~1"
  } else {
    form_full = paste0(form_red, " + group")
  }
  
  #run deseq2
  
  dds = DESeqDataSetFromMatrix(counts_comp, colData = meta_comp, 
                               design = as.formula(form_full))
  dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red)) 
  
  bulk_data$deseq_dataset[[comp_name]] = dds
  
}



#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



#################################################
# Extract DESeq2 results
#################################################

#load(file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

### extract DESeq2 results tables and DEGs 

for (comp1 in names(bulk_data$deseq_dataset)){
  
  dds = bulk_data$deseq_dataset[[comp1]]
  
  meta_comp = bulk_data$meta[bulk_data$meta$sample %in% colnames(dds), ]
  
  t0 = as.data.frame(results(dds))
  t1 = cbind(gene = rownames(t0), t0)
  bulk_data$deseq_results[[comp1]] = t1
  
  t2 = t1[!is.na(t1$padj) & t1$padj <= 0.1 & t1$log2FoldChange > 0,]
  bulk_data$DEGs[[paste0(comp1, "_up")]] = t2$gene
  
  t2 = t1[!is.na(t1$padj) & t1$padj <= 0.1 & t1$log2FoldChange < 0,]
  bulk_data$DEGs[[paste0(comp1, "_down")]] = t2$gene
  
}



### extract vst norm expression matrix (based on all clusters combined), correct for covariates, calculate gene z-scores 

dds = DESeqDataSetFromMatrix(bulk_data$counts, colData = bulk_data$meta, 
                                   design = as.formula(form_full))

vst_mat = assay(vst(dds))

bulk_data$vst_mat_uncorr = vst_mat

#batch-correct vst matrix

meta = bulk_data$meta

vst_mat_corr = vst_mat

for (cov1 in covars){
  if (length(unique(bulk_data$meta[[cov1]]))>1){
    vst_mat_corr = limma::removeBatchEffect(vst_mat_corr, batch = meta[[cov1]])
  } else {vst_mat_corr = vst_mat_corr}
}

bulk_data$vst_mat = vst_mat_corr

#calculate Z-score per gene by pseudobulk (cluster_sample) for all clusters combined (from uncorrected and corrected vst matrix)
cluster_sample_mat = t(apply(vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat)
bulk_data$gene_Z_scores_uncorr[["combined"]] = cluster_sample_mat

cluster_sample_mat = t(apply(vst_mat_corr, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat_corr)
bulk_data$gene_Z_scores[["combined"]] = cluster_sample_mat


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



# save table with DEGs by comparison including TFs and HSA21 genes

l1 = bulk_data$DEGs
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

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "DESeq2_up_down_by_cluster_genes.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "DESeq2_up_down_by_cluster_N_genes.csv"))



#################################################
# summarise all DEGs (all clusters combined)
#################################################

v1 = unique(unlist(bulk_data$DEGs))

t1 = tibble(N_DEGs_all = length(v1),
            N_DEGs_all_Chr21 = length(intersect(v1, GOI$Chr21)),
            N_DEGs_all_TF = length(intersect(v1, GOI$TF)),)

write_csv(t1, file = paste0(out_dir,script_ind, "DEGs_all_clusters_combined_N_genes.csv"))



#################################################
# Plot number of DEG (incl Chr21 genes) by cluster
#################################################

l1 = bulk_data$DEGs
l2 = lapply(l1, intersect, GOI$Chr21)
meta = bulk_data$meta

t1 = tibble(comp_group = names(l1), comp = NA, DS_up_down = NA,
            N_genes = lengths(l1), N_genes_Chr21 = lengths(l2))

for (comp1 in names(comps)){
  t1$comp[grepl(comp1, t1$comp_group)] =comp1
}

t1$up_down[grepl("up", t1$comp_group)] = "up"
t1$up_down[grepl("down", t1$comp_group)] = "down"

t1 = t1[!is.na(t1$comp),]

p1 = ggplot(t1, aes(x = comp, group = up_down))+
  scale_x_discrete(limits = names(comps))+
  scale_fill_manual(limits = c("down", "up"), values = c("blue", "red"))+
  scale_color_manual(limits = c("down", "up"), values = c("blue", "red"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2 = p1 + geom_col(aes(y = N_genes, fill = up_down), width = 0.7, position = position_dodge(width = 0.7))+
  labs(title = "All DEGs")

p3 = p1 + geom_col(aes(y = N_genes_Chr21, fill = up_down), width = 0.7, position = position_dodge(width = 0.7))+
  labs(title = "Chr21 DEGs")

p4 = p1 + geom_col(aes(y = N_genes, color = up_down), width = 0.7, position = position_dodge(width = 0.9),
                   fill = "grey80")+
  geom_col(aes(y = N_genes_Chr21, fill = up_down, color = DS_up_down), 
           width = 0.7, position = position_dodge(width = 0.9))+
  labs(title = "DEG Chr21 vs others")

pdf(file = paste0(out_dir,script_ind,"DEGs_numbers_by_comp.pdf"), 
    width = 5, height = 4)
{
  plot(p2)
  plot(p3)
  plot(p4)
}
dev.off()



#######################################
# plot PCA plot samples (by cluster, based on uncorrected vst matrix)
#######################################

meta = bulk_data$meta

samples = unique(bulk_data$gr_tab$sample)
gr = unique(bulk_data$gr_tab$group)

vst_mat = bulk_data$vst_mat_uncorr

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
var_genes = names(v1[order(-v1)][1:3000])


### plot PCA

pl = list()

m1 = vst_mat[var_genes,]
z_mat= scale(m1)

pc_analysis = prcomp(z_mat)

meta_pca = cbind(meta, pc_analysis$rotation)
meta_pca$sample_label = paste0(meta_pca$sample)

pl[[paste0("PC1_PC2_by_group")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
  scale_color_manual(limits = gr, values = pal(gr))+
  theme_minimal()+labs(title = paste0("PC1_PC2_by_group"))

pl[[paste0("PC1_PC2_by_batch")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = batch))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = batch))+
  scale_color_manual(limits = unique(meta_pca$batch), values = pal(unique(meta_pca$batch)))+
  theme_minimal()+labs(title = paste0("PC1_PC2_by_batch"))

pl[[paste0("PC1_PC2_by_target")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = target))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = target))+
  scale_color_manual(limits = unique(meta_pca$target), values = pal(unique(meta_pca$target)))+
  theme_minimal()+labs(title = paste0("PC1_PC2_by_target"))


pdf(file = paste0(out_dir,script_ind, "PCA_plot_samples_uncorr.pdf"), 
    width = 7, height = 6)
{
  lapply(pl, function(x){x})
}
dev.off()



#######################################
# plot PCA plot samples (by cluster, based on corrected vst matrix)
#######################################

meta = bulk_data$meta

samples = unique(bulk_data$gr_tab$sample)
gr = unique(bulk_data$gr_tab$group)

vst_mat = bulk_data$vst_mat

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
var_genes = names(v1[order(-v1)][1:3000])


### plot PCA

pl = list()

m1 = vst_mat[var_genes,]
z_mat= scale(m1)

pc_analysis = prcomp(z_mat)

meta_pca = cbind(meta, pc_analysis$rotation)
meta_pca$sample_label = paste0(meta_pca$sample)

pl[[paste0("PC1_PC2_by_group")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
  scale_color_manual(limits = gr, values = pal(gr))+
  theme_minimal()+labs(title = paste0("PC1_PC2_by_group"))

pl[[paste0("PC1_PC2_by_batch")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = batch))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = batch))+
  scale_color_manual(limits = unique(meta_pca$batch), values = pal(unique(meta_pca$batch)))+
  theme_minimal()+labs(title = paste0("PC1_PC2_by_batch"))

pl[[paste0("PC1_PC2_by_target")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = target))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = target))+
  scale_color_manual(limits = unique(meta_pca$target), values = pal(unique(meta_pca$target)))+
  theme_minimal()+labs(title = paste0("PC1_PC2_by_target"))


pdf(file = paste0(out_dir,script_ind, "PCA_plot_samples_corr.pdf"), 
    width = 7, height = 6)
{
  lapply(pl, function(x){x})
}
dev.off()



#################################################
# plot heatmap top3000 variable genes by cluster_sample
#################################################

vst_mat = bulk_data$vst_mat

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
pl_genes = names(v1[order(-v1)][1:3000])

pl_mat_X = bulk_data$gene_Z_scores$combined[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

names(meta)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_var_genes_top3000.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "batch", "ASO_design",
                                            "target"),
                        show_rownames = FALSE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 5, cellheight = 0.05, 
                        fontsize = 5, title = paste0("Top3000 variable genes - Z-score"))
  
}
dev.off()



#################################################
# plot heatmap DEGs combined by cluster_sample
#################################################

vst_mat = bulk_data$vst_mat

pl_genes = unique(unlist(bulk_data$DEGs))

pl_mat_X = bulk_data$gene_Z_scores$combined[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

names(meta)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_comb.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "batch", "ASO_design",
                                            "target"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 5, cellheight = 5, 
                        fontsize = 5, title = paste0("Combined DEGs - Z-score"))
  
}
dev.off()







message("\n\n##########################################################################\n",
        "# Completed J01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

