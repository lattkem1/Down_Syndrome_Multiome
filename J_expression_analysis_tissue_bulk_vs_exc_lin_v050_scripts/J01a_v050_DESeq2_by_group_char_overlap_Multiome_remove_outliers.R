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

# Open packages necessary for analysis.
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(colorRamps)
library(viridis)
library(pheatmap)



#specify script/output index as prefix for file names
script_ind = "J01a_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_tissue_bulk_vs_exc_lin_v050/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#load group and file info 

gr_tab = read_csv("A_input/group_tab_tissue_bulk.csv")
rownames(gr_tab) = gr_tab$sample

gr_tab = gr_tab[!(gr_tab$sample %in% c("B12D2")), ]


###load dataset, convert to matrix 

t1 = read_tsv(paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                     "E18_250305_DS_tissue_PCW11_13_bulk_RNAseq/A_nfcore_rnaseq_output/", 
                     "star_salmon/salmon.merged.gene_counts.tsv"))

m1 = as.matrix(t1[,-c(1,2)])
rownames(m1) = t1[[2]]
m2 = m1[,gr_tab$IGF_ID]
colnames(m2) = gr_tab$sample
count_mat = round(m2)



### get marker gene panels

GOI = list()

t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol

#add all DS up/down genes in tissue

load("E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_PCW11_13_v045/E03_bulk_data_w_expr_z_scores.rda")
bulk_data_snMultiome = bulk_data
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

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
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
  
  p1 = pheatmap::pheatmap(pl_mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                          show_rownames = show_rownames, show_colnames = show_colnames,
                          color = color,
                          breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                          annotation_col = annot_col, annotation_colors = annot_colors,
                          border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                          fontsize = fontsize, main = title
  )
  
  return(p1)
  
}





#################################################
# DEG analysis for count data with DESeq2
#################################################

#remove genes from countmatrix with <0.5 cpm in all samples  

m1 = count_mat/colSums(count_mat)*1e6
v1 = apply(m1,1,max)

comp_counts = count_mat[v1>=0.1,]

comp_counts = apply(comp_counts,c(1,2), as.integer)

gr_tab$group = factor(gr_tab$group, levels = unique(gr_tab$group))
gr_tab$sample = factor(gr_tab$sample, levels = gr_tab$sample)
gr_tab = as.data.frame(gr_tab)
rownames(gr_tab) = gr_tab$sample


#run DESeq2

message("\n\n          ***Gen matrix")
dds = DESeqDataSetFromMatrix(comp_counts, colData = gr_tab, 
                             design = ~group, ignoreRank = TRUE)

message("\n\n          ***Run DESeq2")

dds = DESeq(dds)

message("\n\n          ***DESeq2 completed \n\n")


#collect results and save


bulk_data = list(counts = comp_counts,
                 meta = gr_tab,
                 gr_tab = gr_tab)

bulk_data$deseq_dataset[["bulk_PCW11_14"]] = dds


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



#################################################
# Extract DESeq2 results
#################################################

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

#load(file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 


#extract vst normalised expression data, calculate Z-score

dds = bulk_data$deseq_dataset$bulk_PCW11_14

vst_mat = assay(vst(dds))

bulk_data$vst_mat = vst_mat

z_mat = t(apply(vst_mat, 1, scale))
colnames(z_mat) = colnames(vst_mat)

bulk_data$gene_Z_scores$bulk_PCW11_14 = z_mat 



#extract stats and DEGs for each comparison 

t1 = as.data.frame(results(dds, contrast=c("group", "DS", "CON")))

bulk_data$deseq_results[["bulk_PCW11_14"]] = t1

v3 = rownames(t1)[!is.na(t1$padj) & 
                    t1$padj <= 0.10 & t1$log2FoldChange >= +log2(1.2)]
bulk_data$DEGs[[paste0("bulk_PCW11_14_up")]] = v3
v3 = rownames(t1)[!is.na(t1$padj) & 
                    t1$padj <= 0.10 & t1$log2FoldChange <= -log2(1.2)]
bulk_data$DEGs[[paste0("bulk_PCW11_14_down")]] = v3


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



# save table with DEGs by comparison including TFs and HSA21 genes

l1 = bulk_data$DEGs
l3 = l1

for (goi in names(GOI)){
  
  l2 = lapply(l1, function(x){x = x[x%in% GOI[[goi]] ]})
  names(l2) = paste0(names(l1), "_", goi)
  l3 = c(l3, l2)
  
}

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
# Plot number of DEG (incl Chr21 genes) by comp_group
#################################################

l1 = bulk_data$DEGs
l2 = lapply(l1, intersect, GOI$Chr21)
meta = bulk_data$meta

t1 = tibble(comp = NA, comp_group = names(l1), DS_up_down = NA,
            N_genes = lengths(l1), N_genes_Chr21 = lengths(l2))

t1$DS_up_down[grepl("_up", t1$comp_group)] = "up"
t1$DS_up_down[grepl("_down", t1$comp_group)] = "down"

for (i in 1:nrow(t1)){t1$comp[i] = str_remove_all(t1$comp_group[i], paste0("_", t1$DS_up_down[i]))}

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

pdf(file = paste0(out_dir,script_ind,"DEGs_numbers_by_comp_group.pdf"), 
    width = 2, height = 3)
{
  plot(p2)
  plot(p3)
  plot(p4)
}
dev.off()




#######################################
# plot sample distance heatmap and PCA plot
#######################################

meta = bulk_data$meta

samples = unique(bulk_data$gr_tab$sample)
gr = unique(bulk_data$gr_tab$group)


### plot cluster_sample distance/correlation heatmap (by cluster)

vst_mat = bulk_data$vst_mat

z_mat = bulk_data$gene_Z_scores$bulk_PCW11_14

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
var_genes = names(v1[order(-v1)][1:3000])

#calculate sample distance matrix 
sampleDists <- dist(t(vst_mat[var_genes,]))
sampleDistMatrix <- as.matrix(sampleDists)

#calculate sample correlation matrix 
cor_mat = cor(vst_mat[var_genes,])

#plot sample distance and corr matrices


pdf(file = paste0(out_dir,script_ind, "Heatmap_sample_distance_corr_matrix.pdf"), 
    width = 12, height = 12)
{
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=viridis(250),
           cellwidth = 10, cellheight = 10,
           main = paste0("Sample Distance heatmap"))
  
  #plot sample distance matrix 
  
  pheatmap(cor_mat,
           color = colorRampPalette(c("blue", "white", "red"))(250),
           breaks = seq(0, 1, length.out = 251),
           cellwidth = 10, cellheight = 10,
           main = paste0("Sample Correlation heatmap"))
  
}
dev.off()


### plot PCA

pc_analysis = prcomp(z_mat[var_genes,])

meta_pca = cbind(meta, pc_analysis$rotation)
meta_pca$sample_label = paste0(meta_pca$sample)

p1 = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
  scale_color_manual(limits = gr, values = pal(gr))+
  theme_minimal()+labs(title = paste0("PCA"))



pdf(file = paste0(out_dir,script_ind, "PCA_plot_samples.pdf"), 
    width = 5, height = 4)
{
  plot(p1)
}
dev.off()




###########################################################
# plot volcano plots with selected gene sets highlighted
###########################################################

pl = list()

t1 = bulk_data$deseq_results$bulk_PCW11_14
t1 = t1[!is.na(t1$padj),]
t1$log10padj = log10(t1$padj)
t1$DEG = t1$padj <= 0.1 & abs(t1$log2FoldChange)>= log2(1.2)
t1$top_DEGs = abs(t1$log2FoldChange) > 0.5*max(abs(t1$log2FoldChange))|
  -t1$log10padj > 0.5*max(-t1$log10padj)

l1 = lapply(GOI, intersect, GOI$TF)
names(l1) = paste0(names(GOI), "_TFs")

l2 = c(GOI, l1)

for (goi_set in names(l2)){
  
  goi = l2[[goi_set]]
  
  t1$highlight = rownames(t1) %in% goi
  
  t1$gene_cat = "Other_gene"
  t1$gene_cat[t1$DEG] = "DEG"
  t1$gene_cat[rownames(t1) %in% goi] = goi_set
  
  #define plot labels (if more than 30 genes in gene set, keep only each top 10 sign and up/down genes)
  
  t1$plot_label = ""
  
  
  t2 = t1[t1$highlight,]
  
  if (nrow(t2)>30){
    
    t3 = t2[order(t2$log10padj),]
    t4 = t2[order(t2$log2FoldChange),]
    t5 = t2[order(-t2$log2FoldChange),]
    v5 = unique(c(rownames(t3)[1:10], rownames(t4)[1:10], rownames(t5)[1:10]))
    
    t1$plot_label[rownames(t1) %in% v5] = rownames(t1)[rownames(t1) %in% v5]
    
  } else {
    
    t1$plot_label[rownames(t1) %in% rownames(t2)] = rownames(t1)[rownames(t1) %in% rownames(t2)]
  }
  
  t3 = rbind(t1[!t1$highlight,], t1[t1$highlight,])
  
  
  pl[[goi_set]] = ggplot(t3, aes(x = log2FoldChange, y = -log10padj, color = gene_cat))+
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_hline(yintercept = -log10(0.1), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_point(aes(size = highlight), alpha = 0.8)+
    geom_label_repel(aes(label = plot_label), seed = 42, min.segment.length = 0, max.overlaps = Inf,
                     max.time = 5)+
    scale_size_manual(limits = c(FALSE, TRUE), values = c(1, 2))+
    scale_color_manual(limits = c(goi_set, "DEG", "Other_gene"), 
                       values = c("red", "black", "grey"))+
    theme_minimal()+
    labs(title = goi_set)
  
  
  ### plot zoomed plot only logFC -1 to 1, -log10padj <2
  
  t1$plot_label = ""
  
  t2 = t1[t1$highlight & -t1$log10padj <2 & abs(t1$log2FoldChange)<1,]
  
  if (nrow(t2)>30){
    
    t3 = t2[order(t2$log10padj),]
    t4 = t2[order(t2$log2FoldChange),]
    t5 = t2[order(-t2$log2FoldChange),]
    v5 = unique(c(rownames(t3)[1:10], rownames(t4)[1:10], rownames(t5)[1:10]))
    
    t1$plot_label[rownames(t1) %in% v5] = rownames(t1)[rownames(t1) %in% v5]
    
  } else {
    
    t1$plot_label[rownames(t1) %in% rownames(t2)] = rownames(t1)[rownames(t1) %in% rownames(t2)]
  }
  
  t3 = rbind(t1[!t1$highlight,], t1[t1$highlight,])
  t3 = t3[-t3$log10padj <2 & abs(t3$log2FoldChange)<1,]
  
  
  pl[[paste0(goi_set, "_zoomed")]] = ggplot(t3, aes(x = log2FoldChange, y = -log10padj, color = gene_cat))+
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_hline(yintercept = -log10(0.1), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_point(aes(size = highlight), alpha = 0.8)+
    geom_label_repel(aes(label = plot_label), seed = 42, min.segment.length = 0, max.overlaps = Inf,
                     max.time = 5)+
    scale_size_manual(limits = c(FALSE, TRUE), values = c(1, 2))+
    scale_color_manual(limits = c(goi_set, "DEG", "Other_gene"), 
                       values = c("red", "black", "grey"))+
    coord_cartesian(xlim = c(-1,1), ylim = c(0, 2))+
    theme_minimal()+
    labs(title = paste0(goi_set, "_zoomed"))
  
}


pdf(file = paste0(out_dir,script_ind, "Volcano_plots_DEGs_by_GOI_set.pdf"), 
    width = 10, height = 8)
{
  lapply(pl, function(x){x})
}
dev.off()




#######################################
# plot all DEGs from snMultiome
#######################################

pl_genes = unique(c(GOI$DS_up_any_cluster, GOI$DS_down_any_cluster))

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs.pdf"), 
    width = 15, height = 70)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs_smaller.pdf"), 
    width = 15, height = 15)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = FALSE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 0.5, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()


#######################################
# plot all DEGs up in  snMultiome
#######################################

pl_genes = GOI$DS_up_any_cluster

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs_up.pdf"), 
    width = 15, height = 70)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs_up_smaller.pdf"), 
    width = 15, height = 15)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = FALSE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 0.5, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()


#######################################
# plot all DEGs down in  snMultiome
#######################################

pl_genes = GOI$DS_down_any_cluster

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs_down.pdf"), 
    width = 15, height = 70)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs_down_smaller.pdf"), 
    width = 15, height = 15)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = FALSE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 0.5, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()

#######################################
# plot all concordant DEGs from snMultiome and bulk
#######################################

v1 = intersect(bulk_data$DEGs$bulk_PCW11_14_up, GOI$DS_up_any_cluster)
v2 = intersect(bulk_data$DEGs$bulk_PCW11_14_down, GOI$DS_down_any_cluster)

pl_genes = c(v1, v2)

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_concord_bulk_sn_multiome.pdf"), 
    width = 15, height = 20)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()




#######################################
# plot TFs diff expr in snMultiome
#######################################

pl_genes = intersect(GOI$TF, c(GOI$DS_up_any_cluster, GOI$DS_down_any_cluster))

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs_TFs.pdf"), 
    width = 15, height = 15)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()





#######################################
# plot Chr21 genes diff expr in snMultiome
#######################################

pl_genes = intersect(GOI$Chr21, c(GOI$DS_up_any_cluster, GOI$DS_down_any_cluster))

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_sn_multiome_DEGs_Chr21_genes.pdf"), 
    width = 15, height = 15)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()



#######################################
# plot bulk_DEGs not in snMultiome
#######################################

v1 = unlist(bulk_data$DEGs)
pl_genes = v1[!(v1 %in% c(GOI$DS_up_any_cluster, GOI$DS_down_any_cluster))]

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_bulk_DEGs_not_in_sn_multiome_DEGs.pdf"), 
    width = 15, height = 15)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()



#######################################
# plot bulk_DEGs not detected in snMultiome
#######################################

v1 = unlist(bulk_data$DEGs)
pl_genes = v1[!(v1 %in% rownames(bulk_data_snMultiome$vst_mat))]

pl_mat_X = bulk_data$gene_Z_scores$bulk_PCW11_14

pl_genes = intersect(rownames(pl_mat_X), pl_genes)

pl_mat_X = pl_mat_X[pl_genes,]

lims_X = 0.5*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_bulk_DEGs_not_in_sn_multiome_genes.pdf"), 
    width = 15, height = 15)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = bulk_data$meta,
                        pl_genes = pl_genes,
                        x_col = "sample", 
                        meta_annot_cols = c("group", "dev_PCW"),
                        show_rownames = TRUE, show_colnames = TRUE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("snMultiome DEGs - Z-score"))
  
}
dev.off()




###########################################################
# identify overlaps between tissue DEGs and graft DEGs
###########################################################

overlap_mat = matrix(nrow = length(bulk_data_snMultiome$DEGs), ncol = length(bulk_data$DEGs),
                     dimnames = list(names(bulk_data_snMultiome$DEGs), names(bulk_data$DEGs)))

overlap_mat_N = overlap_mat
overlap_mat_fract = overlap_mat

DEG_sets_snMultiome = "RG_s11_DS_down"

for (DEG_sets_snMultiome in names(bulk_data_snMultiome$DEGs)){
  
  for (DEG_sets_bulk in names(bulk_data$DEGs)){
    
    DEGs_snMultiome = bulk_data_snMultiome$DEGs[[DEG_sets_snMultiome]]
    DEGs_bulk = bulk_data$DEGs[[DEG_sets_bulk]]
    
    overlap_mat_N[DEG_sets_snMultiome, DEG_sets_bulk] = length(intersect(DEGs_snMultiome, DEGs_bulk))
    overlap_mat_fract[DEG_sets_snMultiome, DEG_sets_bulk] = length(intersect(DEGs_snMultiome, DEGs_bulk))/
      length(DEGs_snMultiome)
    
  }
}

overlap_mat_fract[is.nan(overlap_mat_fract)] = 0


###plot overlap stats

pdf(file = paste0(out_dir,script_ind, "Overlap_DEGs_sn_Multiome_by_cluster_vs_bulk.pdf"), 
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
           main = "Number of DEGs by cluster vs reg in models")
  
  pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(overlap_mat_N), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Number of DEGs by cluster vs reg in models")
  
  
  pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(na.omit(overlap_mat_fract)), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Fraction of DEGs by cluster vs reg in models")
  
  pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(na.omit(overlap_mat_fract)), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Fraction of DEGs by cluster vs reg in models")
  
}

dev.off()






message("\n\n##########################################################################\n",
        "# Completed J01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

