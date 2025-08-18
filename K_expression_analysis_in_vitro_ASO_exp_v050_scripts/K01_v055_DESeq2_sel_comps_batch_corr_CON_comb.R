message("\n#################################\n",
        "####### Start K01 DESeq2 analysis: ", Sys.time(),
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


#specify script/output index as prefix for file names
script_ind = "K01_"

#specify output directory
out_dir = paste0(main_dir,"K_expression_analysis_in_vitro_ASO_exp_v050/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}


#load group and file info 

gr_tab = read_csv("A_input/group_tab_in_vitro_bulk_ASO_exp.csv")
rownames(gr_tab) = gr_tab$sample

gr_tab$group = paste0(gr_tab$genotype, "_",gr_tab$target, "_KD")
gr_tab$group[gr_tab$target == "none"] = gr_tab$genotype[gr_tab$target == "none"]


###load dataset, convert to matrix 

t1 = read_tsv(paste0("/rds/general/user/mlattke/projects/dsseq23/live/E21_250716_ASO_exp3_BGI/",
  "A_nfcore_rnaseq_output/star_salmon/salmon.merged.gene_counts.tsv"))
m1 = as.matrix(t1[,-c(1,2)])
rownames(m1) = t1[[2]]
m2 = m1[,gr_tab$sample]
count_mat = round(m2)


#define covariates to correct for
covars = c("batch")


### get marker gene panels

GOI = list()

t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol

#define comparisons

unique(gr_tab$group)

comp_groups = list(DS_vs_CON = c("DS", "CON"),
                   PKNOX1_KD_vs_DS = c("DS_PKNOX1_KD", "DS"),
                   BACH1_KD_vs_DS = c("DS_BACH1_KD", "DS"),
                   GABPA_KD_vs_DS = c("DS_GABPA_KD", "DS")
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
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}



#################################################
# DEG analysis for count data with DESeq2
#################################################

#remove genes from countmatrix with <0.5 cpm in all samples  

m1 = count_mat/colSums(count_mat)*1e6
v1 = apply(m1,1,max)

count_mat = count_mat[v1>=0.5,]

count_mat = apply(count_mat,c(1,2), as.integer)

gr_tab$group = factor(gr_tab$group, levels = unique(gr_tab$group))
gr_tab$sample = factor(gr_tab$sample, levels = gr_tab$sample)
gr_tab = as.data.frame(gr_tab)
rownames(gr_tab) = gr_tab$sample

bulk_data = list(gr_tab = gr_tab,
                 meta = gr_tab,
                 counts = count_mat,
                 comp_groups = comp_groups
                )

### run DESeq2 analysis by cluster (LRT test to assess full model vs reduced model to identify sign changes)

message("\n\n          *** Running DESeq2 analysis... ", Sys.time(),"\n\n")

comp = names(comp_groups)[1]

for (comp in names(comp_groups)){
  
  message("\n          *** Running DESeq2 analysis ", comp, " - ", Sys.time(),"\n")
  
  meta_comp = gr_tab[gr_tab$group %in% unlist(comp_groups[[comp]]),]
  counts_comp = count_mat[, meta_comp$sample]
  
  #check whether covariates have multiple levels and generate formula for DESeq2 model with relevant covariates
  
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
  
  bulk_data$deseq_dataset[[comp]] = dds
  
}



#################################################
# Extract DESeq2 results
#################################################

#load(file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

### extract DESeq2 results tables and DEGs 

for (comp in names(bulk_data$deseq_dataset)){
  
  dds = bulk_data$deseq_dataset[[comp]]
  
  meta_comp = bulk_data$meta[bulk_data$meta$cluster_sample %in% colnames(dds), ]
  
  t0 = as.data.frame(results(dds))
  t1 = cbind(gene = rownames(t0), t0)
  bulk_data$deseq_results[[comp]] = t1
  
  t2 = t1[!is.na(t1$padj) & t1$padj <= 0.1 & t1$log2FoldChange > 0,]
  bulk_data$DEGs[[paste0(comp, "_up")]] = t2$gene
  
  t2 = t1[!is.na(t1$padj) & t1$padj <= 0.1 & t1$log2FoldChange < 0,]
  bulk_data$DEGs[[paste0(comp, "_down")]] = t2$gene
  
}



### extract vst norm expression matrix (based on all samples combined), correct for covariates, calculate gene z-scores 

count_mat = bulk_data$counts
meta = bulk_data$meta

dds = DESeqDataSetFromMatrix(count_mat, colData = meta, 
                             design = ~1)

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
bulk_data$gene_Z_scores_uncorr$all_samples = cluster_sample_mat

cluster_sample_mat = t(apply(vst_mat_corr, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat_corr)
bulk_data$gene_Z_scores$all_samples = cluster_sample_mat


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



### save table with DEGs by comparison including TFs and HSA21 genes

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

t1 = tibble(comp = NA, comp_group = names(l1), up_down = NA,
            N_genes = lengths(l1), N_genes_Chr21 = lengths(l2))

t1$up_down[grepl("_up", t1$comp_group)] = "up"
t1$up_down[grepl("_down", t1$comp_group)] = "down"

for (i in 1:nrow(t1)){t1$comp[i] = str_remove_all(t1$comp_group[i], paste0("_", t1$up_down[i]))}

p1 = ggplot(t1, aes(x = comp, group = up_down))+
  scale_x_discrete(limits = unique(t1$comp))+
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
  geom_col(aes(y = N_genes_Chr21, fill = up_down, color = up_down), 
           width = 0.7, position = position_dodge(width = 0.9))+
  labs(title = "DEG Chr21 vs others")

pdf(file = paste0(out_dir,script_ind,"DEGs_numbers_by_comp_group.pdf"), 
    width = 5, height = 4)
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

z_mat = bulk_data$gene_Z_scores$all_samples

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
  labs(title = paste0("PCA"))



pdf(file = paste0(out_dir,script_ind, "PCA_plot_samples.pdf"), 
    width = 9, height = 8)
{
  plot(p1)
}
dev.off()






message("\n\n##########################################################################\n",
        "# Completed K01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

