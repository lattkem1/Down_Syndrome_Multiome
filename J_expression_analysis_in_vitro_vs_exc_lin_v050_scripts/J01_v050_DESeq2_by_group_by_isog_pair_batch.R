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


#specify script/output index as prefix for file names
script_ind = "J01_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_in_vitro_vs_exc_lin_v050/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#load group and file info 

gr_tab = read_csv("A_input/group_tab_in_vitro_bulk_v050.csv")
rownames(gr_tab) = gr_tab$sample


###load dataset, convert to matrix 

t1 = read_tsv(paste0("/rds/general/user/mlattke/projects/dsseq23/live/Singleron_raw_data/",
  "LILAC1_NPC_vs_NEU_lines_comp_singleron_proj1/LILAC1-Pool-R_matrix.tsv.gz"))
m1 = as.matrix(t1[,-c(1)])
rownames(m1) = t1[[1]]
m2 = m1[,colnames(m1) %in% gr_tab$sample]
count_mat = round(m2)

t1 = read_tsv(paste0("/rds/general/user/mlattke/projects/dsseq23/live/Singleron_raw_data/",
  "LILAC3_250306_ASO_exp_singleron_proj3/Pool-6-R/Pool-6-R_matrix.tsv.gz"))
m1 = as.matrix(t1[,-c(1)])
rownames(m1) = t1[[1]]
m2 = m1[,colnames(m1) %in% gr_tab$sample]
comb_genes = intersect(rownames(m2), rownames(count_mat))
count_mat = cbind(count_mat[comb_genes,], round(m2[comb_genes,]))

t1 = read_tsv(paste0("/rds/general/user/mlattke/projects/dsseq23/live/Singleron_raw_data/",
  "LILAC3_250306_ASO_exp_singleron_proj3/Pool-72-R/Pool-72-R_matrix.tsv.gz"))
m1 = as.matrix(t1[,-c(1)])
rownames(m1) = t1[[1]]
m2 = m1[,colnames(m1) %in% gr_tab$sample]
comb_genes = intersect(rownames(m2), rownames(count_mat))
count_mat = cbind(count_mat[comb_genes,], round(m2[comb_genes,]))

count_mat = count_mat[,gr_tab$sample]


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


#define comparisons

gr_tab$comp = paste0(gr_tab$cell_type, "_", gr_tab$isogenic_pair,"_batch_", gr_tab$batch)
gr_tab$comp_group = paste0(gr_tab$comp, "_", gr_tab$group)




###########################################################
# functions
###########################################################




#################################################
# DEG analysis for count data with DESeq2
#################################################

#remove genes from countmatrix with <0.5 cpm in all samples  

m1 = count_mat/colSums(count_mat)*1e6
v1 = apply(m1,1,max)

comp_counts = count_mat[v1>=0.5,]

comp_counts = apply(comp_counts,c(1,2), as.integer)

gr_tab$comp_group = factor(gr_tab$comp_group, levels = unique(gr_tab$comp_group))
gr_tab$sample = factor(gr_tab$sample, levels = gr_tab$sample)
gr_tab = as.data.frame(gr_tab)
rownames(gr_tab) = gr_tab$sample


#run DESeq2

message("\n\n          ***Gen matrix")
dds = DESeqDataSetFromMatrix(comp_counts, colData = gr_tab, 
                             design = ~comp_group, ignoreRank = TRUE)

message("\n\n          ***Run DESeq2")

dds = DESeq(dds)

message("\n\n          ***DESeq2 completed \n\n")


#collect results and save


bulk_data = list(counts = comp_counts,
                meta = gr_tab,
                deseq_dataset = dds)


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



#################################################
# Extract DESeq2 results
#################################################

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

#load(file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 


#extract vst normalised expression data, calculate Z-score

dds = bulk_data$deseq_dataset

vst_mat = assay(vst(dds))

bulk_data$vst_mat = vst_mat

z_mat = t(apply(vst_mat, 1, scale))
colnames(z_mat) = colnames(vst_mat)

bulk_data$gene_Z_scores$combined = z_mat 



#extract stats and DEGs for each comparison 
comps = unique(gr_tab$comp)

comp = comps[1]

for (comp in comps){
  message("  * extracting DESeq results for ", comp)
  
  comp_groups = as.character(unique(gr_tab$comp_group[gr_tab$comp == comp]))
  
  t1 = as.data.frame(results(dds, contrast=c("comp_group", comp_groups[2], comp_groups[1])))
  
  bulk_data$deseq_results[[comp]] = t1
  
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.10 & t1$log2FoldChange >= +log2(1.2)]
  bulk_data$DEGs[[paste0(comp,"_up")]] = v3
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.10 & t1$log2FoldChange <= -log2(1.2)]
  bulk_data$DEGs[[paste0(comp,"_down")]] = v3
}

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
    width = 5, height = 4)
{
  plot(p2)
  plot(p3)
  plot(p4)
}
dev.off()





message("\n\n##########################################################################\n",
        "# Completed J01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

