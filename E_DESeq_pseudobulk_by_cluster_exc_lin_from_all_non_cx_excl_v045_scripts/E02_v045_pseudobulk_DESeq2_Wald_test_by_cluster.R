message("\n\n##########################################################################\n",
        "# Start E02 Pseudobulk DESeq2 analysis: ", Sys.time(),
        "\n##########################################################################\n",
        "\n   contrasts to extract DEGs by cluster",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E14_241219_DS_foetal_brain_grafts_for_man_v02_low_string/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DESeq2)


#specify script/output index as prefix for file names
script_ind = "E02_"


#specify output directory
out_dir = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_v045/")

#load dataset
load(file = paste0(out_dir,"E01_bulk_data.rda")) 


#get marker gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol



###########################################################
# functions
###########################################################




###########################################################
# remove clusters with each <2 pseudobulks per cluster per group
###########################################################

message("\n\n          *** Preparing dataset for DESeq2 analysis... ", Sys.time(),"\n\n")

set.seed(12345)

t1 = bulk_data$meta
t2 = t1 %>% group_by(cluster_name, group) %>% summarise(N_pseudobulks = n())

comp_clusters = NULL

for (cl in unique(t1$cluster_name)){
  t3 = t2[t2$cluster_name == cl,]
  if (nrow(t3)==length(unique(t1$group))){
    if(min(t3$N_pseudobulks >=2)){
      comp_clusters = c(comp_clusters, cl)
    }
  }
}

comp_meta = t1[t1$cluster_name %in% comp_clusters,]

# define grouping variable to extract comparisons by group for each cluster
comp_meta$cluster_group = paste0(comp_meta$cluster_name, "_", comp_meta$group)

rownames(comp_meta) = comp_meta$cluster_sample

bulk_data$meta = comp_meta


#remove genes from countmatrix with <0.1 counts/cell in all pseudobulks  

m1 = bulk_data$counts

pseudobulks = unique(comp_meta$cluster_sample)

N_pb_cells = lengths(bulk_data$cells)

for (pb in pseudobulks){
  m1[,pb] = m1[,pb]/N_pb_cells[pb]
}

keep_genes = rownames(m1)[apply(m1, 1, max)>0.1]

t1 = bulk_data$counts

comp_counts = t1[keep_genes, match(comp_meta$cluster_sample, colnames(t1))]


#################################################
# DEG analysis for count data with DESeq2
#################################################

#run DESeq2 analysis (LRT test to assess full model vs reduced model to identify sign changes)

message("\n\n          *** Running DESeq2 analysis... ", Sys.time(),"\n\n")

dds = DESeqDataSetFromMatrix(comp_counts, colData = comp_meta, 
                             design = ~cluster_group)
dds = DESeq(dds) #standard Wald test

bulk_data$deseq_dataset = dds

#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



#################################################
# Extract DESeq2 results
#################################################

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

### extract vst norm expression matrix, calculate gene z-scores 

meta = bulk_data$meta

vst_mat = assay(vst(bulk_data$deseq_dataset))

bulk_data$vst_mat = vst_mat

#calculate Z-score per gene by pseudobulk (cluster_sample) for all clusters combined
cluster_sample_mat = t(apply(vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat)

bulk_data$gene_Z_scores$clusters_combined = cluster_sample_mat


### extract stats for each comparison (genes signif dereg and with tendency of dereg)

for (cl in comp_clusters){
  message("  * extracting DESeq results for ", cl)
  
  cl_DS = paste0(cl, "_DS")
  cl_CON = paste0(cl, "_CON")
  t1 = as.data.frame(results(dds, contrast=c("cluster_group", cl_DS, cl_CON)))
  
  bulk_data$deseq_results[[paste0(cl, "_DS_vs_CON")]] = t1
  
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.1 & t1$log2FoldChange >= +log2(1.2)]
  bulk_data$DEGs[[paste0(cl, "_DS_up")]] = v3
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.1 & t1$log2FoldChange <= -log2(1.2)]
  bulk_data$DEGs[[paste0(cl, "_DS_down")]] = v3
  
}

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




#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed E02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


