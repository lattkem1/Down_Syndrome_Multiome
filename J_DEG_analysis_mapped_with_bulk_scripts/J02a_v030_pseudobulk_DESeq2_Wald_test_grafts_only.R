message("\n\n##########################################################################\n",
        "# Start E02 Pseudobulk DESeq2 analysis (only grafts, clusters by ref mapping): ", Sys.time(),
        "\n##########################################################################\n",
        "\n   contrasts to extract DEGs by cluster",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DESeq2)


#specify script/output index as prefix for file names
script_ind = "J02a_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_mapped_data_vs_reference_with_neurons_in_vitro_bulk/")


#load group info and pseudobulk dataset for data to mapped on reference

gr_tab = read_csv(paste0(main_dir, "A_input/group_tab_grafts.csv"))

load(file = paste0(out_dir,"J01_bulk_data_mapped_data.rda")) 

bulk_data = bulk_data_mapped

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
  if (v2 < 6){p2 = matlab.like(6)[1:v2]} else {p2 = matlab.like(v2)}
  return(p2)
}



###########################################################
# remove clusters with each <2 pseudobulks per cluster_stage per group
###########################################################

message("\n\n          *** Preparing dataset for DESeq2 analysis... ", Sys.time(),"\n\n")

set.seed(12345)

t1 = bulk_data$meta
t2 = t1 %>% group_by(cluster_name, group) %>% summarise(N_pseudobulks = n())

comp_clusters = NULL

for (cl in bulk_data$grouping$cluster_name){
  t3 = t2[t2$cluster_name == cl,]
  if (nrow(t3)==length(bulk_data$grouping$group)){
    if(min(t3$N_pseudobulks >=2)){
      comp_clusters = c(comp_clusters, cl)
    }
  }
}

comp_meta = t1[t1$cluster_name %in% comp_clusters,]

# define grouping variable to extract comparisons by group for each cluster
comp_meta$cluster_group = paste0(comp_meta$cluster_name, "_", comp_meta$group)

rownames(comp_meta) = comp_meta$pseudobulk

bulk_data$meta = comp_meta


#remove genes from countmatrix with <0.1 counts/cell in all pseudobulks  

m1 = bulk_data$counts

pseudobulks = unique(comp_meta$cluster_sample)

N_pb_cells = lengths(bulk_data$cells)

for (pb in pseudobulks){
  m1[,pb] = m1[,pb]/N_pb_cells[pb]
}

keep_genes = rownames(m1)[apply(m1, 1, max)>0.1]

comp_counts = bulk_data$counts[keep_genes ,comp_meta$pseudobulk]



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

#extract stats for each comparison (genes signif dereg and with tendency of dereg)
deseq_res_list = NULL
genes_reg_list = NULL

for (cl in comp_clusters){
  message("  * extracting DESeq results for ", cl)
  
  cl_DS = paste0(cl, "_DS")
  cl_CON = paste0(cl, "_CON")
  t1 = as.data.frame(results(dds, contrast=c("cluster_group", cl_DS, cl_CON)))
  
  deseq_res_list[[paste0(cl, "_DS_vs_CON")]] = t1
  
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.05 & t1$log2FoldChange >= +log2(1.2)]
  genes_reg_list[[paste0(cl, "_DS_up")]] = v3
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.05 & t1$log2FoldChange <= -log2(1.2)]
  genes_reg_list[[paste0(cl, "_DS_down")]] = v3
  
}

bulk_data$deseq_results = deseq_res_list
bulk_data$DEGs = genes_reg_list


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



# save table with DEGs by comparison including TFs and HSA21 genes

l1 = genes_reg_list
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




message("\n\n##########################################################################\n",
        "# Completed E02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


