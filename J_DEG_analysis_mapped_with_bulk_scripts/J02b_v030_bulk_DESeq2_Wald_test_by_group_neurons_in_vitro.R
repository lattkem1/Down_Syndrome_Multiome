message("\n\n##########################################################################\n",
        "# Start J02: Bulk DESeq2 analysis (only bulk dataset for comaprison with sc): ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(biomaRt)

#specify script/output index as prefix for file names
script_ind = "J02b_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_mapped_data_vs_reference_with_neurons_in_vitro_bulk/")


#load group info and bulk dataset (convert ensembl IDs to gene symbols)

gr_tab_bulk = read_csv(paste0(main_dir, "A_input/group_tab_neurons_in_vitro_bulk_Peter24.csv"))
gr_tab_bulk$cluster_name = "NEU_bulk"
gr_tab_bulk$cluster_sample = paste0(gr_tab_bulk$cluster_name, "_", gr_tab_bulk$sample)
rownames(gr_tab_bulk) = gr_tab_bulk$cluster_sample

bulk_counts = read_csv(paste0(main_dir, "A_input/count_table_neurons_in_vitro_bulk_Peter24.csv"))

bulk_counts$gene_symbol = mapIds(org.Hs.eg.db,keys=bulk_counts$Gene,column="SYMBOL", 
                                 keytype="ENSEMBL", multiVals="first")

bulk_counts_mat = as.matrix(bulk_counts[!(colnames(bulk_counts) %in% c("Gene", "gene_symbol"))])
rownames(bulk_counts_mat) = bulk_counts$gene_symbol
bulk_counts_mat = bulk_counts_mat[!is.na(bulk_counts$gene_symbol),]
bulk_counts_mat = bulk_counts_mat[,gr_tab_bulk$sample]
colnames(bulk_counts_mat) = gr_tab_bulk$cluster_sample


#get marker gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol


#################################################
# DEG analysis for count data with DESeq2
#################################################

#run DESeq2 analysis (LRT test to assess full model vs reduced model to identify sign changes)

message("\n\n          *** Running DESeq2 analysis... ", Sys.time(),"\n\n")

gr_tab = gr_tab_bulk

dds = DESeqDataSetFromMatrix(bulk_counts_mat, colData = gr_tab, 
                             design = ~group)
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds = DESeq(dds) #standard Wald test


bulk_data = list(meta = gr_tab, counts = bulk_counts_mat, deseq_dataset = dds)

#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



#################################################
# Extract DESeq2 results
#################################################

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

#extract stats for each comparison

t1 = as.data.frame(results(dds, contrast=c("group", "DS", "CON")))

bulk_data$deseq_results$DS_vs_CON = t1

bulk_data$DEGs$DS_vs_CON_up = rownames(t1)[!is.na(t1$padj) & 
                                             t1$padj <= 0.10 & 
                                             t1$log2FoldChange > log2(1.2)]

bulk_data$DEGs$DS_vs_CON_down = rownames(t1)[!is.na(t1$padj) & 
                                             t1$padj <= 0.10 & 
                                             t1$log2FoldChange < -log2(1.2)]


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




message("\n\n##########################################################################\n",
        "# Completed J02b ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


