message("\n\n##########################################################################\n",
        "# Start J02: Pseudobulk DESeq2 analysis (combined mapped + ref data): ", Sys.time(),
        "\n##########################################################################\n",
        "\n   contrasts to extract DEGs by cluster split by ref and mapped data",
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
library(org.Hs.eg.db)

#specify script/output index as prefix for file names
script_ind = "J02_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_mapped_data_vs_reference_with_neurons_in_vitro_bulk/")

#load group info and pseudobulk dataset for data to mapped on reference

gr_tab = read_csv(paste0(main_dir, "A_input/group_tab_grafts.csv"))

load(file = paste0(out_dir,"J01_bulk_data_mapped_data.rda")) 



#load group info and bulk dataset (convert ensembl IDs to gene symbols)

gr_tab_bulk = read_csv(paste0(main_dir, "A_input/group_tab_neurons_in_vitro_bulk_Peter24.csv"))
gr_tab_bulk$cluster_name = "NEU_bulk"
gr_tab_bulk$cluster_sample = paste0(gr_tab_bulk$cluster_name, "_", gr_tab_bulk$sample)

bulk_counts = read_csv(paste0(main_dir, "A_input/count_table_neurons_in_vitro_bulk_Peter24.csv"))

bulk_counts$gene_symbol = mapIds(org.Hs.eg.db,keys=bulk_counts$Gene,column="SYMBOL", 
                                 keytype="ENSEMBL", multiVals="first")

bulk_counts_mat = as.matrix(bulk_counts[!(colnames(bulk_counts) %in% c("Gene", "gene_symbol"))])
rownames(bulk_counts_mat) = bulk_counts$gene_symbol
bulk_counts_mat = bulk_counts_mat[!is.na(bulk_counts$gene_symbol),]
bulk_counts_mat = bulk_counts_mat[,gr_tab_bulk$sample]
colnames(bulk_counts_mat) = gr_tab_bulk$cluster_sample

# load reference dataset and cell type annotations

clust_tab_ref = read_csv(paste0(main_dir, "B_basic_analysis/B03_cluster_assignment_all.csv"))

gr_tab_ref = read_csv(file =paste0(main_dir, "B_basic_analysis/B02_gr_tab_filtered.csv"))

load(file = paste0(main_dir,"E_DESeq_pseudobulk_all_by_cluster/E03_bulk_data_w_expr_z_scores.rda")) 
bulk_data_ref = bulk_data



#get marker gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol

t1 = bulk_counts_mat[rownames(bulk_counts_mat) %in% GOI$Chr21,]



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
# combine reference and mapped pseudobulk datasets
###########################################################

###merge metadata, add cell_type for mapped data, add grouping variable for DESeq2 (cluster_name x sample_type x group)

meta_ref = bulk_data_ref$meta
meta_mapped = bulk_data_mapped$meta

t1 = rbind( meta_ref, meta_mapped, gr_tab_bulk)

t1$cell_type = meta_ref$cell_type[match(t1$cluster_name, meta_ref$cluster_name )]
t1$cell_type[t1$sample_type == "in_vitro"] = "NEU_bulk"

t1$cluster_sample = paste0(t1$cluster_name, "_", t1$sample)

t1$comp = paste0(t1$cluster_name, "_", t1$sample_type)
t1$comp_group = paste0(t1$comp, "_", t1$group)

meta_comb = t1

### merge count tables
m1 = cbind(bulk_data_ref$counts, bulk_data_mapped$counts)
m2 = m1[rownames(m1) %in% rownames(bulk_counts_mat),]
m3 = bulk_counts_mat[match(rownames(m2), rownames(bulk_counts_mat)),]
count_mat_comb = cbind(m2, m3)

###create merged bulk_data dataset

bulk_data = list(meta = meta_comb, cells = c(bulk_data_ref$cells, bulk_data_mapped$cells), 
                 counts = count_mat_comb)


###########################################################
# remove clusters with each <2 pseudobulks per cluster per cluster_group (and clusters not in bulk_data_ref)
###########################################################

message("\n\n          *** Preparing dataset for DESeq2 analysis... ", Sys.time(),"\n\n")

t1 = bulk_data$meta[bulk_data$meta$sample_type != "in_vitro",]
t2 = t1 %>% group_by(cluster_name, sample_type, group, comp, comp_group) %>% summarise(N_pseudobulks = n())

comp_clusters = NULL

for (cl in unique(bulk_data_ref$meta$cluster_name)){
  t3 = t2[t2$cluster_name == cl,]
  if (nrow(t3) == length(unique(t2$sample_type))*length(unique(t2$group))  ){
    if(min(t3$N_pseudobulks >=2)){
      comp_clusters = c(comp_clusters, cl)
    }
  }
}

t1 = bulk_data$meta
comp_clusters = c(comp_clusters, unique(t1$cluster_name[t1$sample_type == "in_vitro"]))

comp_meta = t1[t1$cluster_name %in% comp_clusters,]

rownames(comp_meta) = comp_meta$cluster_sample

bulk_data$meta = comp_meta


#remove genes from countmatrix with <0.1 counts/cell in all pseudobulks  

m1 = bulk_data$counts

pseudobulks = na.omit(comp_meta$pseudobulk)

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
                             design = ~comp_group)
dds = DESeq(dds) #standard Wald test

bulk_data$deseq_dataset = dds

#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 



#################################################
# Extract DESeq2 results
#################################################

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

gr = unique(comp_meta$group)

#extract stats for each comparison (genes signif dereg and with tendency of dereg)
deseq_res_list = NULL
genes_reg_list = NULL

for (comp in unique(comp_meta$comp)){
  message("  * extracting DESeq results for ", comp)
  
  comp_gr2 = paste0(comp, "_", gr[2])
  comp_gr1 = paste0(comp, "_", gr[1])
  t1 = as.data.frame(results(dds, contrast=c("comp_group", comp_gr2, comp_gr1)))
  
  deseq_res_list[[paste0(comp, "_", gr[2], "_vs_", gr[1])]] = t1
  
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.1 & t1$log2FoldChange >= +log2(1.2)]
  genes_reg_list[[paste0(comp, "_", gr[2],"_up")]] = v3
  v3 = rownames(t1)[!is.na(t1$padj) & 
                      t1$padj <= 0.1 & t1$log2FoldChange <= -log2(1.2)]
  genes_reg_list[[paste0(comp, "_", gr[2],"_down")]] = v3
  
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
        "# Completed J02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


