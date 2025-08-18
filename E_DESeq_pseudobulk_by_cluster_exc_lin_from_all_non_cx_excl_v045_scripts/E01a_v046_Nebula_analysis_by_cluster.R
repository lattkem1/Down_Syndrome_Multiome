message("\n\n##########################################################################\n",
        "# Start E01a: DEG analysis with Nebula   ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
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
library(nebula)

#specify script/output index as prefix for file names
script_ind = "E01a_"


#specify output directory
out_dir = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_v045/")


#load group and file info
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered_non_cx_excl.csv")

#load dataset

load(file = paste0(main_dir,"C_subsetting_exc_lin_from_all_non_cx_excl/C03_seur_integr_labelled.rda")) 

clust_tab = read_csv(paste0(main_dir,"C_subsetting_exc_lin_from_all_non_cx_excl/C02_subset_cluster_assignment.csv"))


### get marker gene panels

GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol


###########################################################
# create metadata table for cluster_sample (cells grouped by cluster and sample); 
# calculate mean SCT per cluster_sample for visualisations
###########################################################

message("\n\n          *** Extract cluster_sample metadata, calculate mean SCT per cluster_sample  - ", Sys.time(),"\n\n")

#define grouping variables (clusters ordered by number)
gr = unique(gr_tab$group)
samples = as.character(unique(gr_tab$sample))
cell_types = unique(clust_tab$cell_type)
cluster_names = clust_tab$cluster_name


#define combined variable for clusters by sample
seur$cluster_sample = paste0(seur$cluster_name, "_", seur$sample)
Idents(seur) = "cluster_sample"

cluster_samples = unique(seur$cluster_sample[order(match(seur$cluster_name, cluster_names), 
                                                   match(seur$sample, samples))])

#define metadata

t1 = seur@meta.data
t2 = t1 %>% group_by(cell_type, cluster_name, sample, cluster_sample) %>% summarise(N_cells = n())
t2 = t2[t2$cluster_name %in% clust_tab$cluster_name,]
t3 = gr_tab[match(t2$sample, gr_tab$sample), colnames(gr_tab) != "sample"]
t4 = cbind(t2, t3)
t4 = t4[order(match(t2$cluster_sample, cluster_samples)),]
rownames(t4) = t4$cluster_sample

cluster_sample_meta = t4


###########################################################
# Nebula analysis
###########################################################

options(future.globals.maxSize = 30 * 1024^3) #SCTransform exceeds default memory limit for parallelisation with futures

meta = seur@meta.data

nebula_data = list(meta = cluster_sample_meta, gr_tab = gr_tab, clust_tab = clust_tab)


### analysis by cluster

cluster_names = clust_tab$cluster_name

for (cl in cluster_names){
  
  message("\n          *** Nebula analysis cluster ", cl, " - ", Sys.time(),"\n")
  
  seuratdata <- scToNeb(obj = seur[, seur$cluster_name == cl], 
                        assay = "RNA", id = "sample", pred = c("group"), offset="nCount_RNA")
  
  df = model.matrix(~group, data=seuratdata$pred)
  
  re = nebula(seuratdata$count,seuratdata$id,pred=df,offset=seuratdata$offset, ncore = 16)
  
  nebula_data$results[[cl]]$df = df
  nebula_data$results[[cl]]$res = re
  
}

message("\n\n          *** Save Nebula analysis - ", Sys.time(),"\n\n")

save(nebula_data, file = paste0(out_dir,script_ind, "nebula_data.rda")) 



#################################################
# Extract Nebula results
#################################################

message("\n\n          ***Extracting Nebula results... ", Sys.time(),"\n\n")

# load(file = paste0(out_dir,script_ind, "nebula_data.rda")) 

for (comp in names(nebula_data$results)){
  
  l1 = nebula_data$results[[comp]]
  t1 = l1$res$summary
  
  t2 = t1[t1$p_groupDS <= 0.05,]
  
  nebula_data$DEGs[[paste0(comp,"_DS_up")]] = t2$gene[t2$logFC_groupDS >= log2(1.2)]
  nebula_data$DEGs[[paste0(comp,"_DS_down")]] = t2$gene[t2$logFC_groupDS <= -log2(1.2)]
  
}



### calculate mean SCT norm expression per cluster_sample

sct_mat = seur@assays$SCT@counts

#remove lowly expressed genes
v1 = rowSums(sct_mat)
sct_mat = sct_mat[v1/ncol(sct_mat)>=0.001 | rownames(sct_mat) %in% unlist(nebula_data$DEGs),]

cluster_sample_sct_mat = matrix(nrow = nrow(sct_mat), ncol = nrow(cluster_sample_meta), 
                                dimnames = list(rownames(sct_mat), rownames(cluster_sample_meta)))

cluster_samples = colnames(cluster_sample_sct_mat)

for (i in 1:length(cluster_samples)){
  
  cs = cluster_samples[i]
  
  m1 = sct_mat[,seur$cluster_sample == cs]
  
  cluster_sample_sct_mat[,cs] = apply(m1, 1, mean)
  
  if(i%%50 == 0){message("calculated mean SCT for cluster_sample ", i, " of ", length(cluster_samples))}
  
}

nebula_data$sct_mat = cluster_sample_sct_mat


### save bulk_dataset with DESeq results

save(nebula_data, file = paste0(out_dir,script_ind, "nebula_data.rda")) 



### save table with DEGs by comparison including TFs and HSA21 genes

l1 = nebula_data$DEGs
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

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "Nebula_up_down_by_comp_genes.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "Nebula_up_down_by_comp_N_genes.csv"))



#################################################
# Plot number of DEG (incl Chr21 genes) by comp_group
#################################################

comp_groups = nebula_data$clust_tab$cluster_name

l1 = nebula_data$DEGs
l2 = lapply(l1, intersect, GOI$Chr21)
meta = nebula_data$meta

t1 = tibble(comps = names(l1), comp_group = NA, DS_up_down = NA,
            N_genes = lengths(l1), N_genes_Chr21 = lengths(l2))

for (comp_group in comp_groups){
  t1$comp_group[grepl(comp_group, t1$comps)] = comp_group
}

t1$DS_up_down[grepl("up", t1$comps)] = "up"
t1$DS_up_down[grepl("down", t1$comps)] = "down"

p1 = ggplot(t1, aes(x = comp_group, group = DS_up_down))+
  scale_x_discrete(limits = comp_groups)+
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
        "# Completed E01a ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


