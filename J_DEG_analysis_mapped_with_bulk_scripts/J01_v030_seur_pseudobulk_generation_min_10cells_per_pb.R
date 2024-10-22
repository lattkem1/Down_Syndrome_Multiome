message("\n\n##########################################################################\n",
        "# Start J01: Pseudobulk generation for mapped dataset ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
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


#specify script/output index as prefix for file names
script_ind = "J01_"

#specify output directory
out_dir = paste0(main_dir,"J_expression_analysis_mapped_data_vs_reference_with_neurons_in_vitro_bulk/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#load group and file info for data to mapped on reference

gr_tab = read_csv(paste0(main_dir, "A_input/group_tab_grafts.csv"))

#load mapped own dataset

load(file = paste0(main_dir,"I_mapping_to_reference_grafts_to_tissue/I02_seur_filtered_mapped.rda")) 

# load reference dataset and cell type annotations

clust_tab_ref = read_csv(paste0(main_dir, "B_basic_analysis/B03_cluster_assignment_all.csv"))

gr_tab_ref = read_csv(file =paste0(main_dir, "B_basic_analysis/B02_gr_tab_filtered.csv"))

load(file = paste0(main_dir,"E_DESeq_pseudobulk_all_by_cluster/E03_bulk_data_w_expr_z_scores.rda")) 
bulk_data_ref = bulk_data


###########################################################
# define grouping and pseudobulking variables (pseudobulk by cluster_sample)
###########################################################

set.seed(123)

#define grouping variables (clusters ordered by number)
gr = unique(gr_tab$group)
samples = as.character(unique(gr_tab$sample))
cluster_names = clust_tab_ref$cluster_name


#assign cluster_name from mapping to reference
seur$cluster_name = seur$predicted.cluster_name

#define combined variable for clusters by sample
seur$cluster_sample = paste0(seur$cluster_name, "_", seur$sample)
Idents(seur) = "cluster_sample"

cluster_samples = unique(seur$cluster_sample[order(match(seur$cluster_name, cluster_names), 
                                            match(seur$sample, samples))])



###########################################################
# pseudobulking by cluster and sample combination (cluster_sample)
###########################################################

#identify cells for each cluster_sample (keep only bulks with >20 cells/pseudobulk), create bulk_metadata

t1 = seur@meta.data
t2 = t1 %>% group_by(cluster_name, sample, cluster_sample) %>% summarise(N_cells = n())
t2$pseudobulk[t2$N_cells>=10] = t2$cluster_sample[t2$N_cells>=10]
t3 = gr_tab[match(t2$sample, gr_tab$sample), colnames(gr_tab) != "sample"]
t4 = cbind(t2, t3)
t4 = t4[order(match(t2$cluster_sample, cluster_samples)),]

bulk_meta = t4

write_csv(bulk_meta, file = paste0(out_dir, script_ind, "bulk_meta.csv"))

bulk_cell_list = lapply(cluster_samples, function(s1){
  cells =  rownames(t1)[t1$cluster_sample == s1]
  return(cells)
})
names(bulk_cell_list) = cluster_samples

bulk_cell_list = bulk_cell_list[lengths(bulk_cell_list)>=10]


# create list of pseudobulk counts and convert to pseudobulk matrix

bulk_count_list = NULL
sc_counts = seur[["RNA"]]$counts

for (c1 in 1:length(bulk_cell_list)){
  bulk_cells = bulk_cell_list[[c1]]
  bulk_count_list[[names(bulk_cell_list)[c1] ]] =apply(sc_counts[,bulk_cells], 1, sum)
  if(c1%%10 == 0){message("generated pseudobulk ", c1, " of ", length(bulk_cell_list))}
}

bulk_mat = as.matrix(as.data.frame(bulk_count_list))

bulk_meta = bulk_meta[!is.na(bulk_meta$pseudobulk),]


# collect data in bulk_data object
bulk_data = list(meta = bulk_meta, cells = bulk_cell_list, counts = bulk_mat)


### define grouping variables and orders

t2 = list()
t2$group = unique(bulk_meta$group)
t2$sample = unique(gr_tab$sample)
t2$cluster_name = cluster_names
t2$cluster_sample = unique(bulk_meta$cluster_sample)

bulk_data$grouping = t2

bulk_data_mapped = bulk_data

save(bulk_data_mapped, file = paste0(out_dir,script_ind, "bulk_data_mapped_data.rda")) 





message("\n\n##########################################################################\n",
        "# Completed J01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


