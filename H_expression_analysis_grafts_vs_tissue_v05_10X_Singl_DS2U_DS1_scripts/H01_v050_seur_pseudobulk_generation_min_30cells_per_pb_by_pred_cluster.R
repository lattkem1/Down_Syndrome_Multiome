message("\n\n##########################################################################\n",
        "# Start H01a: Pseudobulk generation by reference cluster mapping (grafts)  ", Sys.time(),
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


#specify script/output index as prefix for file names
script_ind = "H01_"

#specify output directory
out_dir = paste0(main_dir,"H_expression_analysis_grafts_vs_tissue_v05_10X_Singl_DS2U_DS1/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}


#load group and file info
gr_tab = read_csv("G_basic_analysis_grafts_map_to_tissue_v05_10X_Singl_DS2U_DS1/G02_gr_tab_filtered.csv")

#load dataset and reference clusters

clust_tab =read_csv("C_subsetting_exc_lin_from_all_non_cx_excl/C02_subset_cluster_assignment.csv")

load(file = "G_basic_analysis_grafts_map_to_tissue_v05_10X_Singl_DS2U_DS1/G04_seur_integr_labelled.rda") 



###########################################################
# define grouping and pseudobulking variables (pseudobulk by cluster_sample)
###########################################################

set.seed(123)

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



###########################################################
# pseudobulking by cluster and sample combination (cluster_sample)
###########################################################

#identify cells for each cluster_sample (keep only bulks with >10 cells/pseudobulk), create bulk_metadata

t1 = seur@meta.data
t2 = t1 %>% group_by(cell_type, cluster_name, sample, cluster_sample) %>% summarise(N_cells = n())
t2 = t2[t2$cluster_name %in% clust_tab$cluster_name,]
t2$pseudobulk[t2$N_cells>=30] = t2$cluster_sample[t2$N_cells>=30]
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

bulk_cell_list = bulk_cell_list[lengths(bulk_cell_list)>=30]


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
bulk_data = list(meta = bulk_meta, gr_tab = gr_tab, clust_tab = clust_tab, 
                 cells = bulk_cell_list, counts = bulk_mat)


save(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.rda")) 



message("\n\n##########################################################################\n",
        "# Completed H01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


