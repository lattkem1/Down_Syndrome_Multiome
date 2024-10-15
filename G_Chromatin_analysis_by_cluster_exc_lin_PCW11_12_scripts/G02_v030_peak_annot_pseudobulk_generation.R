message("\n\n##########################################################################\n",
        "# Start G02: Peak annotation and pseudobulking chromatin accessibility read counts  ", Sys.time(),
        "\n##########################################################################\n",
        "\n    with subsampling of Chr21 counts due remove bias by 3rd copy of Chr21 ",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MatrixGenerics)


#specify script/output index as prefix for file names
script_ind = "G02_"

#specify output directory
out_dir = paste0(main_dir,"G_Chromatin_analysis_by_cluster_exc_lin_PCW11_12/")

#load group and file info
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")
gr_tab = gr_tab[gr_tab$dev_PCW <= 12,]

#load dataset

clust_tab = read_csv(paste0(main_dir,"C_subsetting_exc_lin_PCW11_12/C02_subset_cluster_assignment.csv"))

load(file = paste0(out_dir, "G01_seur_w_peaks_by_cluster_quant.rda"))


###################################################
# annotate peaks create updated tsv and bed file of peakset 
#       (add peak IDs, map peaks to genomic elements and nearest promoters)
###################################################

message("\n\n          *** Annotating peaks... ", Sys.time(),"\n\n")

t1 = seur@assays$peaks_by_cluster
t2 = as.data.frame(peaks)
t2$peak_id = rownames(t1)

t3 = t2[,c("seqnames", "start",    "end",    "peak_id",  "width" ,   "strand", "peak_called_in" )]

peaks_by_cluster = t3

### extract GRanges object for peaks 

t1 = peaks_by_cluster

peaks_GRanges = makeGRangesFromDataFrame(t1)


### extract gene, promoter, exon, intron annotations from Ensembl database 

edb = EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"

ensembl_genes = genes(edb)
t2 = findOverlaps(peaks_GRanges, ensembl_genes)
overlaps = as_tibble(t2)
t3 = as_tibble(ensembl_genes)
peaks_gene_annot = peaks_by_cluster[overlaps$queryHits, ]
peaks_gene_annot$gene_id = t3$gene_id[overlaps$subjectHits]
peaks_gene_annot$gene_name = t3$gene_name[overlaps$subjectHits]
peaks_gene_annot$gene_biotype = t3$gene_biotype[overlaps$subjectHits]

ensembl_promoters = promoters(edb, columns = c("gene_id", "gene_name" ))
t2 = distanceToNearest(peaks_GRanges, ensembl_promoters)
overlaps = as_tibble(t2)
t3 = as_tibble(ensembl_promoters)
peaks_prom_annot = peaks_by_cluster[overlaps$queryHits, ]
peaks_prom_annot$gene_id = t3$gene_id[overlaps$subjectHits]
peaks_prom_annot$closest_promoter = t3$gene_name[overlaps$subjectHits]
peaks_prom_annot$dist_to_closest_promoter = overlaps$distance

ensembl_introns = intronicParts(edb)
t2 = findOverlaps(peaks_GRanges, ensembl_introns)
overlaps = as_tibble(t2)
t3 = as_tibble(ensembl_introns)
peaks_int_annot = peaks_by_cluster[overlaps$queryHits, ]
peaks_int_annot$gene_id = t3$gene_id[overlaps$subjectHits]

ensembl_exons = exonicParts(edb)
t2 = findOverlaps(peaks_GRanges, ensembl_exons)
overlaps = as_tibble(t2)
t3 = as_tibble(ensembl_exons)
peaks_exon_annot = peaks_by_cluster[overlaps$queryHits, ]
peaks_exon_annot$gene_id = t3$gene_id[overlaps$subjectHits]


### add Ensembl annotations to peak set
#     first add closest promoter annotation and distance to promoter
#     then add gene annotation, then intron, then exon (=> if overlapping with intron and exon, assign as potential exon)
#     assign peaks not annotated for gene as promoter (dist_to_closest_promoter = 0) or intergenic (dist_to_closest_promoter>0)

t1 = peaks_prom_annot
peaks_by_cluster$closest_promoter = t1$closest_promoter[match(peaks_by_cluster$peak_id, t1$peak_id)]
peaks_by_cluster$dist_to_closest_promoter = t1$dist_to_closest_promoter[match(peaks_by_cluster$peak_id, t1$peak_id)]

t1 = peaks_gene_annot
peaks_by_cluster$gene_name = t1$gene_name[match(peaks_by_cluster$peak_id, t1$peak_id)]
peaks_by_cluster$gene_biotype = t1$gene_biotype[match(peaks_by_cluster$peak_id, t1$peak_id)]

peaks_by_cluster$peak_type[peaks_by_cluster$peak_id %in% peaks_int_annot$peak_id] = "intron"
peaks_by_cluster$peak_type[peaks_by_cluster$peak_id %in% peaks_exon_annot$peak_id] = "exon"
peaks_by_cluster$peak_type[peaks_by_cluster$dist_to_closest_promoter == 0] = "promoter"
peaks_by_cluster$peak_type[is.na(peaks_by_cluster$peak_type) & 
                             peaks_by_cluster$dist_to_closest_promoter > 0] = "intergenic"

#save as csv and bed file

write_csv(peaks_by_cluster, file = paste0(out_dir,script_ind, "peaks_by_cluster_annotated.csv"))

write_tsv(peaks_by_cluster[,1:6], col_names = FALSE, file = paste0(out_dir, script_ind, "peaks_by_cluster.bed"))



#calculate number of peaks by peak types

peak_stats = peaks_by_cluster %>% group_by(peak_type) %>% summarise(N_peaks_all = n())

write_csv(peak_stats, file = paste0(out_dir,script_ind, "peak_stats.csv"))





###########################################################
# define grouping and pseudobulking variables (pseudobulk by cluster_sample)
###########################################################

message("\n\n          *** Creating pseudobulk dataset... ", Sys.time(),"\n\n")

#define grouping variables (clusters ordered by number)
gr = unique(gr_tab$group)
samples = as.character(unique(gr_tab$sample))
cell_classes = unique(clust_tab$cell_class)
cluster_names = intersect(clust_tab$cluster_name, seur$cluster_name)
stages = unique(seur$stage[order(seur$stage)] )


#define combined variable for clusters by sample
seur$cluster_sample = paste0(seur$cluster_name, "_", seur$sample)
Idents(seur) = "cluster_sample"

cluster_samples = unique(seur$cluster_sample[order(match(seur$cluster_name, cluster_names), 
                                                   match(seur$sample, samples))])

#define combined variable for clusters by stage
seur$cluster_stage = paste0(seur$cluster_name, "_", seur$stage)


###########################################################
# pseudobulking by cluster and sample combination (cluster_sample)
###########################################################

#identify cells for each cluster_sample (keep only bulks with >20 cells/pseudobulk), create bulk_metadata

t1 = seur@meta.data
t2 = t1 %>% group_by(cell_class, cluster_name, sample, stage, cluster_stage, cluster_sample) %>% summarise(N_cells = n())
t2$pseudobulk[t2$N_cells>=10] = t2$cluster_sample[t2$N_cells>=10]
t3 = gr_tab[match(t2$sample, gr_tab$sample), colnames(gr_tab) != "sample"]
t4 = cbind(t2, t3)
t4 = t4[order(match(t2$cluster_sample, cluster_samples)),]

#define scaled stage variable 
t4$dev_PCW_scaled = as.vector(scale(t4$dev_PCW))

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
sc_counts = seur[["peaks_by_cluster"]]$counts

for (c1 in 1:length(bulk_cell_list)){
  bulk_cells = bulk_cell_list[[c1]]
  bulk_count_list[[names(bulk_cell_list)[c1] ]] =apply(sc_counts[,bulk_cells], 1, sum)
  if(c1%%10 == 0){message("generated pseudobulk ", c1, " of ", length(bulk_cell_list))}
}

bulk_mat = as.matrix(as.data.frame(bulk_count_list))

bulk_meta = bulk_meta[!is.na(bulk_meta$pseudobulk),]

# collect data in bulk_data object
bulk_data_atac = list(meta = bulk_meta, cells = bulk_cell_list, counts = bulk_mat)

save(bulk_data_atac, file = paste0(out_dir,script_ind, "bulk_data_atac.rda")) 



#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed G02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


