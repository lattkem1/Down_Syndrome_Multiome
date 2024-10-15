message("\n\n##########################################################################\n",
        "# Start G04: peak annotation, RNA-seq integration, visualisation: ", Sys.time(),
        "\n##########################################################################\n",
        "\n  plot N diff peaks per cluster, identify diff peaks with linked DEGs ",
        "\n  plot as heatmap and genome tracks ",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E12_240806_DS_foetal_brain_grafts_for_man_v01/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(colorRamps)
library(viridis)
library(pheatmap)
library(DESeq2)
library(DEGreport)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#specify script/output index as prefix for file names
script_ind = "G04_"

#specify output directory
out_dir = paste0(main_dir,"G_Chromatin_analysis_by_cluster_exc_lin_PCW11_12/")


#load group and file info
gr_tab = read_csv("B_basic_analysis/B02_gr_tab_filtered.csv")
gr_tab = gr_tab[gr_tab$dev_PCW <= 12,]

#load ATAC DEseq2 dataset
load(file = paste0(out_dir,"G03_bulk_data_with_DESeq_results.rda")) 
bulk_data_atac = bulk_data

#load peakset
peaks_by_cluster = read_csv(paste0(out_dir,"G02_peaks_by_cluster_annotated.csv"))


#load RNA DESeq2 dataset for linking diff peaks and genes
load(file = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_PCW11_12/E03_bulk_data_w_expr_z_scores.rda")) 
bulk_data_rna = bulk_data

#load Seurat dataset
load(paste0(out_dir,"G01_seur_w_peaks_by_cluster_quant.rda"))



#get marker gene panels

GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol



gc()

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


### function: calculate and plot mean expr Z-score for selected metadata variable (mean_z_scores_by) 
#         for each CON and DS samples, and group difference

calculate_mean_Z_by = function(mat_full, meta, mean_z_scores_by){
  
  meta[["mean_var"]]= meta[[mean_z_scores_by]]
  meta2 = meta
  
  mean_var_vals = unique(meta$mean_var)
  
  
  mean_mat_CON = matrix(nrow= nrow(mat_full), ncol = length(mean_var_vals),
                        dimnames = list(rownames(mat_full), mean_var_vals))
  
  mean_mat_DS = mean_mat_CON 
  
  for (mean_var_val in mean_var_vals){
    
    samples_CON = meta$cluster_sample[meta$mean_var == mean_var_val & meta$group == "CON"]
    m1 = mat_full[,samples_CON]
    mean_mat_CON[,mean_var_val] = apply(m1, 1, mean)
    
    samples_DS = meta$cluster_sample[meta$mean_var == mean_var_val & meta$group == "DS"]
    m1 = mat_full[,samples_DS]
    mean_mat_DS[,mean_var_val] = apply(m1, 1, mean)
    
    #keep only meta columns with unique values for each value of the mean_var
    v1 = unlist(lapply(meta2[meta2$mean_var == mean_var_val,], function(x){length(unique(x))}))
    meta2 = meta2[, intersect(colnames(meta2), names(v1)[v1 == 1])]
    
  }
  meta2 = meta2[!duplicated(meta2$mean_var),]
  rownames(meta2) = meta2$mean_var
  
  mean_mat_list = list()
  mean_mat_list$meta = meta2
  mean_mat_list$CON = mean_mat_CON
  mean_mat_list$DS = mean_mat_DS
  mean_mat_list$DELTA = mean_mat_DS - mean_mat_CON
  
  return(mean_mat_list)
  
}




###function: plot bulk gene expression heatmap with annotations

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = NULL,  cellwidth = 15, cellheight = 10, 
                            fontsize = 10, title = "Z-score vst-norm gene expression"){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}
  
  pl_mat = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(pl_meta[[x_col]], colnames(pl_mat), nomatch = 0),]
  
  #define annotation bars
  
  if (!is.null(meta_annot_cols)){
    
    annot_col = data.frame(row.names = pl_meta[[x_col]]) #if not defined cbind converts factor values to factor levels
    
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      v1 = factor(v1, levels = unique(pl_meta[[col1]]))
      annot_col = as.data.frame(cbind(annot_col, v1))
    }
    colnames(annot_col) = meta_annot_cols
    rownames(annot_col) = colnames(pl_mat)
    
    annot_colors = lapply(meta_annot_cols, function(x){
      v1 = pal(unique(annot_col[[x]]) )
      names(v1) = levels(annot_col[[x]])
      return(v1)
    })
    names(annot_colors) = meta_annot_cols
    
  }else{
    annot_col = NULL
    annot_colors = NULL
  }
  
  # define plot limits
  
  if (is.null(lims)){lims = c(-0.7*max(abs(na.omit(pl_mat))), 0.7*max(abs(na.omit(pl_mat))))}
  
  #create plot
  
  p1 = pheatmap(pl_mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                color = color,
                breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                annotation_col = annot_col, annotation_colors = annot_colors,
                border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                fontsize = fontsize, main = title
  )
  
  return(p1)
  
}



#################################################
# Plot number of diff peaks (incl Chr21 peaks) by cluster
#################################################

bulk_data = bulk_data_atac

l1 = bulk_data$DEGs

v1 = unlist(bulk_data$DEGs)
v2 = v1[grepl("chr21-", v1)]

l2 = lapply(l1, intersect, v2)

meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)

t1 = tibble(comps = names(l1), cluster_name = NA, DS_up_down = NA,
            N_peaks = lengths(l1), N_peaks_Chr21 = lengths(l2))

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

p2 = p1 + geom_col(aes(y = N_peaks, fill = DS_up_down), width = 0.7, position = position_dodge(width = 0.7))+
  labs(title = "All diff peaks")

p3 = p1 + geom_col(aes(y = N_peaks_Chr21, fill = DS_up_down), width = 0.7, position = position_dodge(width = 0.7))+
  labs(title = "Chr21 diff peaks")

p4 = p1 + geom_col(aes(y = N_peaks, color = DS_up_down), width = 0.7, position = position_dodge(width = 0.9),
                   fill = "grey80")+
  geom_col(aes(y = N_peaks_Chr21, fill = DS_up_down, color = DS_up_down), 
           width = 0.7, position = position_dodge(width = 0.9))+
  labs(title = "diff peaks Chr21 vs others")

pdf(file = paste0(out_dir,script_ind,"Diff_peaks_numbers_by_cluster.pdf"), 
    width = 5, height = 4)
{
  plot(p2)
  plot(p3)
  plot(p4)
}
dev.off()





###############################################
### identify diff peaks linked to diff genes
###############################################

diff_peaks_linked_to_diff_genes = peaks_by_cluster[(peaks_by_cluster$peak_id %in% unlist(bulk_data_atac$DEGs)) &
                                                     (peaks_by_cluster$closest_promoter %in% unlist(bulk_data_rna$DEGs)),]

write_csv(diff_peaks_linked_to_diff_genes, file = paste0(out_dir,script_ind, "diff_peaks_linked_to_diff_genes.csv"))


### stats all peaks, cleaned peaks, diff peaks, diff peaks linked to diff genes (N and fraction of peaks)

peak_stats = peaks_by_cluster %>% group_by(peak_type) %>% summarise(N_peaks_all = n())

t2 = peaks_by_cluster[peaks_by_cluster$peak_id %in% rownames(bulk_data_atac$deseq_dataset),] %>% 
  group_by(peak_type) %>% summarise(N_peaks_cleaned = n())
peak_stats$N_peaks_cleaned = t2$N_peaks_cleaned[match(t2$peak_type, peak_stats$peak_type)]

t2 = peaks_by_cluster[peaks_by_cluster$peak_id %in% unlist(bulk_data_atac$DEGs),] %>% 
  group_by(peak_type) %>% summarise(N_peaks_diff = n())
peak_stats$N_peaks_diff = t2$N_peaks_diff[match(t2$peak_type, peak_stats$peak_type)]

t2 = diff_peaks_linked_to_diff_genes %>% group_by(peak_type) %>% summarise(N_peaks_diff_linked_diff_genes = n())
peak_stats$N_peaks_diff_linked_diff_genes = t2$N_peaks_diff_linked_diff_genes[match(t2$peak_type, peak_stats$peak_type)]

t2 = peak_stats[,-1]
t3 = as.data.frame(lapply(t2, function(x){x/sum(x)}))
names(t3) = paste("fract_", names(t3))
peak_stats = cbind(peak_stats, t3)

write_csv(peak_stats, file = paste0(out_dir,script_ind, "peak_stats.csv"))


### diff peaks vs diff peaks linked to diff genes, diff genes vs diff_genes with linked diff peaks

t1 = tibble(DEGs = length(unlist(bulk_data_rna$DEGs)), 
            DEGs_w_linked_diff_peaks = length(unique(diff_peaks_linked_to_diff_genes$closest_promoter)),
            diff_peaks = length(unlist(bulk_data_atac$DEGs)),
            diff_peaks_w_linked_DEGs = length(unique(diff_peaks_linked_to_diff_genes$peak_id)))

write_csv(t1, file = paste0(out_dir,script_ind, "diff_peaks_vs_DEGs_stats.csv"))


### add annotations to bulk_data_atac

bulk_data_atac$peaks_by_cluster_annot = peaks_by_cluster
bulk_data_atac$peaks_diff_links_diff_genes = diff_peaks_linked_to_diff_genes

save(bulk_data_atac, bulk_data_rna, file = paste0(out_dir,script_ind, "bulk_data_with_DESeq_results_atac_rna.rda")) 


#######################################
# calculate per peak Z-scores, mean Z-scores per cluster_name per group and group difference (DELTA), add to bulk_data
#######################################

message("\n\n          *** Calculate peak accessibility Z-scores... ", Sys.time(),"\n\n")

bulk_data = bulk_data_atac

meta = bulk_data$meta

vst_mat = assay(vst(bulk_data$deseq_dataset))

bulk_data$vst_mat = vst_mat


#calculate Z-score per gene by pseudobulk (cluster_sample)
cluster_sample_mat = t(apply(vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat)

bulk_data$peak_Z_scores$by_cluster_sample = cluster_sample_mat


#calculate mean Z-scores by cluster_name

l1 = calculate_mean_Z_by(mat_full = bulk_data$peak_Z_scores$by_cluster_sample, 
                         meta = meta, 
                         mean_z_scores_by = "cluster_name")
bulk_data$peak_Z_scores$mean_by_cluster_name = l1

bulk_data_atac = bulk_data

save(bulk_data_atac, file = paste0(out_dir, script_ind, "bulk_data_atac_w_access_z_scores.rda"))






#############################################
#plot heatmaps diff peaks linked to DEGs and linked DEGs (means by cluster_name)
#############################################

# extract matrix for diff peaks linked to diff genes, add linked gene to label

bulk_data = bulk_data_atac

peak_gene_tab = bulk_data$peaks_by_cluster_annot[(bulk_data$peaks_by_cluster_annot$peak_id %in% bulk_data$peaks_diff_links_diff_genes$peak_id),]
peak_gene_tab$peak_labels = paste0(peak_gene_tab$peak_id, " ( > ", peak_gene_tab$closest_promoter,")")
peak_gene_tab$gene_labels = paste0(peak_gene_tab$closest_promoter, " (", peak_gene_tab$peak_id,")")

CON_mat_peaks = bulk_data$peak_Z_scores$mean_by_cluster_name$CON[peak_gene_tab$peak_id,]
rownames(CON_mat_peaks) = peak_gene_tab$peak_labels
DELTA_mat_peaks = bulk_data$peak_Z_scores$mean_by_cluster_name$DELTA[peak_gene_tab$peak_id,]
rownames(DELTA_mat_peaks) = peak_gene_tab$peak_labels

# extract matrix for diff genes linked to diff peaks
m1 = bulk_data_rna$gene_Z_scores$mean_by_cluster_name$CON
CON_mat_genes = m1[match(peak_gene_tab$closest_promoter, rownames(m1), nomatch = 0),]
rownames(CON_mat_genes) = peak_gene_tab$gene_labels

m1 = bulk_data_rna$gene_Z_scores$mean_by_cluster_name$DELTA
DELTA_mat_genes = m1[match(peak_gene_tab$closest_promoter, rownames(m1), nomatch = 0),]
rownames(DELTA_mat_genes) = peak_gene_tab$gene_labels


#plot heatmaps (cluster by gene expression DELTA)

pdf(file = paste0(out_dir,script_ind, "Peak_access_gene_expr_heatmap_diff_peaks_w_diff_genes.pdf"), 
    width = 10, height = 30)
{
  
  p1 = bulkdata_heatmap(pl_mat = DELTA_mat_genes, 
                        pl_meta = bulk_data_rna$gene_Z_scores$mean_by_cluster_name$meta,
                        pl_genes = NULL,
                        x_col = "cluster_name", 
                        meta_annot_cols = c("cluster_name"),
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                        lims = NULL,  cellwidth = 10, cellheight = 10, 
                        fontsize = 10, title = paste0("Gene expression - Z-score mean DELTA"))
  
  gene_labels_ordered = rownames(DELTA_mat_genes)[p1$tree_row$order] 
  peak_labels_ordered = peak_gene_tab$peak_labels[match(gene_labels_ordered, peak_gene_tab$gene_labels)]
  
  bulkdata_heatmap(pl_mat = CON_mat_peaks, 
                   pl_meta = bulk_data$peak_Z_scores$mean_by_cluster_name$meta,
                   pl_genes = peak_labels_ordered,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = viridis(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, title = paste0("Peak accessibility - Z-score mean CON"))
  
  bulkdata_heatmap(pl_mat = DELTA_mat_peaks, 
                   pl_meta = bulk_data$peak_Z_scores$mean_by_cluster_name$meta,
                   pl_genes = peak_labels_ordered,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, title = paste0("Peak accessibility - Z-score mean DELTA"))
  
  bulkdata_heatmap(pl_mat = CON_mat_genes, 
                   pl_meta = bulk_data_rna$gene_Z_scores$mean_by_cluster_name$meta,
                   pl_genes = gene_labels_ordered,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = viridis(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, title = paste0("Gene expression - Z-score mean CON"))
  
  bulkdata_heatmap(pl_mat = DELTA_mat_genes, 
                   pl_meta = bulk_data_rna$gene_Z_scores$mean_by_cluster_name$meta,
                   pl_genes = gene_labels_ordered,
                   x_col = "cluster_name", 
                   meta_annot_cols = c("cluster_name"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                   lims = NULL,  cellwidth = 10, cellheight = 10, 
                   fontsize = 10, title = paste0("Gene expression - Z-score mean DELTA"))
  
}
dev.off()




#############################################
#plot peaks linked to DEGs (manually adaptable plot)
#############################################

DefaultAssay(seur) = "peaks_by_cluster"

#select regions/genes for plotting (DEGs with linked diff accessible regions)

pl_genes = unique(bulk_data_atac$peaks_diff_links_diff_genes$closest_promoter)


#select clusters vs groups and order to plot 
# subset seurat to plot only clusters analysed in pseudobulk analyses
seur$cluster_sample = paste0(seur$cluster_name, "_", seur$sample)
seur = subset(x = seur, subset = cluster_sample %in% bulk_data$meta$cluster_sample)

gc()

save(seur, file = paste0(out_dir, script_ind, "seur_subsetted_for_tests.rda"))

#set cluster x group as identities for plotting, color scheme for plotting 
seur$cluster_group = paste0(seur$cluster_name,"_",seur$group)

t1 = bulk_data$meta
cluster_groups = unique(t1$cluster_group)

Idents(seur) <- factor(x = seur$cluster_group, levels = cluster_groups)

v1 = factor(t1$group[!duplicated(t1$cluster_group)])
cluster_group_colors = pal(levels(v1))[v1]


# combined coverage plot manual; define plot groups as Idents(), specify y_order as idents in coverageplot

pl = lapply(pl_genes, function(g1){
  
  #define region around gene to include all linked peaks; peaks to plot 
  t2 = bulk_data$peaks_by_cluster_annot[bulk_data$peaks_by_cluster_annot$closest_promoter == g1,]
  pl_region = paste(t2$seqnames[1],min(t2$start), max(t2$end), sep = "-")
  pl_peaks_cleaned = makeGRangesFromDataFrame(t2, keep.extra.columns = TRUE)
  pl_peaks_diff = makeGRangesFromDataFrame(t2[t2$peak_id %in% 
                                                bulk_data$peaks_diff_links_diff_genes$peak_id,],
                                           keep.extra.columns = TRUE)
  
  cov_plot <- CoveragePlot(
    object = seur,
    region = pl_region,
    assay = "peaks_by_cluster",
    idents = cluster_groups,
    window = 500,
    extend.upstream = 5000,
    extend.downstream = 5000,
    ymax = "q90", #scale cutoff at 90 percent quantile
    annotation = FALSE,
    peaks = FALSE,
    links = FALSE
  )
  cov_plot = cov_plot + scale_fill_manual(limits = cluster_groups, values = cluster_group_colors)
  
  expr_plot <- ExpressionPlot(
    object = seur,
    features = g1,
    assay = "SCT"
  )+ scale_fill_manual(limits = cluster_groups, values = cluster_group_colors)
  
  peak_plot_all = PeakPlot(
    object = seur,
    region = pl_region,
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  
  peak_plot_cleaned = PeakPlot(
    object = seur,
    region = pl_region,
    peaks = pl_peaks_cleaned,
    color = "blue",
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  
  peak_plot_diff = PeakPlot(
    object = seur,
    region = pl_region,
    peaks = pl_peaks_diff,
    color = "red",
    extend.upstream = 5000,
    extend.downstream = 5000
  ) +geom_text(aes(x = pl_peaks_diff@ranges@start,y = -0.01, label = pl_peaks_diff$peak_id), check_overlap = TRUE)
  
  gene_plot = AnnotationPlot(
    object = seur,
    region = pl_region,
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  p1 = CombineTracks(
    plotlist = list(cov_plot, peak_plot_all, peak_plot_cleaned, peak_plot_diff, gene_plot),
    expression.plot = expr_plot,
    heights = c(15, 0.5,0.5,0.5, 1.5),
    widths = c(10, 1)
  )
  
  message("Plot generated for ", g1)
  return(p1)
})

pdf(file = paste0(out_dir,script_ind, "Genome_tracks_DEGs_w_linked_diff_peaks.pdf"), width = 8, height = 12)
lapply(pl, function(x){x})
dev.off()






#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed G04 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")



