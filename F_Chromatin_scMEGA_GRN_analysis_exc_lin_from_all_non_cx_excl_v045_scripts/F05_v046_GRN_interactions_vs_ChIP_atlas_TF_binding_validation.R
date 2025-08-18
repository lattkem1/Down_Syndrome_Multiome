message("\n\n##########################################################################\n",
        "# Start F05: GRN interactions vs published ChIP peaks: Load peak sets ", Sys.time(),
        "\n##########################################################################\n",
        "\n create/update locally saved collection of bed files from ChIP Atlas database ",
        "\n check locally saved collection, download missing ChIP sets for GRN TFs",
        "\n##########################################################################\n\n")

main_dir = paste0("/rds/general/user/mlattke/projects/dsseq23/live/",
                  "E14_241219_DS_foetal_brain_grafts_for_man_v02_low_string/")
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(RCurl)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(colorRamps)
library(viridis)
library(pheatmap)
library(ggraph)
library(tidygraph)

#specify script/output index as prefix for file names
script_ind = "F05_"

#specify output directory
out_dir = paste0(main_dir,"F_Chromatin_scMEGA_GRN_analysis_exc_lin_from_all_non_cx_excl_v045/")

#load Seurat dataset with GRN
load(file = paste0(out_dir,"F03_seur_with_GRN.rda")) 

#load dataset pseudobulk dataset
load(file = paste0(out_dir,"F04_bulk_data_with GRN_mean_vst.rda")) 

#load annotated peak set and network nodes
peaks = read_csv(file = paste0(out_dir,"F02_peaks_annotated.csv"))

net_nodes = bulk_data$GRN_DS_reg$nodes
net_TFs = unique(net_nodes$gene[net_nodes$N_targets>0])

#load previously downloaded ChIP peak sets

ChIP_atlas_dir = "/rds/general/user/mlattke/projects/dsseq23/live/ChIP_Atlas_peak_beds/"

if (dir.exists(ChIP_atlas_dir)){
  peaks_files_local = list.files(ChIP_atlas_dir)
  TFs_saved_local = str_remove_all(peaks_files_local, ".csv")
  TFs_to_load = net_TFs[!(net_TFs %in% TFs_saved_local)]
} else {
  dir.create(ChIP_atlas_dir)  
  TFs_to_load = net_TFs
  }

net_TFs = intersect(net_TFs, c(TFs_saved_local, TFs_to_load))
net_TFs = net_TFs[!(net_TFs %in% c("CTCF"))]


#get marker gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol
t1 = read_tsv(paste0(main_dir,"A_input/Genomics_EnglandPanelApp_Intellectual_disability_v8.243.tsv"))
GOI$ID_genes = unique(t1$`Gene Symbol`)



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



### function: plot network graph 
# default: node-shape: gene type, node-size: mean Z-score all CON clusters, node-colour: mean DELTA Z-score all clusters
# nodes/edges: full network
# plot_nodes: nodes to plot 
# node_size: name of column specifying node size
# node_color: name of column specifying node color

grn_plot = function(nodes, edges, plot_genes = NULL, layout='fr', 
                    circular = FALSE, 
                    node_size_col = NULL, node_color_col = NULL,
                    node_size_limits = NULL, node_color_limits = NULL, 
                    edge_alpha = 1, edgewidth_range = c(0.1, 1),
                    plot_unconnected = FALSE, title = "GRN plot"){
  
  if (is.null(plot_genes)){plot_genes = nodes$gene} 
  if (is.null(node_size_col)){node_size_col = "mean_z_CON"}
  if (is.null(node_color_col)){node_color_col = "mean_z_DELTA"}
  if (is.null(node_size_limits)){node_size_limits = c(0, max(nodes[[node_size_col]]))}
  if (is.null(node_color_limits)){node_color_limits = c(-max(abs(nodes[[node_color_col]])), 
                                                        +max(abs(nodes[[node_color_col]])))}
  
  plot_nodes = nodes[nodes$gene %in% plot_genes,]
  plot_nodes$node_size = plot_nodes[[node_size_col]]
  plot_nodes$node_color = plot_nodes[[node_color_col]]
  
  plot_edges = edges[edges$from %in% plot_genes & edges$to %in% plot_genes,]
  
  if (!plot_unconnected){
    plot_nodes = plot_nodes[plot_nodes$gene %in% c(plot_edges$from, plot_edges$to),]
  }
  
  g1 = tbl_graph(nodes = plot_nodes, 
                 edges = plot_edges, 
                 node_key = "gene")
  
  gene_classes = c("TF", "TF_Chr21", "Chr21", "other")
  gene_class_shapes = c(22, 23, 24, 21)
  
  set.seed(123)
  
  if (circular == TRUE){
    p1 = ggraph(g1, layout=layout, circular == TRUE)
  } else {p1 = ggraph(g1, layout=layout)}
  
  p1 = p1 +  
    labs(title = title) + 
    geom_edge_hive(aes(linetype=correlation<0, edge_width = abs(correlation)), strength = 0,
                   edge_color = "grey30",
                   arrow = arrow(angle = 10, length = unit(3, "mm"),
                                 ends = "last", type = "closed"),
                   start_cap = circle(1, 'mm'),
                   end_cap = circle(1, 'mm'), alpha = edge_alpha) +
    scale_edge_width_continuous(limits = c(0,1), range = edgewidth_range)+
    scale_edge_linetype_manual(limits = c(FALSE, TRUE), values = c("solid", "longdash"))+
    geom_node_point(aes(fill = node_color, shape = gene_class, size = node_size), stroke = 0.01) +
    scale_size_continuous(limits = node_size_limits, range = c(0.05 , 5), name = node_size_col)+
    geom_node_text(aes(label=labels), size=3, repel=T, fontface = "bold") +
    scale_fill_gradient2(limits = node_color_limits, low = "blue", mid = "white", high = "red", 
                         midpoint = 0, name = node_color_col)+
    scale_shape_manual(limits = gene_classes, values = gene_class_shapes)+
    theme_void()
  
  return(p1)
  
}





###############################################
# download missing bed files for TF ChIP from ChIP-Atlas 
################################################

message("\n\n     *** Download missing peak sets... (", Sys.time(), ") \n\n")

if(length(TFs_to_load)>0){
  
  for (i in 1:length(TFs_to_load)){
    tf = TFs_to_load[i]
    message("downloading ChIP peaks for ", tf, " (TF ", i, " of ", length(TFs_to_load), ")")
    
    file1 = paste0("https://chip-atlas.dbcls.jp/data/hg38/assembled/Oth.ALL.05.",
                   tf, ".AllCell.bed")
    
    if (url.exists(file1)){
      
      t1 = read_lines(file1)
      t2 = t1[-1]
      t2 = str_split(t2, "\t")
      ChIP_peaks = as.data.frame(do.call(rbind, t2))
      colnames(ChIP_peaks) = c("seqnames", "start", "end")
      write_csv(ChIP_peaks, file = paste0(ChIP_atlas_dir, tf, ".csv"))
    } 
  }
}


###############################################
# load  ChIP-Atlas peak sets for overlapping with network peaks/TFs
################################################

message("\n\n     *** Load peak sets... (", Sys.time(), ") \n\n")

#load names of TFs with ChIP peak sets available for analysis, keep network TFs with available peaks for analysis

peaks_files_local = list.files(ChIP_atlas_dir)
net_TFs = intersect(net_TFs, str_remove_all(peaks_files_local, ".csv"))

ChIP_peaks_list = list()

for (i in 1:length(net_TFs)){
  
  tf = net_TFs[i]
  
  message("load ChIP peaks for ", tf, " (TF ", i, " of ", length(net_TFs), ")")
  
  t1 = read_csv(file = paste0(ChIP_atlas_dir, tf, ".csv"))
  
  ChIP_peaks_list[[tf]] = t1
  
}



###############################################
# overlap targets from GRN with published targets
################################################

#get network edges and nodes
net_nodes = bulk_data$GRN_DS_reg$nodes
net_edges = bulk_data$GRN_DS_reg$edges

#select TFs for validation
net_TFs = unique(net_nodes$gene[net_nodes$N_targets>0])
net_TFs = intersect(net_TFs, names(ChIP_peaks_list))

#filter ChIP peak set for TFs in GRN analysis
ChIP_peaks_list = ChIP_peaks_list[names(ChIP_peaks_list) %in% net_TFs]


overlap_list = list()
peaks_tf_in_published_list = list()
N_overlaps_tab = NULL

for (tf in names(ChIP_peaks_list)){
  
  message("  ** Extracting overlap predicted regulatory elements - ChIP peaks for TF - ", tf)
  
  #identify overlap of published ChIP peaks for TF with ATAC peaks
  
  ChIP_peaks_GRanges = makeGRangesFromDataFrame(ChIP_peaks_list[[tf]])
  peaks_GRanges = makeGRangesFromDataFrame(peaks)
  
  t1 = findOverlaps(peaks_GRanges, ChIP_peaks_GRanges)
  overlaps = as_tibble(t1)
  peaks_ChIP_overlap = peaks[overlaps$queryHits, ]
  peaks[[paste0("ChIP_overlap_", tf)]][peaks$peak_id %in% peaks_ChIP_overlap$peak_id] = tf
  
  #quantify how many of predicted targets have regulatory elements binding to TF in published ChIP data 
  
  t1 = motif.matching[,tf]
  peaks_tf_motif = names(t1)[t1]
  
  #add information whether peaks have TF motif, show TF ChIP peak and are predicted to regulate gene 
  peaks_to_gene_w_tf_motifs = df.p2g
  peaks_to_gene_w_tf_motifs$motif[peaks_to_gene_w_tf_motifs$peak %in% peaks_tf_motif] = tf
  
  genes_tf_targets = df.grn[df.grn$tf == tf & df.grn$gene %in% peaks_to_gene_w_tf_motifs$gene &!is.na(df.grn$fdr),]
  
  peaks_to_gene_w_tf_motifs$ChIP_tf_binding[peaks_to_gene_w_tf_motifs$peak %in% peaks_ChIP_overlap$peak_id] = tf
  peaks_to_gene_w_tf_motifs$grn_tf_target[peaks_to_gene_w_tf_motifs$gene %in% genes_tf_targets$gene] = tf
  
  
  # identify whether target genes (and non-targets as background) 
  #     have published ChIP peak in predicted regulatory element, 
  #     collect gene sets for each tf in overlap_list 
  
  t2 = peaks_to_gene_w_tf_motifs
  peaks_tf_targets_peaks_in_published = t2[!is.na(t2$grn_tf_target)&!is.na(t2$ChIP_tf_binding)&!is.na(t2$motif),]
  
  tf_targets_peaks_in_published = unique(peaks_tf_targets_peaks_in_published$gene)
  tf_targets_no_peaks_in_published = unique(t2$gene[!is.na(t2$grn_tf_target)&
                                                      !(t2$gene %in% tf_targets_peaks_in_published)])
  non_targets_peaks_in_published = unique(t2$gene[is.na(t2$grn_tf_target)&!is.na(t2$ChIP_tf_binding)&
                                                    !is.na(t2$motif)])
  non_targets_no_peaks_in_published = unique(t2$gene[is.na(t2$grn_tf_target)&
                                                       !(t2$gene %in% tf_targets_peaks_in_published)])
  
  overlap_list[[tf]] = list(tf_targets_peaks_in_published = tf_targets_peaks_in_published,
                            tf_targets_no_peaks_in_published = tf_targets_no_peaks_in_published,
                            non_targets_peaks_in_published  = non_targets_peaks_in_published,
                            non_targets_no_peaks_in_published = non_targets_no_peaks_in_published)
  
  #collect validated targets with validated TF binding regulatory elements
  
  if (nrow(peaks_tf_targets_peaks_in_published)>0){
    peaks_tf_in_published_list[[tf]] = cbind(tf = tf, peaks_tf_targets_peaks_in_published)
  }
  
}


### save overlapping (validated) targets in table

t1 = bind_rows(peaks_tf_in_published_list)

write_csv(as_tibble(t1), file = paste0(out_dir, script_ind, "Overlap_GRN_predictions_published_datasets_genes.csv"))




###############################################
# calculate stats for overlap targets from GRN with published targets
################################################

for (tf in names(overlap_list)){
  
  l1 = overlap_list[[tf]]
  
  #build table with gene counts
  
  t1 = data.frame(t(lengths(overlap_list[[tf]])))
  
  #check statistical siginificance of enrichment (Fisher-test)
  
  t2 = data.frame(in_published = c(t1$tf_targets_peaks_in_published[1], 
                                   t1$non_targets_peaks_in_published[1]),
                  not_in_published = c(t1$tf_targets_no_peaks_in_published[1],
                                       t1$non_targets_no_peaks_in_published[1])
  )
  rownames(t2) = c("targets_tf", "non_targets")
  
  t3 = fisher.test(t2)
  
  t4 = as_tibble(cbind(tf = tf, t1, 
                       p = t3$p.value, odds_ratio = t3$estimate, 
                       conf95_lower = t3$conf.int[1], conf95_upper = t3$conf.int[2]))
  
  N_overlaps_tab = rbind(N_overlaps_tab, t4)
  
}

N_overlaps_tab$fract_tf_targets_peaks_in_published = N_overlaps_tab$tf_targets_peaks_in_published/
  (N_overlaps_tab$tf_targets_peaks_in_published + N_overlaps_tab$tf_targets_no_peaks_in_published)

N_overlaps_tab$fract_non_targets_peaks_in_published = N_overlaps_tab$non_targets_peaks_in_published/
  (N_overlaps_tab$non_targets_peaks_in_published + N_overlaps_tab$non_targets_no_peaks_in_published)

N_overlaps_tab$padj = p.adjust(N_overlaps_tab$p, method = "BH")

write_csv(N_overlaps_tab, file = paste0(out_dir, script_ind, "Overlap_GRN_predictions_published_datasets_stats.csv"))



###############################################
# plot overlap stats
################################################

# N_overlaps_tab = read_csv(file = paste0(out_dir, script_ind, "Overlap_GRN_predictions_published_datasets_stats.csv"))

t1 = N_overlaps_tab

p1 = ggplot(t1)+geom_point(aes(x = tf, y = fract_tf_targets_peaks_in_published, 
                               size = odds_ratio, color = padj<=0.05))+
  scale_color_manual(limits = c(FALSE, TRUE), values = c("grey", "red"))+
  scale_x_discrete(limits = net_TFs)+
  scale_size_continuous(range = c(0, 3)) +
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pdf(file = paste0(out_dir,script_ind, "Overlap_GRN_predictions_published_datasets_stats.pdf"), 
    width = 5, height = 3)
{
  p1
}

dev.off()



###############################################
# extract high confidence ChIP-validated interactions (consistent cor and expr Z-scores)
################################################

t1 = bind_rows(peaks_tf_in_published_list)

t2 = bulk_data$GRN_DS_reg$edges
net_edges = t2[paste0(t2$tf, t2$target) %in% paste0(t1$tf, t1$gene),]


#update net_nodes (filter only nodes with high confidence interactions, update N_targets)

t2 = bulk_data$GRN_DS_reg$nodes
t2 = t2[t2$gene %in% net_edges$tf | t2$gene %in% net_edges$target,]

t3 = net_edges %>% group_by(tf) %>% dplyr::summarise(N_targets = n())

t2$N_targets = t3$N_targets[match(t2$gene, t3$tf)]
t2$N_targets[is.na(t2$N_targets)] = 0
t2 = t2[order(-t2$N_targets),]

t2$labels = ""
t2$labels[t2$N_targets>0] = paste0(t2$gene," (", t2$N_targets,")")[t2$N_targets>0]

net_nodes = t2

bulk_data$GRN_ChIP_val$nodes = net_nodes
bulk_data$GRN_ChIP_val$edges = net_edges

write_csv(net_nodes, file = paste0(out_dir, script_ind, "network_nodes_ChIP_val.csv"))
write_csv(net_edges, file = paste0(out_dir, script_ind, "network_edges_ChIP_val.csv"))

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_with_GRN_ChIP_val.rda"))



#########################################
#count and classify DEGs and Network
#########################################

v1 = unique(unlist(bulk_data$DEGs))
stats_degs = c(N_genes = length(v1), N_TFs = length(intersect(v1, GOI$TF)), 
               N_Chr21_genes = length(intersect(v1, GOI$Chr21)),
               N_Chr21_TFs = length(intersect(v1, intersect(GOI$Chr21, GOI$TF))),
               N_network_interactions = NA,
               N_network_TFs = NA,
               N_Chr21_TF_targets = NA)

n1 = bulk_data$GRN_unfiltered$nodes
e1 = bulk_data$GRN_unfiltered$edges
stats_unfiltered = c(N_genes = length(n1$gene), N_TFs = length(intersect(n1$gene, GOI$TF)), 
                     N_Chr21_genes = length(intersect(n1$gene, GOI$Chr21)),
                     N_Chr21_TFs = length(intersect(n1$gene, intersect(GOI$Chr21, GOI$TF))),
                     N_network_interactions = nrow(e1),
                     N_network_TFs = length(n1$gene[n1$N_targets>0]),
                     N_Chr21_TF_targets = length(unique(e1$target[e1$tf %in% GOI$Chr21]))
)

n1 = bulk_data$GRN_DS_reg$nodes
e1 = bulk_data$GRN_DS_reg$edges
stats_DS_reg = c(N_genes = length(n1$gene), N_TFs = length(intersect(n1$gene, GOI$TF)), 
                 N_Chr21_genes = length(intersect(n1$gene, GOI$Chr21)),
                 N_Chr21_TFs = length(intersect(n1$gene, intersect(GOI$Chr21, GOI$TF))),
                 N_network_interactions = nrow(e1),
                 N_network_TFs = length(n1$gene[n1$N_targets>0]),
                 N_Chr21_TF_targets = length(unique(e1$target[e1$tf %in% GOI$Chr21]))
)


n1 = bulk_data$GRN_ChIP_val$nodes
e1 = bulk_data$GRN_ChIP_val$edges
stats_ChIP_val = c(N_genes = length(n1$gene), N_TFs = length(intersect(n1$gene, GOI$TF)), 
                   N_Chr21_genes = length(intersect(n1$gene, GOI$Chr21)),
                   N_Chr21_TFs = length(intersect(n1$gene, intersect(GOI$Chr21, GOI$TF))),
                   N_network_interactions = nrow(e1),
                   N_network_TFs = length(n1$gene[n1$N_targets>0]),
                   N_Chr21_TF_targets = length(unique(e1$target[e1$tf %in% GOI$Chr21]))
)

grn_stats = cbind(DEGs_all = stats_degs, GRN_unfiltered = stats_unfiltered,
                  GRN_DS_reg = stats_DS_reg, GRN_ChIP_val = stats_ChIP_val)

grn_stats = as_tibble(cbind(stats = rownames(grn_stats), grn_stats))

write_csv(grn_stats , file = paste0(out_dir, script_ind, "network_stats.csv"))




####################################
# plot full network 
####################################
# node-shape: gene type, node-size: mean vst-norm expression CON clusters, 
# node-colour: mean DELTA Z-score all clusters

net_nodes = bulk_data$GRN_ChIP_val$nodes
net_edges = bulk_data$GRN_ChIP_val$edges

### plot complete network for all clusters

pdf(file = paste0(out_dir,script_ind, "Network_plot_complete_ChIP_val.pdf"), 
    width = 6, height = 7)
{
  
  plot(grn_plot(nodes = net_nodes, edges = net_edges, 
                plot_genes = NULL, 
                node_size_col = "mean_vst_CON", 
                node_color_col = "mean_z_DELTA",
                edge_alpha = 0.5, edgewidth_range = c(0.01, 0.1),
                node_size_limits = c(min(net_nodes$mean_vst_CON), max(net_nodes$mean_vst_CON)*2),
                title = paste0("Full Network (node size ~ mean vst all CON clusters)")))
  
}

dev.off()



###########################################
# plot network (only TFs)
###########################################

p1 = grn_plot(nodes = net_nodes, edges = net_edges, 
              plot_genes = net_nodes$gene[net_nodes$N_targets>0], 
              node_size_col = "mean_vst_CON", 
              node_color_col = "mean_z_DELTA",
              edge_alpha = 0.5, edgewidth_range = c(0.05, 0.5),
              title = paste0("Mean all clusters - TF Network"))


pdf(file = paste0(out_dir,script_ind, "Network_plot_TFs_only_ChIP_val.pdf"), 
    width = 5, height = 6.5)
{
  plot(p1)
}

dev.off()



########################################################################
# plot network (only edges Chr21 TFs => direct targets) 
########################################################################

edges_sub = net_edges[net_edges$tf %in% GOI$Chr21,]

nodes_sub = net_nodes[net_nodes$gene %in% c(edges_sub$tf, edges_sub$target),]
nodes_sub$labels = nodes_sub$gene


p1 = grn_plot(nodes = nodes_sub, edges = edges_sub, layout = "fr",
              plot_genes = NULL, 
              node_size_col = "mean_vst_CON", 
              node_color_col = "mean_z_DELTA",
              node_size_limits = NULL,
              node_color_limits = NULL,
              edge_alpha = 0.5, edgewidth_range = c(0.05, 0.5),
              title = paste0("Mean Z all clusters - Chr21-TFs and direct targets"))



pdf(file = paste0(out_dir,script_ind, "Network_plot_Chr21_TFs_dir_targets_ChIP_val.pdf"), 
    width = 10, height = 10)
{
  plot(p1)
}

dev.off()




########################################################################
# plot network (only edges Chr21 TFs => direct target TFs)
########################################################################

v1 = net_nodes$gene[net_nodes$gene_class == "TF_Chr21"]
v2 = net_nodes$gene[net_nodes$N_targets >0]
edges_sub = net_edges[net_edges$tf %in% v1 & net_edges$target %in% v2,]

nodes_sub = net_nodes[net_nodes$gene %in% c(edges_sub$tf, edges_sub$target),]

p1 = grn_plot(nodes = nodes_sub, edges = edges_sub, layout = "fr",
              plot_genes = NULL, 
              node_size_col = "mean_vst_CON", 
              node_color_col = "mean_z_DELTA",
              node_size_limits = NULL,
              node_color_limits = NULL,
              title = paste0("Mean Z all clusters - Chr21-TFs and direct target TFs"))



pdf(file = paste0(out_dir,script_ind, "Network_plot_Chr21_TFs_dir_target_TFs_interactions_ChIP_val.pdf"), 
    width = 4, height = 6.5)
{
  plot(p1)
}

dev.off()




########################################################################
# plot TF targets vs top 20 diff GO term genes
########################################################################

t1 = bulk_data$GO_results$full
if (nrow(t1)>20){t1 = t1[1:20,]}

go_genes_list = str_split(t1$geneID, "/")
names(go_genes_list) = paste0(t1$Description, " (", t1$ID, ")")

TFs = unique(net_nodes$gene[net_nodes$N_targets>0])

tf_go_mat = matrix(nrow = length(go_genes_list), ncol = length(TFs),
                   dimnames = list(names(go_genes_list), TFs))


for (tf in TFs){
  
  targets = unique(net_edges$target[net_edges$tf == tf])
  
  N_targets = as_vector(lapply(names(go_genes_list), function(go){
    length(intersect(go_genes_list[[go]], targets))
  }))
  
  tf_go_mat[,tf] = N_targets
  
}

colnames(tf_go_mat) = net_nodes$labels[match(colnames(tf_go_mat), net_nodes$gene)]


### plot number of targets in GO Term DEGs for each TF

pdf(file = paste0(out_dir,script_ind, "TF_targets_by_top_20_GO_term_ChIP_val_heatmap.pdf"), width = 20, height = 10)
{
  
  pheatmap::pheatmap(tf_go_mat, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = TRUE, cluster_cols = TRUE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(tf_go_mat), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "N TF targets vs GO_terms"
  )
  
}
dev.off()



#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed F05 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


