message("\n\n##########################################################################\n",
        "# Start H03: scMEGA network visualisation ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
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
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(igraph)
library(ggraph)
library(tidygraph)
library(ArchR)
library(scMEGA)
library(GenomeInfoDb)
library(JASPAR2024)
library(TFBSTools)
library(ComplexHeatmap)


#specify script/output index as prefix for file names
script_ind = "H03_"

#specify output directory
out_dir = paste0(main_dir,"H_scMEGA_GRN_analysis_exc_lin_PCW11_12/")


#load Seurat dataset with GRN

load(file = paste0(out_dir,"H02_seur_with_GRN.rda")) 


#load pseudobulk dataset

load(file = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_PCW11_12/",
                   "E03_bulk_data_w_expr_z_scores.rda"))


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




####################################
# extract network nodes and edges
####################################

###save peak-gene-links

write_csv(df.p2g, file = paste0(out_dir, script_ind, "GRN_peak_gene_links.csv"))


### manually extract network interactions

#keep interactions with padj <0.05 and |TF gene corr| >= 0.3, only for Chr21 and DEG TFs 
t1 = df.grn[!is.na(df.grn$fdr) & df.grn$fdr <= 0.05 & abs(df.grn$correlation) >= 0.3 ,]
t1 = dplyr::rename(t1, target = gene)
t1 = t1[t1$tf %in% unlist(bulk_data$DEGs),]

# generate table with network edges
t1$from = t1$tf
t1$to = t1$target

#keep only edges with expressed genes
m1 = bulk_data$gene_Z_scores$by_cluster_sample
t1 = t1[t1$from %in% rownames(m1) & t1$to %in% rownames(m1),]

net_edges = t1


#extract node/gene metadata (TF, Chr21 genes, number of targets, label for plots "[gene]([N_targets])" for TFs, no label for target genes )

t2 = tibble(gene = unique(c(net_edges$tf, net_edges$target)), TF = "", Chr21 = "", 
            gene_class = "other", N_targets = 0, labels = "", DEG = "")

t2$TF[t2$gene %in% c(net_edges$tf, GOI$TF)] = "TF"
t2$Chr21[t2$gene %in% GOI$Chr21] = "Chr21"
t2$gene_class[t2$TF == "TF"] = "TF"
t2$gene_class[t2$Chr21 == "Chr21"] = "Chr21"
t2$gene_class[t2$TF == "TF" & t2$Chr21 == "Chr21"] = "TF_Chr21"
t2$DEG[t2$gene %in% unlist(bulk_data$DEGs)] = "DEG"

t3 = net_edges[!(duplicated(paste0(net_edges$tf,net_edges$target))),] %>% group_by(tf) %>% dplyr::summarise(N_targets = n())

t2$N_targets = t3$N_targets[match(t2$gene, t3$tf)]
t2$N_targets[is.na(t2$N_targets)] = 0
t2 = t2[order(-t2$N_targets),]
t2$labels[t2$TF !=""] = paste0(t2$gene," (", t2$N_targets,")")[t2$TF !=""]

net_nodes = t2

net_edges = net_edges[order(match(net_edges$tf, net_nodes$gene)),]

bulk_data$GRN_unfiltered$nodes = net_nodes
bulk_data$GRN_unfiltered$edges = net_edges


###count and classify DEGs and Network

v1 = unique(unlist(bulk_data$DEGs))
t1 = c(N_genes = length(v1), N_TFs = length(intersect(v1, GOI$TF)), 
       N_Chr21_genes = length(intersect(v1, GOI$Chr21)),
       N_Chr21_TFs = length(intersect(v1, intersect(GOI$Chr21, GOI$TF))),
       N_network_interactions = NA,
       N_network_TFs = NA,
       N_Chr21_TF_targets = NA)

v2 = net_nodes$gene
t2 = c(N_genes = length(v2), N_TFs = length(intersect(v2, GOI$TF)), 
       N_Chr21_genes = length(intersect(v2, GOI$Chr21)),
       N_Chr21_TFs = length(intersect(v2, intersect(GOI$Chr21, GOI$TF))),
       N_network_interactions = nrow(net_edges),
       N_network_TFs = length(net_nodes$gene[net_nodes$N_targets>0]),
       N_Chr21_TF_targets = length(unique(net_edges$target[net_edges$tf %in% GOI$Chr21]))
       )

grn_stats = cbind(all_DEGs = t1, GRN_unfiltered = t2)


####################################
# plot full network (unfiltered)
####################################
# node-shape: gene type, node-size: mean vst-norm expression CON clusters, 
# node-colour: mean DELTA Z-score all clusters


### calculate mean Z-score for all CON clusters and DELTA Z-score DS-CON

meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)

m1 = bulk_data$gene_Z_scores$by_cluster_sample
m_CON = m1[net_nodes$gene, meta$cluster_sample[meta$group == "CON"]]
m_DS = m1[net_nodes$gene, meta$cluster_sample[meta$group == "DS"]]
mean_z_CON = apply(m_CON, 1, mean)
mean_z_DS = apply(m_DS, 1, mean)
mean_z_DELTA = mean_z_DS-mean_z_CON

net_nodes$mean_z_CON = mean_z_CON
net_nodes$mean_z_DS = mean_z_DS
net_nodes$mean_z_DELTA = mean_z_DELTA

bulk_data$GRN$nodes = net_nodes


### calculate mean vst for all control clusters combined and individually (to visualise absolute expression by node size)

m3 = bulk_data$vst_mat[, meta$cluster_sample[meta$group == "CON"]]
bulk_data$vst_mean$CON_all = apply(m3, 1, mean)

l1 = lapply(cluster_names, function(cl){
  m3 = bulk_data$vst_mat[, meta$cluster_sample[meta$group == "CON" & meta$cluster_name == cl] ]
  m3 = apply(m3, 1, mean)
  })
names(l1) = cluster_names
t1 = as.matrix(as.data.frame(l1))
bulk_data$vst_mean$CON_by_cluster_name = t1

net_nodes$mean_z_CON = mean_z_CON
net_nodes$mean_z_DS = mean_z_DS
net_nodes$mean_z_DELTA = mean_z_DELTA
net_nodes$mean_vst_CON = bulk_data$vst_mean$CON_all[net_nodes$gene]

bulk_data$GRN$nodes = net_nodes

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_with GRN_mean_vst.rda"))


### plot complete network

pdf(file = paste0(out_dir,script_ind, "Network_plot_complete_unfiltered.pdf"), 
    width = 10, height = 10)
{
  
  plot(grn_plot(nodes = net_nodes, edges = net_edges, 
                plot_genes = NULL, 
                node_size_col = "mean_vst_CON", 
                node_color_col = "mean_z_DELTA",
                edge_alpha = 0.5, edgewidth_range = c(0.01, 0.1),
                node_size_limits = c(min(net_nodes$mean_vst_CON), max(net_nodes$mean_vst_CON)*3),
                title = paste0("Full Network (node size ~ mean vst all CON clusters)")))
  
}

dev.off()



###########################################
# plot network (only TFs) (unfiltered)
###########################################

p1 = grn_plot(nodes = net_nodes, edges = net_edges, 
              plot_genes = net_nodes$gene[net_nodes$N_targets>0], 
              node_size_col = "mean_vst_CON", 
              node_color_col = "mean_z_DELTA",
              node_size_limits = NULL,
              node_color_limits = NULL,
              edge_alpha = 0.5, edgewidth_range = c(0.05, 0.5),
              title = paste0("Mean all clusters - TF Network"))


pdf(file = paste0(out_dir,script_ind, "Network_plot_TFs_unfiltered.pdf"), 
    width = 5.5, height = 6)
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
              edge_alpha = 0.5,
              title = paste0("Mean Z all clusters - Chr21-TFs and direct targets"))



pdf(file = paste0(out_dir,script_ind, "Network_plot_Chr21_TFs_dir_targets_unfiltered.pdf"), 
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
                                 edge_alpha = 0.7,
                                 title = paste0("Mean Z all clusters - Chr21-TFs and direct target TFs"))



pdf(file = paste0(out_dir,script_ind, "Network_plot_Chr21_TFs_dir_target_TFs_interactions_unfiltered.pdf"), 
    width = 4, height = 6)
{
  plot(p1)
}

dev.off()




########################################################################
# plot TF targets vs top 20 diff GO term genes (unfiltered)
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

pdf(file = paste0(out_dir,script_ind, "TF_targets_by_top_20_GO_term_heatmap_unfiltered.pdf"), width = 20, height = 10)
{
  pheatmap::pheatmap(tf_go_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(tf_go_mat), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "N TF targets vs GO_terms"
  )
  pheatmap::pheatmap(tf_go_mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(tf_go_mat), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "N TF targets vs GO_terms"
  )
  
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


### plot networks of GO genes with upstream TFs

pl_list = list()
pl_list_tree = list()

for (go in names(go_genes_list)){
  
  go_genes = go_genes_list[[go]]
  
  upstream_TFs = unique(net_edges$tf[net_edges$target %in% go_genes])
  
  #include only interactions TFs => GO genes
  edges_sub = net_edges[net_edges$tf %in% upstream_TFs &
                          (net_edges$target %in% go_genes),]
  
  nodes_sub = net_nodes[net_nodes$gene %in% c(edges_sub$tf, edges_sub$target),]
  nodes_sub$labels[nodes_sub$labels == ""] = nodes_sub$gene[nodes_sub$labels == ""]
  
  l1 = bulk_data$gene_Z_scores$mean_by_cluster_name
  
  net_size_limits = range(bulk_data$vst_mean$CON_by_cluster_name)
  net_color_limits = c(-max(abs(l1$DELTA)), max(abs(l1$DELTA)) )
  
  pl_list[[go]] = grn_plot(nodes = nodes_sub, edges = edges_sub, layout = "fr",
                                   plot_genes = NULL, 
                                   node_size_col = "mean_vst_CON", 
                                   node_color_col = "mean_z_DELTA",
                                   node_size_limits = net_size_limits,
                                   node_color_limits = net_color_limits, edge_alpha = 1,
                                   title = paste0(go, " - Mean Z all clusters - with upstream TFs"))
  
}



pdf(file = paste0(out_dir,script_ind, "Network_plot_top20_GO_genes_with_upstream_TFs_unfiltered.pdf"), 
    width = 12, height = 8)
{
  lapply(pl_list, function(x){x})
}

dev.off()





########################################################################
# filter network to retain only interactions consistent with DS-dependent regulation
########################################################################

net_edges$z_DELTA_tf = net_nodes$mean_z_DELTA[match(net_edges$tf, net_nodes$gene)]
net_edges$z_DELTA_target = net_nodes$mean_z_DELTA[match(net_edges$target, net_nodes$gene)]
net_edges$DS_reg = net_edges$correlation*net_edges$z_DELTA_tf*net_edges$z_DELTA_target>0

write_csv(net_edges, file = paste0(out_dir, script_ind, "network_edges.csv"))

net_edges = net_edges[net_edges$DS_reg,]


#update net_nodes (filter only nodes with high confidence interactions, update N_targets)

t2 = net_nodes
t2$N_targets_unfiltered = t2$N_targets

t3 = net_edges %>% group_by(tf) %>% dplyr::summarise(N_targets = n())

t2$N_targets = t3$N_targets[match(t2$gene, t3$tf)]
t2$N_targets[is.na(t2$N_targets)] = 0
t2 = t2[order(-t2$N_targets),]

t2$labels = ""
t2$labels[t2$N_targets>0] = paste0(t2$gene," (", t2$N_targets,")")[t2$N_targets>0]

t2$DS_reg = t2$gene %in% c(net_edges$tf, net_edges$target)

net_nodes = t2

write_csv(net_nodes, file = paste0(out_dir, script_ind, "network_nodes.csv"))

net_nodes = t2[t2$gene %in% c(net_edges$tf, net_edges$target),]

bulk_data$GRN_DS_reg$nodes = net_nodes
bulk_data$GRN_DS_reg$edges = net_edges

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_with GRN_mean_vst.rda"))


###count and classify DEGs and Network

v2 = net_nodes$gene
t2 = c(N_genes = length(v2), N_TFs = length(intersect(v2, GOI$TF)), 
       N_Chr21_genes = length(intersect(v2, GOI$Chr21)),
       N_Chr21_TFs = length(intersect(v2, intersect(GOI$Chr21, GOI$TF))),
       N_network_interactions = nrow(net_edges),
       N_network_TFs = length(net_nodes$gene[net_nodes$N_targets>0]),
       N_Chr21_TF_targets = length(unique(net_edges$target[net_edges$tf %in% GOI$Chr21]))
)

grn_stats2 = as_tibble(cbind(stats = rownames(grn_stats), grn_stats, GRN_DS_reg = t2))

write_csv(grn_stats2 , file = paste0(out_dir, script_ind, "network_stats.csv"))


#############################################
# plot complete network (only interactions consistent with DS-dependent regulation)
#############################################

pdf(file = paste0(out_dir,script_ind, "Network_plot_complete_DS_reg.pdf"), 
    width = 10, height = 10)
{
  
  plot(grn_plot(nodes = net_nodes, edges = net_edges, 
                plot_genes = NULL, 
                node_size_col = "mean_vst_CON", 
                node_color_col = "mean_z_DELTA",
                edge_alpha = 0.5, edgewidth_range = c(0.01, 0.1),
                node_size_limits = c(min(net_nodes$mean_vst_CON), max(net_nodes$mean_vst_CON)*5),
                title = paste0("Full Network (node size ~ mean vst all CON clusters)")))
  
}

dev.off()



###########################################
# plot network (only TFs) (only interactions consistent with DS-dependent regulation)
###########################################

p1 = grn_plot(nodes = net_nodes, edges = net_edges, 
              plot_genes = net_nodes$gene[net_nodes$N_targets>0], 
              node_size_col = "mean_vst_CON", 
              node_color_col = "mean_z_DELTA",
              node_size_limits = NULL,
              node_color_limits = NULL,
              edge_alpha = 0.5, edgewidth_range = c(0.05, 0.5),
              title = paste0("Mean all clusters - TF Network"))


pdf(file = paste0(out_dir,script_ind, "Network_plot_TFs_DS_reg.pdf"), 
    width = 5, height = 6)
{
  plot(p1)
}

dev.off()



########################################################################
# plot network (only edges Chr21 TFs => direct targets) (only interactions consistent with DS-dependent regulation)
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
              edge_alpha = 0.5,
              title = paste0("Mean Z all clusters - Chr21-TFs and direct targets"))



pdf(file = paste0(out_dir,script_ind, "Network_plot_Chr21_TFs_dir_targets_DS_reg.pdf"), 
    width = 10, height = 10)
{
  plot(p1)
}

dev.off()


########################################################################
# plot network (only edges Chr21 TFs => direct target TFs) (only interactions consistent with DS-dependent regulation)
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
              edge_alpha = 0.7,
              title = paste0("Mean Z all clusters - Chr21-TFs and direct target TFs"))



pdf(file = paste0(out_dir,script_ind, "Network_plot_Chr21_TFs_dir_target_TFs_interactions_DS_reg.pdf"), 
    width = 4, height = 6)
{
  plot(p1)
}

dev.off()



########################################################################
# plot TF targets vs top 20 diff GO term genes (only interactions consistent with DS-dependent regulation)
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

pdf(file = paste0(out_dir,script_ind, "TF_targets_by_top_20_GO_term_heatmap_DS_reg.pdf"), width = 20, height = 10)
{
  pheatmap::pheatmap(tf_go_mat, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(tf_go_mat), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "N TF targets vs GO_terms"
  )
  pheatmap::pheatmap(tf_go_mat, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = TRUE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(tf_go_mat), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "N TF targets vs GO_terms"
  )
  
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


### plot networks of GO genes with upstream TFs

pl_list = list()
pl_list_tree = list()

for (go in names(go_genes_list)){
  
  go_genes = go_genes_list[[go]]
  
  upstream_TFs = unique(net_edges$tf[net_edges$target %in% go_genes])
  
  #include only interactions TFs => GO genes
  edges_sub = net_edges[net_edges$tf %in% upstream_TFs &
                          (net_edges$target %in% go_genes),]
  
  nodes_sub = net_nodes[net_nodes$gene %in% c(edges_sub$tf, edges_sub$target),]
  nodes_sub$labels[nodes_sub$labels == ""] = nodes_sub$gene[nodes_sub$labels == ""]
  
  l1 = bulk_data$gene_Z_scores$mean_by_cluster_name
  
  net_size_limits = range(bulk_data$vst_mean$CON_by_cluster_name)
  net_color_limits = c(-max(abs(l1$DELTA)), max(abs(l1$DELTA)) )
  
  pl_list[[go]] = grn_plot(nodes = nodes_sub, edges = edges_sub, layout = "fr",
                           plot_genes = NULL, 
                           node_size_col = "mean_vst_CON", 
                           node_color_col = "mean_z_DELTA",
                           node_size_limits = net_size_limits,
                           node_color_limits = net_color_limits, edge_alpha = 1,
                           title = paste0(go, " - Mean Z all clusters - with upstream TFs"))
  
}



pdf(file = paste0(out_dir,script_ind, "Network_plot_top20_GO_genes_with_upstream_TFs_DS_reg.pdf"), 
    width = 12, height = 8)
{
  lapply(pl_list, function(x){x})
}

dev.off()





#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed H03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")



