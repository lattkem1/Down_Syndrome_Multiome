message("\n\n##########################################################################\n",
        "# Start H06: PPI interactions with scMEGA network core TFs ", Sys.time(),
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
library(ggraph)
library(tidygraph)
library(jsonlite)
library(httr)
library(ArchR)
library(scMEGA)

#specify script/output index as prefix for file names
script_ind = "H06_"

#specify output directory
out_dir = paste0(main_dir,"H_scMEGA_GRN_analysis_exc_lin_PCW11_12/")


#load Seurat dataset with GRN

load(file = paste0(out_dir,"H02_seur_with_GRN.rda")) 


#load pseudobulk dataset

load(file = paste0(out_dir,"H03_bulk_data_with GRN_mean_vst.rda")) 



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
    p2 = c("grey", "blue")
  } else if (v2 ==3){
    p2 = c("blue", "grey", "orange")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}



#function: extract network interactions between 2 groups of genes via max one intermediate interactor 

extract_from_x_to_interactions = function(edges, from, to){
  
  from_x_edges = edges[edges$from %in% from, ]
  from_x_edges$x = from_x_edges$to
  x_to_edges = edges[edges$to %in% to, ]
  x_to_edges$x = x_to_edges$from
  from_x_to = merge(from_x_edges[,c("from", "x")], x_to_edges[,c("x", "to")], by = "x")
  from_x_to = from_x_to[,c("from", "x", "to")]
  from_x_to_edges = edges[(edges$from %in% from & edges$to %in% to)|
                            (edges$interaction %in% c(paste0(from_x_to$from, "_",from_x_to$x),
                                                      paste0(from_x_to$x, "_",from_x_to$to)
                            ) ),]
  
  return(from_x_to_edges)
  
}


### function: plot network graph 
# default: node-shape: gene type, node-size: mean Z-score all CON clusters, node-colour: mean DELTA Z-score all clusters
# nodes/edges: full network
# plot_nodes: nodes to plot 
# node_size: name of column specifying node size
# node_color: name of column specifying node color


ppi_tf_act_plot = function(nodes, edges, plot_genes = NULL, layout='fr', 
                           node_size_col = NULL, node_color_col = NULL, node_fill_col = NULL,
                           node_size_limits = NULL, node_color_limits = NULL, node_fill_limits = NULL,
                           edge_alpha = 1, edgewidth_range = c(0.1, 1),
                           plot_unconnected = FALSE, title = "PPI plot"){
  
  if (is.null(plot_genes)){plot_genes = nodes$gene} 
  if (is.null(node_size_col)){node_size_col = "mean_z_CON"}
  if (is.null(node_color_col)){node_color_col = "mean_z_DELTA"}
  if (is.null(node_fill_col)){node_fill_col = "mean_DELTA_act_tf"}
  if (is.null(node_size_limits)){node_size_limits = c(0, max(nodes[[node_size_col]]))}
  if (is.null(node_color_limits)){node_color_limits = c(-max(abs(na.omit(nodes[[node_color_col]]))), 
                                                        +max(abs(na.omit(nodes[[node_color_col]]))))}
  if (is.null(node_fill_limits)){node_fill_limits = c(-max(abs(na.omit(nodes[[node_fill_col]]))), 
                                                      +max(abs(na.omit(nodes[[node_fill_col]]))))}
  plot_nodes = nodes[nodes$gene %in% plot_genes,]
  plot_nodes$node_size = plot_nodes[[node_size_col]]
  plot_nodes$node_color = plot_nodes[[node_color_col]]
  plot_nodes$node_fill = plot_nodes[[node_fill_col]]
  
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
  
  p1 = ggraph(g1, layout=layout) +  
    labs(title = title) + 
    geom_edge_hive(aes(linetype=cor<0, edge_width = abs(cor)), strength = 0,
                   edge_color = "grey30",
                   arrow = arrow(angle = 10, length = unit(3, "mm"),
                                 ends = "last", type = "closed"),
                   start_cap = circle(1, 'mm'),
                   end_cap = circle(1, 'mm'), alpha = edge_alpha) +
    scale_edge_width_continuous(limits = c(0,1), range = edgewidth_range)+
    scale_edge_linetype_manual(limits = c(FALSE, TRUE), values = c("solid", "longdash"))+
    geom_node_point(aes(fill = node_fill, color = node_color, shape = gene_class, 
                        size = node_size, stroke = 0.5 + node_size/10)) +
    scale_size_continuous(limits = node_size_limits, range = c(0.05 , 5), name = node_size_col)+
    geom_node_text(aes(label=label), size=3, repel=T, fontface = "bold") +
    scale_fill_gradient2(limits = node_fill_limits, low = "blue", mid = "white", high = "red", 
                         midpoint = 0, na.value = "white", name = node_fill_col)+
    scale_color_gradient2(limits = node_color_limits, low = "blue", mid = "grey80", high = "red", 
                          midpoint = 0, name = node_color_col)+
    scale_shape_manual(limits = gene_classes, values = gene_class_shapes)+
    theme_void()
  
  return(p1)
  
}




####################################
# extract network nodes and edges
####################################

net_edges = bulk_data$GRN_DS_reg$edges

net_nodes = bulk_data$GRN_DS_reg$nodes



###############################################
# get TF and HSA21 interactors from BioGrid (via REST API) 
################################################

Chr21_DEGs = intersect(unlist(bulk_data$DEGs), GOI$Chr21)
net_TFs = net_nodes$gene[net_nodes$N_targets > 0]


query_genes = unique(c(net_TFs, Chr21_DEGs))
N_query_genes = length(query_genes)

interact_tab = NULL

for (i in 1:N_query_genes){
  
  if(i%%10 == 1){
    message("\n    *** Querying gene ",i," of ", N_query_genes)
  }
  
  url2 = paste0("https://webservice.thebiogrid.org/interactions/",
                "?searchNames=true&geneList=",query_genes[i],
                "&taxId=9606&includeInteractors=true&includeInteractorInteractions=false",
                "&selfInteractionsExcluded=true",
                "&accesskey=06c072eaaddbf340747a5eabea18481d&format=json")
  
  resp <- GET(url2)
  l2 <- content(resp, as = "text") %>% fromJSON()
  
  l3 = lapply(l2, function(x){
    x = x[names(x)!= "ONTOLOGY_TERMS"] # ONTOLOGY_TERMS is list => generates complications in tables
    x = unlist(x)
  })
  
  t4 = as_tibble(t(as.data.frame(l3)))
  
  if (nrow(t4) > 0){
    interact_tab = rbind(interact_tab, t4)
  }
  
}

interact_tab$date_database_access = Sys.time()

interact_tab$from = interact_tab$OFFICIAL_SYMBOL_A
interact_tab$to = interact_tab$OFFICIAL_SYMBOL_B
interact_tab$interaction = paste0(interact_tab$from,"_",interact_tab$to)

write_csv(interact_tab, file = paste0(out_dir, script_ind,"BioGrid_core_all_interactions_net_TFs_Chr21_DEGs.csv"))


###identify direct interactions of differential Chr21 genes and network TFs

t1 = interact_tab

t2 = t1[(t1$to %in% net_TFs & t1$from %in% Chr21_DEGs),]

Chr21_TF_interactions = t2

write_csv(Chr21_TF_interactions, file = paste0(out_dir, script_ind,"BioGrid_direct_Chr21_TF_interactions.csv"))


### identify indirect interactions between TFs and differential HSA21 genes via single intermediate factors

Chr21_X_TF_interactions = extract_from_x_to_interactions(edges = interact_tab, from = Chr21_DEGs, to = net_TFs)

write_csv(Chr21_X_TF_interactions, file = paste0(out_dir, script_ind,"BioGrid_Chr21_X_TF_interactions.csv"))


###############################################
#extract PPI network nodes and edges
###############################################

###prepare node table: take TFs, HSA21 genes and interactors from BioGrid analysis, 
#       add number of TF targets and mean CON vst (for node size) /DELTA Z-scores (for node color)

t2 = tibble(gene = unique(c(Chr21_X_TF_interactions$from, 
                            Chr21_X_TF_interactions$to)), 
            TF = "", Chr21 = "", gene_class = "other", N_targets = 0, label = "")

t2$TF[t2$gene %in% c(net_edges$tf, GOI$TF)] = "TF"
t2$Chr21[t2$gene %in% GOI$Chr21] = "Chr21"
t2$gene_class[t2$TF == "TF"] = "TF"
t2$gene_class[t2$Chr21 == "Chr21"] = "Chr21"
t2$gene_class[t2$TF == "TF" & t2$Chr21 == "Chr21"] = "TF_Chr21"
t2$N_targets = net_nodes$N_targets[match(t2$gene, net_nodes$gene)]
t2$N_targets[is.na(t2$N_targets)] = 0
t2$label = t2$gene
t2$label[t2$N_targets>0] = paste0(t2$gene, " (",t2$N_targets,")")[t2$N_targets>0]

t3 = bulk_data$vst_mean$CON_all
t2$mean_vst_CON = t3[t2$gene]

PPI_nodes = t2

meta = bulk_data$meta

m1 = bulk_data$gene_Z_scores$by_cluster_sample
m_CON = m1[match(PPI_nodes$gene, rownames(m1)), meta$cluster_sample[meta$group == "CON"]]
m_DS = m1[match(PPI_nodes$gene, rownames(m1)), meta$cluster_sample[meta$group == "DS"]]
mean_z_CON = apply(m_CON, 1, mean)
mean_z_DS = apply(m_DS, 1, mean)
mean_z_DELTA = mean_z_DS-mean_z_CON

PPI_nodes$mean_z_CON = mean_z_CON
PPI_nodes$mean_z_DS = mean_z_DS
PPI_nodes$mean_z_DELTA = mean_z_DELTA

PPI_nodes = PPI_nodes[!is.na(PPI_nodes$mean_vst_CON),]


###prepare edges (with number of experiments demonstrating interaction, type of interaction)

t1 = Chr21_X_TF_interactions
t2 = t1 %>% group_by(from, to, interaction) %>% 
  dplyr::summarise (N_evidence = n(),
                    interaction_type = paste0(unique(EXPERIMENTAL_SYSTEM_TYPE), collapse = "_"),
                    modification = paste0(unique(MODIFICATION), collapse = "_"))
#fix interaction_type and modification naming
t2$interaction_type[t2$interaction_type == "genetic_physical"] = "physical_genetic"
t2$modification = str_remove(t2$modification, "-_")
t2$modification = str_remove(t2$modification, "_-")
t2$modification[t2$modification == "-"] = t2$interaction_type[t2$modification == "-"]

t2$edge_class = t2$modification

PPI_edges = t2

save(net_nodes, net_edges, PPI_nodes, PPI_edges, file = paste0(out_dir, script_ind, "GRN_PPI_network.rda"))


###############################################
#for each network TF, identify all potential upstream Chr21 genes
###############################################

load(file = paste0(out_dir, script_ind, "GRN_PPI_network.rda"))


Chr21_net_TF_links = NULL

for (tf in net_TFs){
  
  t1 = PPI_edges[PPI_edges$to == tf,] #TF interactors
  Chr21_tf_links = tibble(linked_Chr21_gene = t1$from[t1$from %in% GOI$Chr21], interactor = "", tf = tf)
  t2 = PPI_edges[PPI_edges$to %in% t1$from,] #interactions of TF interactors with Chr21 genes
  Chr21_X_tf_links = tibble(linked_Chr21_gene = t2$from[t2$from %in% GOI$Chr21], 
              interactor = t2$to[t2$from %in% GOI$Chr21], tf = tf)
  
  Chr21_net_TF_links = rbind(Chr21_net_TF_links, Chr21_tf_links, Chr21_X_tf_links)
  
}



###############################################
#calculate correlation of TF activity with expression of linked Chr21 genes
###############################################

### get chromvar TF activity data

chrom_var_mat <- seur[["chromvar"]]$data

#convert rownames (motif name to TF name)
identical(rownames(chrom_var_mat), names(seur@assays$peaks_by_cluster@motifs@motif.names))
rownames(chrom_var_mat) = seur@assays$peaks_by_cluster@motifs@motif.names

#create trajectory bin matrix and calculate mean TF activity per bin
meta = seur@meta.data
meta = meta[order(meta$Trajectory),]
meta$Trajectory_bin = ceiling(meta$Trajectory)

chrom_var_traj_mat = matrix(nrow = nrow(chrom_var_mat), ncol = 100, 
                            dimnames = list(rownames(chrom_var_mat), paste0("T", 1:100)))

for (i in 1:100){
  m1 = chrom_var_mat[, rownames(meta)[meta$Trajectory_bin == i]]
  m2 = apply(m1, 1, mean)
  chrom_var_traj_mat[,i] = m2
}

#keep only TFs with Chr21 interactions
chrom_var_traj_mat = chrom_var_traj_mat[unique(Chr21_net_TF_links$tf),]


### get expression data for all Chr21 genes linked to TFs

t1 <- seur[["SCT"]]$data
expr_mat = t1[rownames(t1) %in% Chr21_net_TF_links$linked_Chr21_gene,]

expr_traj_mat = matrix(nrow = nrow(expr_mat), ncol = 100, 
                       dimnames = list(rownames(expr_mat), paste0("T", 1:100)))

for (i in 1:100){
  m1 = expr_mat[, rownames(meta)[meta$Trajectory_bin == i]]
  m2 = apply(m1, 1, mean)
  expr_traj_mat[,i] = m2
}

#calculate expr z-score per gene
t1 = apply(expr_traj_mat, 1, scale)
t2 = t(t1)
colnames(t2) = colnames(expr_traj_mat)
expr_traj_mat = t2


### calculate corr Chr21 gene expression with TF activity for all interactions

t1 = Chr21_net_TF_links
t1 = t1[t1$linked_Chr21_gene %in% rownames(expr_traj_mat) & t1$tf %in% rownames(chrom_var_traj_mat),]

for (i in 1:nrow(t1)){
  t2 = cor.test(expr_traj_mat[t1$linked_Chr21_gene[i],], chrom_var_traj_mat[t1$tf[i],])
  t1$t[i] = t2$statistic
  t1$p[i] = t2$p.value
  t1$cor[i] = t2$estimate
  t1$conf_min[i] = t2$conf.int[1]
  t1$conf_max[i] = t2$conf.int[2]
}

t1$padj = p.adjust(t1$p, method = "BH")

Chr21_net_TF_links_cor = t1

write_csv(Chr21_net_TF_links_cor, file = paste0(out_dir, script_ind,"BioGrid_Chr21_X_TF_interactions_with_cor.csv"))


### extract interactions with significant correlation (t1$padj <=0.05 & abs(t1$cor)>0.2, and |z-score|>0.1)
#       and deregulation in DS consistent with direction of correlation

# calculate Z-score DS vs CON cells for expression Chr21 genes and activity TFs

t3 = t(apply(expr_mat, 1, scale))
colnames(t3) = colnames(expr_mat)
mean_DELTA_expr_Chr21 = apply(t3[, rownames(meta)[meta$group == "DS"]], 1, mean)-
  apply(t3[, rownames(meta)[meta$group == "CON"]], 1, mean)

t3 = t(apply(chrom_var_mat, 1, scale))
colnames(t3) = colnames(chrom_var_mat)
mean_DELTA_act_TFs = apply(t3[, rownames(meta)[meta$group == "DS"]], 1, mean)-
  apply(t3[, rownames(meta)[meta$group == "CON"]], 1, mean)

# add TF activity and Chr21 gene expression Z-scores to correlation analysis results
# filter for sign correlation and consistent changes (cor * TF z-score * Chr21 z-score >0 )

t1 = Chr21_net_TF_links_cor
t1$mean_DELTA_expr_Chr21 = mean_DELTA_expr_Chr21[t1$linked_Chr21_gene]
t1$mean_DELTA_act_tf = mean_DELTA_act_TFs[t1$tf]

t2 = t1[(t1$padj <=0.05 & abs(t1$cor)>=0.2 & abs(t1$mean_DELTA_act_tf)>= 0.01 &
          abs(t1$mean_DELTA_expr_Chr21)>= 0.1 &
          t1$cor*t1$mean_DELTA_expr_Chr21*t1$mean_DELTA_act_tf > 0),]

Chr21_net_TF_links_high_conf = t2

write_csv(Chr21_net_TF_links_high_conf, file = paste0(out_dir, script_ind,"BioGrid_Chr21_X_TF_interactions_high_conf.csv"))



###############################################
#keep only PPI interactions with significant correlation between Chr21 expr and TF activity
#     collapse direct/indirect interactions to one interaction per Chr21-TF pair 
#     add interactors as labels of Chr21 genes
###############################################

PPI_nodes_all = PPI_nodes
PPI_edges_all = PPI_edges

### filter edges

t1 = Chr21_net_TF_links_high_conf
t1$interaction_Chr21_tf = paste0(t1$linked_Chr21_gene, "_", t1$tf)
t1$interaction_Chr21_int = paste0(t1$linked_Chr21_gene, "_", t1$interactor)
t1$interaction_int_tf = paste0(t1$interactor, "_", t1$tf)
t1$from = paste0(t1$linked_Chr21_gene, "(=> ", t1$interactor, " =>) ")
t1$to = t1$tf

#collapse all links from each Chr21-TF pair and create network edges table

t2 = t1 %>% group_by(linked_Chr21_gene, tf, cor, interaction_Chr21_tf) %>% 
  dplyr::summarise(interactors = paste0(interactor, collapse = "/"))

t2$from = t2$linked_Chr21_gene
t2$to = t2$tf

PPI_edges = t2

###update nodes

PPI_nodes = PPI_nodes[(PPI_nodes$gene %in% PPI_edges$from)|(PPI_nodes$gene %in% PPI_edges$to),]

PPI_nodes$mean_DELTA_act_tf = mean_DELTA_act_TFs[PPI_nodes$gene]



############################################################
#summarise network interactions
############################################################

t1 = tibble(total_interactions = nrow(Chr21_net_TF_links), 
            interactions_DS_reg = nrow(Chr21_net_TF_links_high_conf),
            Chr21_gene_TF_links_DS_reg = nrow(PPI_edges),
            linked_Chr21_genes_DS_reg = length(unique(PPI_edges$linked_Chr21_gene)),
            linked_TFs_DS_reg = length(unique(PPI_edges$tf)))


write_csv(t1, file = paste0(out_dir, script_ind,"BioGrid_Chr21_X_TF_interactions_stats.csv"))


###############################################
# plot HSA21-X-TF PPI networks (with correlation and TF activity, without interactors)
################################################


pdf(file = paste0(out_dir,script_ind, "PPI_Network_plot_Chr21_TF_links.pdf"), 
    width = 6, height = 8.5)
{
  
  plot(ppi_tf_act_plot(nodes = PPI_nodes, edges = PPI_edges, 
                plot_genes = NULL, 
                node_size_col = "mean_vst_CON", 
                node_color_col = "mean_z_DELTA",
                node_fill_col = "mean_DELTA_act_tf",
                node_size_limits = c(0, max(PPI_nodes$mean_vst_CON)),
                edge_alpha = 0.5, edgewidth_range = c(0.1, 1),
                title = paste0("Full Network")))
  
}

dev.off()




###############################################
# plot HSA21-X-TF PPI networks (without correlation, with TF activity and interactors)
################################################

PPI_nodes = PPI_nodes_all[PPI_nodes_all$gene %in% 
                            c(Chr21_net_TF_links_high_conf$linked_Chr21_gene,
                              Chr21_net_TF_links_high_conf$interactor,
                              Chr21_net_TF_links_high_conf$tf),]
PPI_nodes$mean_DELTA_act_tf = Chr21_net_TF_links_high_conf$mean_DELTA_act_tf[match(PPI_nodes$gene,
                                                                                   Chr21_net_TF_links_high_conf$tf)]

PPI_edges = Chr21_X_TF_interactions[(Chr21_X_TF_interactions$from %in% PPI_nodes$gene) & 
                            (Chr21_X_TF_interactions$to %in% PPI_nodes$gene),]
PPI_edges$cor = 0.5 #dummy correlation required for plotting function

### plot complete network for all clusters

pdf(file = paste0(out_dir,script_ind, "PPI_Network_plot_with_interactors.pdf"), 
    width = 8, height = 9)
{
  
  plot(ppi_tf_act_plot(nodes = PPI_nodes, edges = PPI_edges, 
                       plot_genes = NULL, 
                       node_size_col = "mean_vst_CON", 
                       node_color_col = "mean_z_DELTA",
                       node_fill_col = "mean_DELTA_act_tf",
                       node_size_limits = c(0, max(PPI_nodes$mean_vst_CON)),
                       edge_alpha = 0.5, edgewidth_range = c(0.05, 0.5),
                       title = paste0("Full Network")))
  
}

dev.off()




###############################################
# plot HSA21-X-TF PPI networks by HSA21 gene (without correlation, with TF activity and interactors)
################################################

pl= list()

for (Chr21_gene in unique(Chr21_net_TF_links_high_conf$linked_Chr21_gene)){
  
  t1 = Chr21_net_TF_links_high_conf[Chr21_net_TF_links_high_conf$linked_Chr21_gene == Chr21_gene,]
  
  PPI_nodes = PPI_nodes_all[PPI_nodes_all$gene %in% 
                              c(t1$linked_Chr21_gene,
                                t1$interactor,
                                t1$tf),]
  PPI_nodes$mean_DELTA_act_tf = t1$mean_DELTA_act_tf[match(PPI_nodes$gene, t1$tf)]
  
  PPI_edges = Chr21_X_TF_interactions[(Chr21_X_TF_interactions$from %in% PPI_nodes$gene) & 
                                        (Chr21_X_TF_interactions$to %in% PPI_nodes$gene),]
  PPI_edges$cor = 0.5 #dummy correlation required for plotting function
  
  pl[[Chr21_gene]] = ppi_tf_act_plot(nodes = PPI_nodes, edges = PPI_edges, 
                                     plot_genes = NULL, 
                                     node_size_col = "mean_vst_CON", 
                                     node_color_col = "mean_z_DELTA",
                                     node_fill_col = "mean_DELTA_act_tf",
                                     node_size_limits = c(0, max(PPI_nodes$mean_vst_CON)),
                                     edge_alpha = 0.5, edgewidth_range = c(0.05, 0.5),
                                     title = paste0(Chr21_gene, " - TF interactions with interactors"))
  
}


# plot interactions for all Chr21 genes

pdf(file = paste0(out_dir,script_ind, "PPI_Network_plot_with_interactors_by Chr21_gene.pdf"), 
    width = 6, height = 9)
{
 lapply(pl, function(x){x}) 
}

dev.off()



# plot interactions for all Chr21 genes

pdf(file = paste0(out_dir,script_ind, "PPI_Network_plot_with_interactors_by Chr21_gene_larger.pdf"), 
    width = 5, height = 7)
{
  lapply(pl, function(x){x}) 
}

dev.off()












#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed H06 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


