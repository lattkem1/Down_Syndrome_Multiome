message("\n\n##########################################################################\n",
        "# Start E04: DEG analysis overlap with other datasets  ", Sys.time(),
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
library(VennDiagram)
library(colorRamps)
library(viridis)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(msigdbr)
library(ggrepel)

#specify script/output index as prefix for file names
script_ind = "E03a_"


#specify output directory
out_dir = paste0(main_dir,"E_DESeq_pseudobulk_by_cluster_exc_lin_from_all_non_cx_excl_v045/")


#load DESeq2 results 

load(file = paste0(out_dir, "E03_bulk_data_w_expr_z_scores.rda"))


#load Nebula dataset for comparison

load(file = paste0(out_dir,"E02a_nebula_data.rda")) 



### get genes of interest panels

GOI = list()
t1 = read_csv(paste0(main_dir,"A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol
t1 = read_csv(paste0(main_dir,"A_input/HSA21_genes_biomaRt_conversion.csv"))
GOI$Chr21 = t1$hgnc_symbol
t1 = read_tsv(paste0(main_dir,"A_input/Genomics_EnglandPanelApp_Intellectual_disability_v8.243.tsv"))
GOI$ID_genes = unique(t1$`Gene Symbol`)

t1 = bulk_data$GO_results$full
l1 = str_split(t1$geneID[t1$Description == "forebrain development"], "/")
GOI$GO_forebrain_development = unlist(l1)
l1 = str_split(t1$geneID[t1$Description == "axonogenesis"], "/")
GOI$GO_axonogenesis = unlist(l1)


### load published DEG sets

pub_DEG_list = list()

t1 = read_csv(paste0(main_dir,"A_input/DEGs_publ_data/Rastogi24_TableS2_DEGs_Cor_genes.csv"))
t2 = t1[t1$FDR<0.1,]
pub_DEG_list[["Rastogi24_DEGs_Cor_up"]] = unique(t2$Symbol[t2$logFC>0])
pub_DEG_list[["Rastogi24_DEGs_Cor_down"]] = unique(t2$Symbol[t2$logFC<0])



t1 = read_csv(paste0(main_dir,"A_input/DEGs_publ_data/Rastogi24_TableS2_DEGs_iNeu_genes.csv"))
t2 = t1[t1$FDR<0.1,]
pub_DEG_list[["Rastogi24_DEGs_iNeu_genes_up"]] = unique(t2$SYMBOL[t2$logFC>0])
pub_DEG_list[["Rastogi24_DEGs_iNeu_genes_down"]] = unique(t2$SYMBOL[t2$logFC<0])


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



###function: plot bulk gene expression heatmap with annotations

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
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
  
  p1 = pheatmap::pheatmap(pl_mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                          show_rownames = show_rownames, show_colnames = show_colnames,
                          color = color,
                          breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                          annotation_col = annot_col, annotation_colors = annot_colors,
                          border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                          fontsize = fontsize, main = title
  )
  
  return(p1)
  
}



###########################################################
# plot volcano plots with selected gene sets highlighted
###########################################################

pl = list()

for (comp in names(bulk_data$deseq_results)){
  
  
  t1 = bulk_data$deseq_results[[comp]]
  t1$Chr21 = rownames(t1) %in% GOI$Chr21
  t1$neuro_dev = rownames(t1) %in% c(GOI$ID_genes, GOI$GO_forebrain_development, GOI$GO_axonogenesis)
  t1$DEG = !is.na(t1$padj) & t1$padj <= 0.1 & !is.na(t1$log2FoldChange) & abs(t1$log2FoldChange)>= log2(1.2)
  
  t1$gene_type = "other_genes"
  t1$gene_type[t1$DEG] = "other_DEGs"
  t1$gene_type[t1$neuro_dev & t1$DEG] = "neurodevelopmental_DEGs"
  t1$gene_type[t1$Chr21] = "Chr21_genes"
  
  t1$plot_label = ""
  t1$plot_label[t1$gene_type == "neurodevelopmental_DEGs"|(t1$Chr21 & t1$DEG)] = 
    rownames(t1)[t1$gene_type == "neurodevelopmental_DEGs"|(t1$Chr21 & t1$DEG)]
  
  pl[[comp]] = ggplot(t1, aes(x = log2FoldChange, y = -log10(padj), color = gene_type))+
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_hline(yintercept = -log10(0.1), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_point(aes(size = gene_type), alpha = 0.8)+
    geom_text_repel((aes(label = plot_label)))+
    scale_size_manual(limits = c("Chr21_genes","neurodevelopmental_DEGs", "other_DEGs", "other_genes"), 
                      values = c(2, 2, 1, 1))+
    scale_color_manual(limits = c("Chr21_genes","neurodevelopmental_DEGs", "other_DEGs", "other_genes"), 
                       values = c("red", "purple","black", "grey"))+
    theme_minimal()+
    labs(title = comp)
  
}


pdf(file = paste0(out_dir,script_ind, "Volcano_plots_DEGs_by_cluster.pdf"), 
    width = 16, height = 12)
{
  lapply(pl, function(x){x})
}
dev.off()




###########################################################
# stats overlaps between deseq DEGs and nebula DEGs 
###########################################################

l1 = bulk_data$DEGs
l2 = nebula_data$DEGs

t1 = tibble(comps = names(l1), N_DEGs_deseq = lengths(l1), N_DEGs_nebula = lengths(l2[match(names(l1), names(l2))]))

for (i in 1:nrow(t1)){
  v1 = intersect(l1[[t1$comps[i] ]], l2[[t1$comps[i] ]])
  t1$N_DEGs_deseq_in_nebula[i] = length(v1)
}

write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_deseq2_vs_nebula.csv"))



###########################################################
# identify overlaps between DESeq2 DEGs and Nebula DEGs (overlap matrices)
###########################################################

overlap_mat = matrix(nrow = length(nebula_data$DEGs), ncol = length(bulk_data$DEGs),
                     dimnames = list(names(nebula_data$DEGs), names(bulk_data$DEGs)))

overlap_mat_N = overlap_mat
overlap_mat_fract = overlap_mat

for (DEG_sets_nebula in names(nebula_data$DEGs)){
  
  for (DEG_sets_deseq in names(bulk_data$DEGs)){
    
    DEGs_nebula = nebula_data$DEGs[[DEG_sets_nebula]]
    DEGs_deseq = bulk_data$DEGs[[DEG_sets_deseq]]
    
    overlap_mat_N[DEG_sets_nebula, DEG_sets_deseq] = length(intersect(DEGs_nebula, DEGs_deseq))
    overlap_mat_fract[DEG_sets_nebula, DEG_sets_deseq] = length(intersect(DEGs_nebula, DEGs_deseq))/
      length(DEGs_deseq)
  }
}

dimnames(overlap_mat_N) = list(paste0(rownames(overlap_mat), " (Nebula)"), 
                               paste0(colnames(overlap_mat), " (DESeq2)"))
dimnames(overlap_mat_fract) = dimnames(overlap_mat_N)
                               
t1 = as_tibble(cbind(comp = rownames(overlap_mat_N), overlap_mat_N))
write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_cluster_deseq_val_by_nebula.csv"))

t1 = as_tibble(cbind(comp = rownames(overlap_mat_fract), overlap_mat_fract))
write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_cluster_deseq_val_by_nebula_fract.csv"))



###plot overlap stats

pdf(file = paste0(out_dir,script_ind, "Overlap_tissue_DEGs_by_cluster_deseq_val_by_nebula.pdf"), 
    width = 10, height = 10)
{
  
  pheatmap::pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(overlap_mat_N), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Number of DEGs (DESeq2) by cluster recalled with Nebula")
  
  pheatmap::pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = TRUE, cluster_cols = TRUE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(overlap_mat_N), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Number of DEGs (DESeq2) by cluster recalled with Nebula")
  
  
  pheatmap::pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(overlap_mat_fract), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Fraction of DEGs (DESeq2) by cluster recalled with Nebula")
  
  pheatmap::pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = TRUE, cluster_cols = TRUE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(overlap_mat_fract), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Fraction of DEGs (DESeq2) by cluster recalled with Nebula")
  
}

dev.off()



###########################################################
# identify overlaps between tissue DEGs and published DEGs (table genes vs comparisons)
###########################################################

l1 = bulk_data$DEGs

gene_comp_tab = tibble(gene = unique(unlist(l1))) 
t1 = gene_comp_tab
t1$x = ""

###identify genes by comparison for DESeq2 analysis 

for (comp in names(l1)){
  
  v1 = l1[[comp]]
  t1$x = ""
  t1$x[match(v1, t1$gene)] = v1
  
  gene_comp_tab[[comp]] = t1$x
  
}


###identify genes by comparison for Nebula analysis 

l1 = nebula_data$DEGs

for (comp in names(l1)){
  
  v1 = intersect(l1[[comp]], t1$gene)
  t1$x = ""
  t1$x[match(v1, t1$gene)] = v1
  
  gene_comp_tab[[paste0("nebula_", comp)]] = t1$x
  
}

###identify genes by comparison for published data 

l1 = pub_DEG_list

for (comp in names(l1)){
  
  v1 = intersect(l1[[comp]], t1$gene)
  t1$x = ""
  t1$x[match(v1, t1$gene)] = v1
  
  gene_comp_tab[[comp]] = t1$x
  
}

### reorder by up vs downregulated genes

t1 = gene_comp_tab
gene_comp_tab = cbind(t1[,"gene"], t1[,grepl("_up", colnames(t1))],t1[,grepl("_down", colnames(t1))])

write_csv(gene_comp_tab, file = paste0(out_dir, script_ind, "DEGs_deseq2_vs_nebula_publ_datasets.csv"))



###########################################################
# identify/plot overlaps between tissue DEGs and published DEGs (DESeq2 DEGs)
###########################################################

overlap_mat = matrix(nrow = length(bulk_data$DEGs), ncol = length(pub_DEG_list),
                     dimnames = list(names(bulk_data$DEGs), names(pub_DEG_list)))

overlap_mat_N = overlap_mat
overlap_mat_fract = overlap_mat

for (DEG_sets_deseq in names(bulk_data$DEGs)){
  
  for (DEG_sets_publ in names(pub_DEG_list)){
    
    DEGs_deseq = bulk_data$DEGs[[DEG_sets_deseq]]
    DEGs_publ = pub_DEG_list[[DEG_sets_publ]]
    
    overlap_mat_N[DEG_sets_deseq, DEG_sets_publ] = length(intersect(DEGs_deseq, DEGs_publ))
    overlap_mat_fract[DEG_sets_deseq, DEG_sets_publ] = length(intersect(DEGs_deseq, DEGs_publ))/
      length(DEGs_deseq)
    
  }
}

t1 = as_tibble(cbind(comp = rownames(overlap_mat_N), overlap_mat_N))
write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_cluster_vs_published_datasets.csv"))

t1 = as_tibble(cbind(comp = rownames(overlap_mat_fract), overlap_mat_fract))
write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_cluster_vs_published_datasets_fract.csv"))



###plot overlap stats

pdf(file = paste0(out_dir,script_ind, "Overlap_tissue_DEGs_by_cluster_vs_publ_datasets.pdf"), 
    width = 10, height = 10)
{
  
  pheatmap::pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(overlap_mat_N), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Number of DEGs by cluster vs reg in published data")
  
  pheatmap::pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, max(overlap_mat_N), length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Number of DEGs by cluster vs reg in published data")
  
  
  pheatmap::pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, 1, length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Fraction of DEGs by cluster vs reg in published data")
  
  pheatmap::pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,  
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 10, treeheight_col = 10,
           color = colorRampPalette(c("white", "blue"))(250),
           breaks = seq(0, 1, length.out = 251),
           border_color = NA, fontsize = 10,
           cellwidth = 10, cellheight = 10,
           main = "Fraction of DEGs by cluster vs reg in published data")
  
  #plot ordered by up vs down
  
  m1 = overlap_mat_fract
  m1 = m1[order(grepl("up", rownames(m1))) ,order(grepl("up", colnames(m1)))]
  
  pheatmap::pheatmap(m1, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, 1, length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Fraction of DEGs by cluster vs reg in published data - ordered up/down")
  
}

dev.off()




###########################################################
# identify/plot overlaps between tissue DEGs and published DEGs (Nebula DEGs)
###########################################################

overlap_mat = matrix(nrow = length(nebula_data$DEGs), ncol = length(pub_DEG_list),
                     dimnames = list(names(nebula_data$DEGs), names(pub_DEG_list)))

overlap_mat_N = overlap_mat
overlap_mat_fract = overlap_mat

for (DEG_sets_deseq in names(nebula_data$DEGs)){
  
  for (DEG_sets_publ in names(pub_DEG_list)){
    
    DEGs_deseq = nebula_data$DEGs[[DEG_sets_deseq]]
    DEGs_publ = pub_DEG_list[[DEG_sets_publ]]
    
    overlap_mat_N[DEG_sets_deseq, DEG_sets_publ] = length(intersect(DEGs_deseq, DEGs_publ))
    overlap_mat_fract[DEG_sets_deseq, DEG_sets_publ] = length(intersect(DEGs_deseq, DEGs_publ))/
      length(DEGs_deseq)
    
  }
}

t1 = as_tibble(cbind(comp = rownames(overlap_mat_N), overlap_mat_N))
write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_Nebula_cluster_vs_published_datasets.csv"))

t1 = as_tibble(cbind(comp = rownames(overlap_mat_fract), overlap_mat_fract))
write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_Nebula_cluster_vs_published_datasets_fract.csv"))



###plot overlap stats

pdf(file = paste0(out_dir,script_ind, "Overlap_tissue_DEGs_Nebula_by_cluster_vs_publ_datasets.pdf"), 
    width = 10, height = 10)
{
  
  pheatmap::pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(overlap_mat_N), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Number of DEGs by cluster vs reg in published data")
  
  pheatmap::pheatmap(overlap_mat_N, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = TRUE, cluster_cols = TRUE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(overlap_mat_N), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Number of DEGs by cluster vs reg in published data")
  
  
  pheatmap::pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, 1, length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Fraction of DEGs by cluster vs reg in published data")
  
  pheatmap::pheatmap(overlap_mat_fract, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = TRUE, cluster_cols = TRUE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, 1, length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Fraction of DEGs by cluster vs reg in published data")
  
  #plot ordered by up vs down
  
  m1 = overlap_mat_fract
  m1 = m1[order(grepl("up", rownames(m1))) ,order(grepl("up", colnames(m1)))]
  
  pheatmap::pheatmap(m1, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, 1, length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = "Fraction of DEGs by cluster vs reg in published data - ordered up/down")
  
}

dev.off()



#########################################
#plot genes validated in any published dataset by GO_term (max top20 GO terms) (mean DELTA CON)
#########################################

###identify genes regulated in same direction in DESeq2 analysis and any published dataset
l1 = bulk_data$DEGs
l2 = l1[grepl("_up", names(l1))]
for (comp in names(l2)){l2[[comp]] = intersect(l2[[comp]], 
                                               unlist(pub_DEG_list[grepl("_up", names(pub_DEG_list))]))}
l3 = l1[grepl("_down", names(l1))]
for (comp in names(l3)){l3[[comp]] = intersect(l3[[comp]], 
                                               unlist(pub_DEG_list[grepl("_down", names(pub_DEG_list))]))}
DEGs_validated = unique(unlist(c(l2, l3)))

t1 = tibble(N_DEGs_deseq_all = length(unique(unlist(l1))), 
            N_DEGs_deseq_reg_same_dir_any_dataset = length(DEGs_validated))

write_csv(t1, file = paste0(out_dir, script_ind, "N_DEGs_validated_in_any_published_dataset.csv"))


###identify overlap with GO genes

go_res = bulk_data$GO_results$full
if (nrow(go_res) > 20){go_res = go_res[1:20,]}

go_genes_list = str_split(go_res$geneID, "/")
names(go_genes_list) = paste0(go_res$Description," (", go_res$ID,")")


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_Nebula_DEGs_validated_z_score_mean_CON_DELTA.pdf"), 
    width = 12, height = 20)
{
  for (go in names(go_genes_list)){
    
    pl_genes = intersect(go_genes_list[[go]], DEGs_validated)
    
    pl_mat_CON = bulk_data$gene_Z_scores_cluster$CON[pl_genes,]
    lims_CON = 0.7*c(-max(abs(pl_mat_CON)), max(abs(pl_mat_CON)))
    
    pl_mat_DELTA = bulk_data$gene_Z_scores_cluster$DELTA[pl_genes,]
    lims_DELTA = 0.4*c(-max(abs(pl_mat_DELTA)), max(abs(pl_mat_DELTA)))
    
    if (length(pl_genes)>1){
      
      p1 = bulkdata_heatmap(pl_mat = pl_mat_DELTA, 
                            pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                            pl_genes = pl_genes,
                            x_col = "cluster_name", 
                            meta_annot_cols = c("cluster_name"),
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = lims_DELTA,  cellwidth = 10, cellheight = 10, 
                            fontsize = 10, title = paste0(go, " - DEGs - Expr DS vs CON (cluster Z-score)"))
      
      pl_genes_clust = pl_genes[p1$tree_row$order]
      
      bulkdata_heatmap(pl_mat = pl_mat_CON, 
                       pl_meta = bulk_data$gene_Z_scores_cluster$meta,
                       pl_genes = pl_genes_clust,
                       x_col = "cluster_name", 
                       meta_annot_cols = c("cluster_name"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = TRUE, show_colnames = TRUE,
                       color = viridis(250),
                       lims = lims_CON,  cellwidth = 10, cellheight = 10, 
                       fontsize = 10, title = paste0(go, " - DEGs - Expr CON (cluster Z-score)"))
    }
  }
  
}
dev.off()











#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed E02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


