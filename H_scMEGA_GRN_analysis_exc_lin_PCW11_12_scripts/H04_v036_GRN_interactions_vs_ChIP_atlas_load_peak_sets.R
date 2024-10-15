message("\n\n##########################################################################\n",
        "# Start H04: scMEGA network interactions vs published ChIP peaks: Load peak sets ", Sys.time(),
        "\n##########################################################################\n",
        "\n create/update locally saved collection of bed files from ChIP Atlas database ",
        "\n check locally saved collection, download missing ChIP sets for GRN TFs",
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



#specify script/output index as prefix for file names
script_ind = "H04_"

#specify output directory
out_dir = paste0(main_dir,"H_scMEGA_GRN_analysis_exc_lin_PCW11_12/")


#load network TFs
net_nodes = read_csv(file = paste0(out_dir, "H03_network_nodes.csv"))
net_nodes = net_nodes[net_nodes$DS_reg,] 

#load previously downloaded ChIP peak sets
ChIP_atlas_path = "/rds/general/user/mlattke/projects/dsseq23/live/ChIP_Atlas_peak_bed_list.rda"

if (file.exists(ChIP_atlas_path)){load(ChIP_atlas_path)}else{ChIP_peaks_list = list()}


###specify TFs for which to download peak bed files from ChIP-Atlas (from all cell types)

#select TFs for validation 
net_TFs = unique(net_nodes$gene[net_nodes$N_targets>0])
TFs_to_load = net_TFs[!(net_TFs %in% names(ChIP_peaks_list))]


###############################################
# load bed files for TF ChIP from ChIP-Atlas 
################################################

#ChIP_peaks_list = list()

if(length(TFs_to_load)>0){
  
  for (i in 1:length(TFs_to_load)){
    tf = TFs_to_load[i]
    message("load ChIP peaks for ", tf, " (TF ", i, " of ", length(TFs_to_load), ")")
    
    file1 = paste0("https://chip-atlas.dbcls.jp/data/hg38/assembled/Oth.ALL.05.",
                   tf, ".AllCell.bed")
    
    if (url.exists(file1)){
      
      t1 = read_lines(file1)
      t2 = t1[-1]
      t2 = str_split(t2, "\t")
      ChIP_peaks = as.data.frame(do.call(rbind, t2))
      colnames(ChIP_peaks) = c("seqnames", "start", "end")
      ChIP_peaks_list[[tf]] = ChIP_peaks
      
    } else {ChIP_peaks_list[[tf]] = NULL}
    
    ChIP_peaks_list[[tf]] = ChIP_peaks
    
  }
  
  save(ChIP_peaks_list, file = ChIP_atlas_path)
  
}


#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed H04 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


