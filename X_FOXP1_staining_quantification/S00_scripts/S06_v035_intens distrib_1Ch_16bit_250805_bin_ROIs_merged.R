### analysis single channel nuclear intensities by cell

# go to main directory (parent directory of scripts)
if (basename(getwd())== "S00_scripts"){setwd("../.")} 

message("Current analysis directory: ", getwd())

#load packages

library(tidyverse)
library(colorRamps)

#define working directories

in_dir = "./S05_R_input/"
out_dir = "./S06_intens_distrib_1Ch/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}


#### load data and set parameters 

load(file = paste0(in_dir,"/Raw_data_set_R.rda"))



### load image groups defined in csv sheet

group_tab_files = read_csv(paste0(in_dir,"/image_groups.csv"))

#add layer/bin to group_tab_files
t1 = gsub(".*layer_","",group_tab_files$file)
t2 = str_remove_all(t1, ".tif_results.csv")
group_tab_files$bin = paste0("bin",t2)

group_tab_files = group_tab_files[order(group_tab_files$bin),]


#add sample metadata

t1 =  read_csv("S00_input/group_tab.csv")

for (i in 1:nrow(t1)){
  group_tab_files$sample[grepl(t1$sample[i], group_tab_files$file)] = t1$sample[i]
}

t2 = group_tab_files[order(group_tab_files$bin, match(group_tab_files$repl, t1$sample)),]

t3 = left_join(t2, t1[,-1], by = "sample")


#group by bin and genotype

t3$genotype = t3$group.y
t3$group = t3$genotype
t3$repl = paste0(t3$group, "_", t3$sample)

group_tab_files = t3



#####################################################################
#  update and organise dataset   
#####################################################################

#define channel for analysis
Ch = set0$channels

#define analysis table order (Channels, Group order) from original table (int.tab)

#select groups and group order and color scheme
gr = unique(group_tab_files$group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}



### extract relevant channels and files, add metadata

replics = unique(group_tab_files$repl)

t1 = set0$intens_tab
t1 = t1[t1$file %in% group_tab_files$file,c("file", Ch)]
t2 = left_join(t1, group_tab_files, by = "file")


#reorder cells by group table order

t2 = t2[order(match(t2$file, group_tab_files$file)),]

#update ROI set, add group information
set1 = set0
set1$groups = gr
set1$replicates = replics
set1$files = group_tab_files$file
set1$group_tab = group_tab_files
set1$intens_tab = t2



#############################################
#calculate normalised intensities vs DAPI
#############################################

#normalise data by mean(DAPI(nuclear)) by file
t3 = set1$intens_tab

for (file1 in unique(t2$file)){
  
  t4 = t3[t3$file == file1,]
  
  mean_DAPI = mean(t4$DAPI)
  
  for (Ch1 in Ch){
    t4[[Ch1]] = t4[[Ch1]]/mean_DAPI
  }
  
  t3[t3$file == file1,] = t4
  
}


#scale each channel to max
for (Ch1 in Ch){
  t3[[Ch1]] = t3[[Ch1]]/max(t3[[Ch1]])
}


set1$intens_tab_DAPI_norm_by_file = t3



#############################################
#calculate normalised intensities scaled by sample/original image
#############################################

#normalise data by mean(DAPI(nuclear)) by file
t3 = set1$intens_tab

for (s1 in unique(group_tab_files$sample)){
  
  t4 = t3[t3$sample == s1,]
  
  for (Ch1 in Ch){
    v1 = t4[[Ch1]]
    quant = quantile(v1, probs = c(0.1, 0.2, 0.4))
    v2 = (v1-quant[2])/(quant[3]-quant[1])
    t4[[Ch1]] = v2
  }
  
  t3[t3$sample == s1,] = t4
  
}


set1$intens_tab_scaled_by_sample = t3


#save processed data for channel for further analysis
save(set1,
     file = paste(out_dir,"/dataset processed.RData",sep=""))





################################################################
### Plot intensity histogrammes by image/replicate/group for all channels
################################################################

#plots by group

t1 = as_tibble(set1$intens_tab)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 100)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(group, levels = gr), ncol = 1)+
    labs(title = paste0("Intensity distribution by group - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})


file.name = paste0(out_dir, "/Intensity_histogrammes_by_group.pdf")

pdf(file = file.name, width = 4, height = length(gr)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()


#plots by replicate

t1 = set1$intens_tab

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 100)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(repl, levels = replics), ncol = 1)+
    labs(title = paste0("Intensity distribution by replicate - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})



file.name = paste0(out_dir, "/Intensity_histogrammes_by_repl.pdf")

pdf(file = file.name, width = 4, height = length(replics)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()


#plots by file

t1 = set1$intens_tab

files = unique(t1$file)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 100)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(file, levels = files), ncol = 1)+
    labs(title = paste0("Intensity distribution by replicate - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})

file.name = paste0(out_dir, "/Intensity_histogrammes_by_file.pdf")

pdf(file = file.name, width = 4, height = length(files)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()






################################################################
# Plot intensity histogrammes by image/replicate/group for all channels (DAPI norm)
################################################################

#plots by group

t1 = as_tibble(set1$intens_tab_DAPI_norm_by_file)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 0.002)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(group, levels = gr), ncol = 1)+
    labs(title = paste0("Intensity distribution by group - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})


file.name = paste0(out_dir, "/Intensity_histogrammes_by_group_DAPI_norm.pdf")

pdf(file = file.name, width = 5, height = length(gr)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()


#plots by replicate

t1 = as_tibble(set1$intens_tab_DAPI_norm_by_file)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 0.01)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(repl, levels = replics), ncol = 1)+
    labs(title = paste0("Intensity distribution by replicate - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})



file.name = paste0(out_dir, "/Intensity_histogrammes_by_repl_DAPI_norm.pdf")

pdf(file = file.name, width = 5, height = length(replics)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()


#plots by file

t1 = as_tibble(set1$intens_tab_DAPI_norm_by_file)

files = unique(t1$file)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 0.01)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(file, levels = files), ncol = 1)+
    labs(title = paste0("Intensity distribution by replicate - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})

file.name = paste0(out_dir, "/Intensity_histogrammes_by_file_DAPI_norm.pdf")

pdf(file = file.name, width = 5, height = length(files)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()






################################################################
# Plot intensity histogrammes by image/replicate/group for all channels (scaled by sample/orig image)
################################################################

#plots by group

t1 = as_tibble(set1$intens_tab_scaled_by_sample)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 0.05)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(group, levels = gr), ncol = 1)+
    labs(title = paste0("Intensity distribution by group - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})


file.name = paste0(out_dir, "/Intensity_histogrammes_by_group_sample_scaled.pdf")

pdf(file = file.name, width = 5, height = length(gr)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()


#plots by replicate

t1 = as_tibble(set1$intens_tab_DAPI_norm_by_file)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 0.05)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(repl, levels = replics), ncol = 1)+
    labs(title = paste0("Intensity distribution by replicate - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})



file.name = paste0(out_dir, "/Intensity_histogrammes_by_repl_sample_scaled.pdf")

pdf(file = file.name, width = 5, height = length(replics)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()


#plots by file

t1 = as_tibble(set1$intens_tab_DAPI_norm_by_file)

files = unique(t1$file)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 0.05)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(file, levels = files), ncol = 1)+
    labs(title = paste0("Intensity distribution by replicate - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})

file.name = paste0(out_dir, "/Intensity_histogrammes_by_file_sample_scaled.pdf")

pdf(file = file.name, width = 5, height = length(files)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()





#plots by group (range -2 to +5)

t1 = as_tibble(set1$intens_tab_scaled_by_sample)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1)+geom_histogram(aes(x = .data[[Ch1]], fill = group), binwidth = 0.05)+
    theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(group, levels = gr), ncol = 1)+
    xlim(-2,5)+
    labs(title = paste0("Intensity distribution by group - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})


file.name = paste0(out_dir, "/Intensity_histogrammes_by_group_sample_scaled_zoomed.pdf")

pdf(file = file.name, width = 5, height = length(gr)*1.5+1)
lapply(pl_file, function(x){x})
dev.off()


