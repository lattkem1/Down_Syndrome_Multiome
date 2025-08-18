### analysis nuclear intensities single channel with different thresholds
#         with t-tests vs first group, and, if more images than replicates, lmerTest analysis

# go to main directory (parent directory of scripts)
if (basename(getwd())== "S00_scripts"){setwd("../.")} 

#load packages
library(tidyverse)
library(colorRamps)
library(lmerTest)

#define working directories

in_dir = "./S06_intens_distrib_1Ch/"
out_dir = "./S07_intens_thresh_analysis_1Ch/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}


#### load data and set parameters 

load(file = paste0(in_dir,"/dataset processed.RData"))

#define channel and thresholds for analysis
set1$channels

Ch = "FOXP1"

thresh_raw = c(10000, 15000, 20000)

thresh_norm = c(0.15, 0.25, 0.35)

thresh_scaled = c(0.15, 0.25, 0.35)



#############################################
#create summary statistics for each threshold (raw intens)
#############################################

t1 = as_tibble(set1$intens_tab)
t1$Ch = t1[,Ch]

thresh = thresh_raw

replics = unique(t1$repl)
gr = unique(t1$group)


#stats by file (with mean/sd by replicate and group)

stats_list_by_file = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t2$fract_pos = t2$N_pos/t2$N_cells
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells_repl_mean = mean(N_cells), 
                                                  N_cells_repl_sd = sd(N_cells), 
                                                  fract_pos_repl_mean = mean(fract_pos),
                                                  fract_pos_repl_sd = sd(fract_pos))
  
  t4 = t2 %>% group_by(group) %>% summarise(N_cells_group_mean = mean(N_cells), 
                                            N_cells_group_sd = sd(N_cells), 
                                            fract_pos_group_mean = mean(fract_pos),
                                            fract_pos_group_sd = sd(fract_pos))
  
  t5 = cbind(t2, t3[match(t2$repl, t3$repl),-c(1,2)], t4[match(t2$group, t4$group),-1])
  
  t5 = t5[order(match(t5$repl, replics)),]
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_image_", Ch, "_raw_int_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_file) = paste0(Ch, "_thresh_", thresh)



#stats by replicate (with mean/sd by group)

thresh1 = thresh[1]

stats_list_by_repl = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells = mean(N_cells), N_pos = mean(N_pos))
  t3$fract_pos = t3$N_pos/t3$N_cells
  
  t4 = t3 %>% group_by(group) %>% summarise(N_cells_mean = mean(N_cells), 
                                            N_cells_sd = sd(N_cells), 
                                            fract_pos_mean = mean(fract_pos),
                                            fract_pos_sd = sd(fract_pos) )
  
  t5 = cbind(t3, t4[match(t3$group, t4$group),-1])
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_repl_", Ch, "_raw_int_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_repl) = paste0(Ch, "_thresh_", thresh)





#############################################
#plot summary statistics for each threshold (raw intens)
#############################################

thresh = thresh_raw

#define plot colors (by group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}



### plot by replicate (individual images)

pl_repl = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = repl, y = fract_pos_repl_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = repl, ymin = fract_pos_repl_mean-fract_pos_repl_sd,
                      y = fract_pos_repl_mean,
                      ymax =  fract_pos_repl_mean+fract_pos_repl_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = repl, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = replics)+
    geom_hline(yintercept = 0)+
    labs(title = paste0("Fraction ", Ch, " pos (thresh = ", thresh1,"; individual images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_repl) = paste0(Ch, "_thresh_", thresh)


# add plot with cell numbers by replicate (individual images)

t1 = stats_list_by_file[[1]]

pl_repl[["N_cells"]] = ggplot(data = t1)+
  geom_col(aes(x = repl, y = N_cells_repl_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = repl, ymin = N_cells_repl_mean-N_cells_repl_sd,
                    y = N_cells_repl_mean,
                    ymax =  N_cells_repl_mean+N_cells_repl_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = repl, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = replics)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#save plots by replicate

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_replicate_raw_int.pdf")

pdf(file = file.name, width = length(replics)/4+2, height = 3)
lapply(pl_repl, function(x){x})
dev.off()




### plots by group (individual images/replicates)

# threshold plots individual images by group

pl_gr = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_group_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_group_mean-fract_pos_group_sd,
                      y = fract_pos_group_mean,
                      ymax =  fract_pos_group_mean+fract_pos_group_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0( Ch, " (thresh = ", thresh1,"; indiv images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_gr) = paste0(Ch, "_thresh_", thresh, "_indiv_images")


# add plot with cell numbers individual images by group

t1 = stats_list_by_file[[1]]

pl_gr[["N_cells_indiv_images"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_group_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_group_mean-N_cells_group_sd,
                    y = N_cells_group_mean,
                    ymax =  N_cells_group_mean+N_cells_group_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



# create threshold plots individual replicates by group

l1 = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_repl[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_mean-fract_pos_sd,
                      y = fract_pos_mean,
                      ymax =  fract_pos_mean+fract_pos_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0(Ch, " (thresh = ", thresh1,"; indiv replicates)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(l1) = paste0(Ch, "_thresh_", thresh, "_indiv_replicates")


# add plot with cell numbers individual replicates by group

t1 = stats_list_by_repl[[1]]

l1[["N_cells_indiv_repl"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_mean-N_cells_sd,
                    y = N_cells_mean,
                    ymax =  N_cells_mean+N_cells_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual replicates)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#merge plots of individual images and replicates and save plots by group
pl_gr = c(pl_gr, l1)

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_group_raw_int.pdf")

pdf(file = file.name, width = length(gr)/4+2, height = 3)
lapply(pl_gr, function(x){x})
dev.off()





#############################################
#create summary statistics for each threshold (DAPI norm intens)
#############################################

t1 = as_tibble(set1$intens_tab_DAPI_norm_by_file)
t1$Ch = t1[,Ch]

thresh = thresh_norm

replics = unique(t1$repl)
gr = unique(t1$group)


#stats by file (with mean/sd by replicate and group)

stats_list_by_file = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t2$fract_pos = t2$N_pos/t2$N_cells
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells_repl_mean = mean(N_cells), 
                                                  N_cells_repl_sd = sd(N_cells), 
                                                  fract_pos_repl_mean = mean(fract_pos),
                                                  fract_pos_repl_sd = sd(fract_pos))
  
  t4 = t2 %>% group_by(group) %>% summarise(N_cells_group_mean = mean(N_cells), 
                                            N_cells_group_sd = sd(N_cells), 
                                            fract_pos_group_mean = mean(fract_pos),
                                            fract_pos_group_sd = sd(fract_pos))
  
  t5 = cbind(t2, t3[match(t2$repl, t3$repl),-c(1,2)], t4[match(t2$group, t4$group),-1])
  
  t5 = t5[order(match(t5$repl, replics)),]
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_image_", Ch, "_DAPI_norm_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_file) = paste0(Ch, "_thresh_", thresh)



#stats by replicate (with mean/sd by group)

thresh1 = thresh[1]

stats_list_by_repl = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells = mean(N_cells), N_pos = mean(N_pos))
  t3$fract_pos = t3$N_pos/t3$N_cells
  
  t4 = t3 %>% group_by(group) %>% summarise(N_cells_mean = mean(N_cells), 
                                            N_cells_sd = sd(N_cells), 
                                            fract_pos_mean = mean(fract_pos),
                                            fract_pos_sd = sd(fract_pos) )
  
  t5 = cbind(t3, t4[match(t3$group, t4$group),-1])
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_repl_", Ch, "_DAPI_norm_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_repl) = paste0(Ch, "_thresh_", thresh)





#############################################
#plot summary statistics for each threshold (DAPI norm intens)
#############################################

thresh = thresh_norm

#define plot colors (by group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}



### plot by replicate (individual images)

pl_repl = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = repl, y = fract_pos_repl_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = repl, ymin = fract_pos_repl_mean-fract_pos_repl_sd,
                      y = fract_pos_repl_mean,
                      ymax =  fract_pos_repl_mean+fract_pos_repl_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = repl, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = replics)+
    geom_hline(yintercept = 0)+
    labs(title = paste0("Fraction ", Ch, " pos (thresh = ", thresh1,"; individual images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_repl) = paste0(Ch, "_thresh_", thresh)


# add plot with cell numbers by replicate (individual images)

t1 = stats_list_by_file[[1]]

pl_repl[["N_cells"]] = ggplot(data = t1)+
  geom_col(aes(x = repl, y = N_cells_repl_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = repl, ymin = N_cells_repl_mean-N_cells_repl_sd,
                    y = N_cells_repl_mean,
                    ymax =  N_cells_repl_mean+N_cells_repl_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = repl, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = replics)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#save plots by replicate

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_replicate_DAPI_norm.pdf")

pdf(file = file.name, width = length(replics)/4+2, height = 3)
lapply(pl_repl, function(x){x})
dev.off()




### plots by group (individual images/replicates)

# threshold plots individual images by group

pl_gr = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_group_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_group_mean-fract_pos_group_sd,
                      y = fract_pos_group_mean,
                      ymax =  fract_pos_group_mean+fract_pos_group_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0( Ch, " (thresh = ", thresh1,"; indiv images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_gr) = paste0(Ch, "_thresh_", thresh, "_indiv_images")


# add plot with cell numbers individual images by group

t1 = stats_list_by_file[[1]]

pl_gr[["N_cells_indiv_images"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_group_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_group_mean-N_cells_group_sd,
                    y = N_cells_group_mean,
                    ymax =  N_cells_group_mean+N_cells_group_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



# create threshold plots individual replicates by group

l1 = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_repl[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_mean-fract_pos_sd,
                      y = fract_pos_mean,
                      ymax =  fract_pos_mean+fract_pos_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0(Ch, " (thresh = ", thresh1,"; indiv replicates)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(l1) = paste0(Ch, "_thresh_", thresh, "_indiv_replicates")


# add plot with cell numbers individual replicates by group

t1 = stats_list_by_repl[[1]]

l1[["N_cells_indiv_repl"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_mean-N_cells_sd,
                    y = N_cells_mean,
                    ymax =  N_cells_mean+N_cells_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual replicates)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#merge plots of individual images and replicates and save plots by group
pl_gr = c(pl_gr, l1)

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_group_DAPI_norm.pdf")

pdf(file = file.name, width = length(gr)/4+2, height = 3)
lapply(pl_gr, function(x){x})
dev.off()







#############################################
#create summary statistics for each threshold (scaled intens by sample/orig image)
#############################################

t1 = as_tibble(set1$intens_tab_scaled_by_sample)
t1$Ch = t1[,Ch]

thresh = thresh_scaled

replics = unique(t1$repl)
gr = unique(t1$group)


#stats by file (with mean/sd by replicate and group)

stats_list_by_file = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t2$fract_pos = t2$N_pos/t2$N_cells
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells_repl_mean = mean(N_cells), 
                                                  N_cells_repl_sd = sd(N_cells), 
                                                  fract_pos_repl_mean = mean(fract_pos),
                                                  fract_pos_repl_sd = sd(fract_pos))
  
  t4 = t2 %>% group_by(group) %>% summarise(N_cells_group_mean = mean(N_cells), 
                                            N_cells_group_sd = sd(N_cells), 
                                            fract_pos_group_mean = mean(fract_pos),
                                            fract_pos_group_sd = sd(fract_pos))
  
  t5 = cbind(t2, t3[match(t2$repl, t3$repl),-c(1,2)], t4[match(t2$group, t4$group),-1])
  
  t5 = t5[order(match(t5$repl, replics)),]
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_image_", Ch, "_scaled_by_sample_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_file) = paste0(Ch, "_thresh_", thresh)



#stats by replicate (with mean/sd by group)

thresh1 = thresh[1]

stats_list_by_repl = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells = mean(N_cells), N_pos = mean(N_pos))
  t3$fract_pos = t3$N_pos/t3$N_cells
  
  t4 = t3 %>% group_by(group) %>% summarise(N_cells_mean = mean(N_cells), 
                                            N_cells_sd = sd(N_cells), 
                                            fract_pos_mean = mean(fract_pos),
                                            fract_pos_sd = sd(fract_pos) )
  
  t5 = cbind(t3, t4[match(t3$group, t4$group),-1])
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_repl_", Ch, "_scaled_by_sample_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_repl) = paste0(Ch, "_thresh_", thresh)





#############################################
#plot summary statistics for each threshold (scaled intens by sample/orig image)
#############################################

thresh = thresh_scaled

#define plot colors (by group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}



### plot by replicate (individual images)

pl_repl = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = repl, y = fract_pos_repl_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = repl, ymin = fract_pos_repl_mean-fract_pos_repl_sd,
                      y = fract_pos_repl_mean,
                      ymax =  fract_pos_repl_mean+fract_pos_repl_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = repl, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = replics)+
    geom_hline(yintercept = 0)+
    labs(title = paste0("Fraction ", Ch, " pos (thresh = ", thresh1,"; individual images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_repl) = paste0(Ch, "_thresh_", thresh)


# add plot with cell numbers by replicate (individual images)

t1 = stats_list_by_file[[1]]

pl_repl[["N_cells"]] = ggplot(data = t1)+
  geom_col(aes(x = repl, y = N_cells_repl_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = repl, ymin = N_cells_repl_mean-N_cells_repl_sd,
                    y = N_cells_repl_mean,
                    ymax =  N_cells_repl_mean+N_cells_repl_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = repl, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = replics)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#save plots by replicate

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_replicate_scaled_by_sample.pdf")

pdf(file = file.name, width = length(replics)/4+2, height = 3)
lapply(pl_repl, function(x){x})
dev.off()




### plots by group (individual images/replicates)

# threshold plots individual images by group

pl_gr = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_group_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_group_mean-fract_pos_group_sd,
                      y = fract_pos_group_mean,
                      ymax =  fract_pos_group_mean+fract_pos_group_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0( Ch, " (thresh = ", thresh1,"; indiv images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_gr) = paste0(Ch, "_thresh_", thresh, "_indiv_images")


# add plot with cell numbers individual images by group

t1 = stats_list_by_file[[1]]

pl_gr[["N_cells_indiv_images"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_group_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_group_mean-N_cells_group_sd,
                    y = N_cells_group_mean,
                    ymax =  N_cells_group_mean+N_cells_group_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



# create threshold plots individual replicates by group

l1 = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_repl[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_mean-fract_pos_sd,
                      y = fract_pos_mean,
                      ymax =  fract_pos_mean+fract_pos_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0(Ch, " (thresh = ", thresh1,"; indiv replicates)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(l1) = paste0(Ch, "_thresh_", thresh, "_indiv_replicates")


# add plot with cell numbers individual replicates by group

t1 = stats_list_by_repl[[1]]

l1[["N_cells_indiv_repl"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_mean-N_cells_sd,
                    y = N_cells_mean,
                    ymax =  N_cells_mean+N_cells_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual replicates)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#merge plots of individual images and replicates and save plots by group
pl_gr = c(pl_gr, l1)

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_group_scaled_by_sample.pdf")

pdf(file = file.name, width = length(gr)/4+2, height = 3)
lapply(pl_gr, function(x){x})
dev.off()





