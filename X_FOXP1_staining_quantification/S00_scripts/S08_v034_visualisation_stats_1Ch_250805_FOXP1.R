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
out_dir = "./S08_visualisation_stats_1Ch_FOXP1/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}


#### load data and set parameters 

load(file = paste0(in_dir,"/dataset processed.RData"))

group_tab_samples = read_csv("S00_input/group_tab.csv")


#define channel and thresholds for analysis
set1$channels

Ch = "FOXP1"

thresh = 10000


#################################
#functions
#################################


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


#############################################
#create summary statistics for each threshold (norm intens)
#############################################

### add metadata to normalised intensity table

t1 = set1$intens_tab
t1$Ch = t1[,Ch]

names(t1) = str_replace(names(t1), "stage_PCW.x" , "stage_PCW")


int_norm_with_meta = t1

set1$int_norm_with_meta = int_norm_with_meta



###quantify cells per file

t1 = int_norm_with_meta

t2 = t1 %>% group_by(group, repl, file, genotype, stage_PCW) %>% summarise(N_cells = n(), 
                                                      N_pos = length(Ch[Ch>thresh]))

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

summary_by_file = t5

set1$summary_by_file = summary_by_file

write_csv(summary_by_file, file = paste0(out_dir, "/summary_by_file.csv"))




###quantify cells per replicate 

t1 = int_norm_with_meta

t2 = t1 %>% group_by(group, repl, genotype, sample, stage_PCW) %>% summarise(N_cells = n(), 
                                                                                      N_pos = length(Ch[Ch>thresh]))

t2$fract_pos = t2$N_pos/t2$N_cells


t3 = t2 %>% group_by(group) %>% summarise(N_cells_group_mean = mean(N_cells), 
                                          N_cells_group_sd = sd(N_cells), 
                                          fract_pos_group_mean = mean(fract_pos),
                                          fract_pos_group_sd = sd(fract_pos))

t4 = cbind(t2, t3[match(t2$group, t3$group),-1])

summary_by_repl = t4

set1$summary_by_repl = summary_by_repl

write_csv(summary_by_repl, file = paste0(out_dir, "/summary_by_repl.csv"))




##################################################################
#plot fraction pos cells in replicate by genotype
##################################################################

t2 = summary_by_repl

genotypes = unique(t2$genotype)

p1 = ggplot()+
  geom_col(data = t2, aes(x = genotype, y = fract_pos_group_mean, color = genotype), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.5)+
  geom_errorbar(data = t2, aes(x = genotype,
                               ymin = fract_pos_group_mean-fract_pos_group_sd,
                               y = fract_pos_group_mean,
                               ymax = fract_pos_group_mean+fract_pos_group_sd,
                               color =  genotype),
                position = position_dodge(width = 0.5), width = 0.3, lwd = 0.5)+
  geom_point(data = t2, aes(x = genotype, y = fract_pos, color = genotype), 
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), 
             size = 1, stroke = 0.5)+
  geom_hline(yintercept = 0)+
  scale_color_manual(limits = genotypes, values = pal(genotypes))+
  scale_fill_manual(limits = genotypes, values = pal(genotypes))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,"Fract_",Ch,"_pos_by_genotype.pdf"), 
    width = 2, height = 2)
plot(p1)
dev.off()




##################################################################
#statistical analysis t-tests
##################################################################

t1 = summary_by_repl

bins = unique(t1$bin)

t3 = pairwise.t.test(t1$fract_pos, t1$genotype, p.adjust.method = "none")

t4 = tibble(comp = t3$data.name, method = t3$method, p = t3$p.value[1,1])


write_csv(t4, paste0(out_dir,"Fract_",Ch,"_pos_by_genotype_t-Test.csv"))





#############################################################################
# plot fraction of positive cells by PCW and group
#############################################################################

t2 = summary_by_repl


p1 = ggplot(t2, aes(x = stage_PCW, y = fract_pos, color = group))+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm")+
  geom_point()+
  scale_color_manual(limits = genotypes, values = pal(genotypes))+
  labs(title = "Fraction of cells positive vs developmental stage (PCW)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,"fract_",Ch,"pos_by_repl_vs_PCW.pdf"), 
    width = 4, height = 3)
plot(p1)
dev.off()






