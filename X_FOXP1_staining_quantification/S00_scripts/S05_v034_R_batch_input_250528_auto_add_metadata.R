###############################################################
# S05: Load FIJI results into R
###############################################################
#
#
###############################################################


# go to main directory (parent directory of scripts)
if (basename(getwd())== "S00_scripts"){setwd("../.")} 

message("Current analysis directory: ", getwd())

#load packages

library("tidyverse")


#define working directories
main_dir = getwd()
exp.name = basename(getwd())
in_dir = "./S04_FIJI_results_nucl_intens"
out_dir = "./S05_R_input"
if (!dir.exists(out_dir)) {dir.create(out_dir)}
setwd(main_dir)

### specify data

#input filelist
flist = list.files(in_dir)

#define channel names
Ch.names = c("DAPI", "SATB2", "FOXP1", "CTIP2")


#load sample/raw image metadata
gr_tab_samples = read_csv("S00_input/group_tab.csv")



#########################
# load data
#########################

########### read in data

int.tab = NULL

#perform for each file (row in group.tab)
i=1
  for (i in 1:length(flist)){
    
    #read file and define ordering parameters
    fname = paste(in_dir,"/",flist[i], sep ="")
    rawdata = read_csv(fname)
    
    #generate table of intensities for each cell for each channel, add group, replicate, file, color in first columns
    int.tab.fl = filter(rawdata, rawdata[,"Ch"] == 1)[,"Mean"]
    if (length(Ch.names) > 1){
      for (h in 2:length(Ch.names)){
        int.tab.fl = cbind(int.tab.fl, filter(rawdata, rawdata[,"Ch"] == h)[,"Mean"])
      }
    }
    int.tab.fl = cbind(group = NA, repl = NA, file = flist[i], color = NA, int.tab.fl)
    
    #add to table with all cells of all files
    int.tab = rbind(int.tab, int.tab.fl)
  }


names(int.tab)[-c(1:4)] = Ch.names

#convert file names from factor to character (else problems with downstream analysis) 
int.tab$file = as.character(int.tab$file) 

#create table of files with group assignment, and list of groups
gr_tab_files = tibble(file = flist, group ="", repl="")


#create collect data in ROI set
set0 = list(channels = Ch.names, 
            groups = NULL, 
            replicates = NULL,
            files = flist,
            group_tab = gr_tab_files,
            intens_tab = int.tab
            )
class(set0) = "ROI_set"


#save processed input data for further analysis
save(set0, file = paste0(out_dir,"/Raw_data_set_R.rda"))


####################################
# add sample metadata
####################################

gr_tab_samples$image = str_remove_all(gr_tab_samples$image, ".tif")

for (i in 1:nrow(gr_tab_samples)){
  
  gr_tab_files$repl[grepl(gr_tab_samples$image[i], gr_tab_files$file)] = gr_tab_samples$sample[i] 
  
}

for (meta1 in colnames(gr_tab_samples)){
  
  gr_tab_files[[meta1]] = gr_tab_samples[match(gr_tab_files$repl, 
                                               gr_tab_samples$sample),][[meta1]] 
}

gr_tab_files = gr_tab_files[!is.na(gr_tab_files$sample),]

t1 = gr_tab_files
t2 = t1[order(match(t1$sample, gr_tab_samples$sample)),]
gr_tab_files = t2



#save group table (keep copy of old version, if existing)
if (file.exists(paste0(out_dir,"/image_groups.csv"))){
  write_csv(gr_tab_files, file = paste0(out_dir,"/image_groups_new.csv"))
} else {write_csv(gr_tab_files, file = paste0(out_dir,"/image_groups.csv"))}


