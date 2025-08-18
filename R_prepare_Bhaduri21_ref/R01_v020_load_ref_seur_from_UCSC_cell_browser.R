getwd()
setwd("/rds/general/user/mlattke/projects/dsseq23/live/E02_230608_foetal_brain_ref_atlas/")

message("####### \n ####### Start R01: ", Sys.time(), "\n #######")

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(data.table)
library(BPCells)

set.seed(1234)

# define output folder
out_dir = "./R_out_Seur5/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}




#####################################
# get count data for reference datasets
#####################################

#meta <- data.frame(fread("https://cells.ucsc.edu/dev-brain-regions/wholebrain/meta.tsv"), row.names=1)
# or after manual download 
meta <- data.frame(fread("./R_in/meta.tsv"), row.names=1)

meta[1:10,]
message("    metadata loaded ", Sys.time())

#mat <- fread("https://cells.ucsc.edu/dev-brain-regions/wholebrain/exprMatrix.tsv.gz")
# or after manual download 
mat <- fread("./R_in/exprMatrix.tsv.gz")

message("    mat loaded ", Sys.time())

genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
mat = as(mat, "sparseMatrix") 
write_matrix_dir(mat = mat, dir = './R_in/counts_Bhaduri21')

message("    start creating seurat object ", Sys.time())

# merge count matrices and metadata 
seur_ref = CreateSeuratObject(counts = open_matrix_dir(dir = './R_in/counts_Bhaduri21'), 
                        meta.data = meta)

save(seur_ref, file = paste0(out_dir,"R01_seur_ref.rda")) 


message("####### End R01: ", Sys.time())



