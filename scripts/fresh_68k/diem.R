#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 12
#$ -l h_data=4G,h_vmem=32G,h_rt=8:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -o diem.R.log

setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "fresh_68k"
method <- "diem"
dir10X <- "data/raw/fresh_68k/matrices_mex/hg19/"
lab_ids <- "fresh_68k"

#=========================================
#=========================================

dp <- paste0("data/processed/", label, "/", method, "/")
fn <- paste0(dp, label, ".diem_sce.rds")

counts <- diem::read_10x(dir10X, clip_end = FALSE)
sce <- diem_pipe1(counts, 
                  dir_label = label, 
                  project = label, 
                  k_init = 20, 
                  threads = 12) 

sce <- diem_pipe2(sce, 
                  dir_label = label, 
                  project = label)

# Call cells by cluster
sce <- call_targets(sce, thresh = NULL, clusters = 1)

# Run Seurat
filtered <- get_clean_ids(sce)
counts <- sce@counts[,filtered]
md <- sce@test_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)
