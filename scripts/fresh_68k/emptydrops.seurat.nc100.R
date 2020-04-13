#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=2:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -o emptydrops.seurat.nc100.R.log

setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/standard_seurat.R")
library(future)
plan("multiprocess", workers = 8)


#=========================================
# Set variables
#=========================================

label <- "fresh_68k"
method <- "emptydrops.nc100"
dir10X <- "data/raw/fresh_68k/matrices_mex/hg19/"
lab_ids <- "fresh_68k"

#=========================================
#=========================================

# Make directories

dp <- paste0("data/processed/", label, "/", method, "/")

ed <- readRDS("data/processed/fresh_68k/emptydrops/fresh_68k.emptydrops_counts.rds")

e.out <- ed$e.out[ ed$e.out$Total >= 100,]
keep <- rownames(e.out)[!is.na(e.out$FDR)]
e.out <- e.out[keep,]
e.out.pass <- e.out[e.out$FDR < 0.05,]
ed_counts <- ed$counts[,rownames(e.out.pass)]


seurat_pipe_single(ed_counts, 
                   min.features = 0, 
                   dir_label=label, 
                   project=label, 
                   method=method)

