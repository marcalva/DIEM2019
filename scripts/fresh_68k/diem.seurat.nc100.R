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
#$ -o diem.seurat.nc100.R.log

setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/standard_seurat.R")
library(future)
plan("multiprocess", workers = 8)

options(future.globals.maxSize = 4092 * 1024^2)

#=========================================
# Set variables
#=========================================

label <- "fresh_68k"
method <- "diem.nc100"
dir10X <- "data/raw/fresh_68k/matrices_mex/hg19/"
lab_ids <- "fresh_68k"

#=========================================
#=========================================

dp <- paste0("data/processed/", label, "/", method, "/")

ed <- readRDS("data/processed/fresh_68k/emptydrops/fresh_68k.emptydrops_counts.rds")
sce <- readRDS("data/processed/fresh_68k/diem/fresh_68k.diem_sce.rds")

e.out <- ed$e.out[ ed$e.out$Total >= 100,]
testdat <- sce@test_data

keep <- rownames(e.out)[!is.na(e.out$FDR)]
e.out <- e.out[keep,]
sce@test_data <- sce@test_data[keep,]
testdat <- testdat[keep,]
testdat.pass <- testdat[testdat$Cluster != 1,]
diem_counts <- sce@counts[,rownames(testdat.pass)]

sce <- call_targets(sce, thresh = NULL, min_genes = 0, clusters = 1)

# Run Seurat
filtered <- get_clean_ids(sce)
counts <- sce@counts[,filtered]
md <- sce@test_data[filtered,]

seurat_pipe_single(counts, 
                   min.features = 0, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)

