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
#$ -o emptydrops.R.log
setwd("../../")

library(diem)
source("scripts/common/empty_drops_pipe.R")
source("scripts/common/standard_seurat.R")

library(future)
plan("multiprocess", workers = 8)

#=========================================
# Set variables
#=========================================

label <- "fresh_68k"
method <- "emptydrops"
dir10X <- "data/raw/fresh_68k/matrices_mex/hg19/"
lab_ids <- "fresh_68k"

#=========================================
#=========================================

# Make directories

counts <- diem::read_10x(dir10X, clip_end=FALSE)
# ed_pipe_out <- empty_drops_pipe(counts, dir_label=label, project=label, retain=Inf, lower=150)
ed_pipe_out <- empty_drops_pipe(counts, dir_label=label, project=label)

# Run Seurat
seurat_pipe_single(ed_pipe_out$counts, dir_label=label, project=label, method=method)

