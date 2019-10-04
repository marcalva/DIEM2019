
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

counts <- diem::read_10x(dir10X, clip_end=FALSE)
diem_pipe_out <- diem_pipe(counts, dir_label=label, project=label, cluster_n=30e3, top_n=80e4)

# Run Seurat
filtered <- get_clean_ids(diem_pipe_out)
counts <- diem_pipe_out@counts[,filtered]
md <- diem_pipe_out@droplet_data[filtered,]

seurat_pipe_single(counts, dir_label=label, project=label, method=method)