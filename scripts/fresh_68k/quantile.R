
# Run analysis using inflection point (by DropletUtils) for filtering

setwd("../../")

source("scripts/common/standard_seurat.R")
source("scripts/common/quantile_pipe.R")

#=========================================
# Set variables
#=========================================


label <- "fresh_68k"
method <- "quantile"
dir10X <- "data/raw/fresh_68k/matrices_mex/hg19/"
lab_ids <- "fresh_68k"

#=========================================
#=========================================

counts <- diem::read_10x(dir10X, clip_end=FALSE)
quantile_pipe_out <- quantile_pipe(counts, dir_label=label, project=label)

# Seurat
seurat_pipe_single(quantile_pipe_out$counts, dir_label=label, project=label, method=method)
