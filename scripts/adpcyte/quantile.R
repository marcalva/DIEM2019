
# Run analysis using inflection point (by DropletUtils) for filtering

setwd("../../")

source("scripts/common/standard_seurat.R")
source("scripts/common/quantile_pipe.R")

#=========================================
# Set variables
#=========================================

label <- "adpcyte"
method <- "quantile"
dir10X <- "data/raw/adpcyte/"
lab_ids <- "ADCY"

#=========================================
#=========================================

counts <- diem::read_10x(dir10X)
quantile_pipe_out <- quantile_pipe(counts, dir_label=label, project=label)

# Seurat
seurat_pipe_single(quantile_pipe_out$counts, dir_label=label, project=label, method=method)
