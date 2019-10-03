
setwd("../../")

library(diem)
source("scripts/common/empty_drops_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "adpcyte"
method <- "emptydrops"
dir10X <- "data/raw/adpcyte/"
lab_ids <- "ADCY"

#=========================================
#=========================================

# Make directories

counts <- diem::read_10x(dir10X)
# ed_pipe_out <- empty_drops_pipe(counts, dir_label=label, project=label, retain=Inf, lower=150)
ed_pipe_out <- empty_drops_pipe(counts, dir_label=label, project=label)

# Run Seurat
seurat_pipe_single(ed_pipe_out$counts, dir_label=label, project=label, method=method)
