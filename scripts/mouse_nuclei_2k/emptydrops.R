
setwd("../../")

library(diem)
source("scripts/common/empty_drops_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "mouse_nuclei_2k"
method <- "emptydrops"
dir10X <- "data/raw/mouse_nuclei_2k/raw_gene_bc_matrices/mm10/"
lab_ids <- "mouse_brain"

#=========================================
#=========================================

# Make directories

counts <- diem::read_10x(dir10X)
ed_pipe_out <- empty_drops_pipe(counts, dir_label=label, project=label)

# Run Seurat
seurat_pipe_single(ed_pipe_out$counts, dir_label=label, project=label, method=method)
