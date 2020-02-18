
# Run analysis using inflection point (by DropletUtils) for filtering

setwd("../../")

source("scripts/common/standard_seurat.R")
source("scripts/common/quantile_pipe.R")

#=========================================
# Set variables
#=========================================

label <- "mouse_nuclei_2k"
method <- "quantile"
dir10X <- "data/raw/mouse_nuclei_2k/raw_gene_bc_matrices/mm10/"
lab_ids <- "mouse_brain"

#=========================================
#=========================================

counts <- diem::read_10x(dir10X)
quantile_pipe_out <- quantile_pipe(counts, dir_label=label, project=label)

# Seurat
seurat_pipe_single(quantile_pipe_out$counts, dir_label=label, project=label, method=method)

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

dp <- paste0("data/processed/", label, "/", method, "/")
seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

seur$SpliceFrctn <- 100 * sf[colnames(seur),1]

saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

