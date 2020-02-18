
# Get sample IDs

setwd("../../")

library(Matrix)
library(diem)

# indir <- "data/raw/mouse_nuclei_2k/raw_gene_bc_matrices/mm10/"
indir <- "data/raw/mouse_nuclei_2k/raw_feature_bc_matrix/"
counts <- read_10x(indir)
sizes <- colSums(counts)
o <- order(sizes, decreasing=TRUE)
droplet_ids <- colnames(counts)[o][1:1e4]
odir <- "data/processed/mouse_nuclei_2k/velocyto/"
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
ofile <- paste0(odir, "top10K_ids.txt")
write.table(droplet_ids, ofile, row.names = FALSE, col.names = FALSE, quote = FALSE)

