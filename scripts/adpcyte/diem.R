
setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "adpcyte"
method <- "diem"
dir10X <- "data/raw/adpcyte/"
lab_ids <- "ADCY"

#=========================================
#=========================================

counts <- diem::read_10x(dir10X)
diem_pipe_out <- diem_pipe(counts, dir_label=label, project=label)

# Run Seurat
filtered <- get_clean_ids(diem_pipe_out)
counts <- diem_pipe_out@counts[,filtered]
md <- diem_pipe_out@droplet_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label=label, 
                   project=label, 
                   meta.data = md, 
                   method=method)
