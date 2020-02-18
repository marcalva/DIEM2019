
setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem_integrated.fltr.1"
dir10X <- c("data/raw/AT1/", 
            "data/raw/AT2/",
            "data/raw/AT3/",
            "data/raw/AT4/",
            "data/raw/AT5/", 
            "data/raw/AT6/")

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

counts <- lapply(dir10X, diem::read_10x)
counts <- do.call(cbind, counts)

diem_pipe_out <- diem_pipe(counts, 
                           dir_label = label, 
                           project = "atsn", 
                           method = method, 
                           k_init = 30, 
                           fltr = 0.1, 
                           top_n = NULL)

# Run Seurat
filtered <- get_clean_ids(diem_pipe_out)
counts <- diem_pipe_out@counts[,filtered]
md <- diem_pipe_out@droplet_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)

