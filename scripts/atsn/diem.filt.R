
setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

library(future)
plan("multiprocess", workers = 4)

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem"
dir10X <- c("data/raw/AT1/", 
            "data/raw/AT2/",
            "data/raw/AT3/",
            "data/raw/AT4/",
            "data/raw/AT5/", 
            "data/raw/AT6/")

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

args = commandArgs(TRUE)
i = as.integer(args[1])

for (j in 1:6){
    counts <- diem::read_10x(dir10X[j])

    dp <- paste0("data/processed/", label, "/", method, "/")
    dir.create(dp, recursive=TRUE, showWarnings=FALSE)
    d_out <- paste0(dp, "filt/")
    dir.create(d_out, recursive=TRUE, showWarnings=FALSE)

    k_init <- 30
    filt_vals <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

    f <- filt_vals[i]

    n <- paste0(lab_ids, ".", as.character(f))
    sce <- create_SCE(counts, name=n)
    sce <- diem(sce, k_init = k_init, fltr = f)
    meta.data <- sce@droplet_data

    write.table(meta.data, 
                paste0(d_out, "meta.data.", lab_ids[j], ".", as.character(f), ".txt"), 
                row.names = TRUE, 
                col.names = NA, 
                quote = FALSE, 
                sep = "\t")
}

