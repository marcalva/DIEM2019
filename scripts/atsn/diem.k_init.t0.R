
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
    d_out <- paste0(dp, "k_init.t0/")
    dir.create(d_out, recursive=TRUE, showWarnings=FALSE)

    flt <- 0
    k_vals <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    k_init <- k_vals[i]

    n <- paste0(lab_ids, ".", as.character(k_init))
    sce <- create_SCE(counts, name=n)
    sce <- diem(sce, k_init = k_init, fltr = flt)
    meta.data <- sce@droplet_data

    write.table(meta.data, 
                paste0(d_out, "meta.data.", lab_ids[j], ".", as.character(k_init), ".txt"), 
                row.names = TRUE, 
                col.names = NA, 
                quote = FALSE, 
                sep = "\t")
}

