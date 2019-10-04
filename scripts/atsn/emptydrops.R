
setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/empty_drops_pipe.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "emptydrops"
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

counts <- diem::read_10x(dir10X[i])
ed_pipe_out <- empty_drops_pipe(counts, dir_label=label, project=lab_ids[i])

