
# Run analysis using inflection point (by DropletUtils) for filtering

setwd("../../")

library(diem)
source("scripts/common/standard_seurat.R")
source("scripts/common/quantile_pipe.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "quantile"
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

if (is.na(i)){
	stop("Specify integer ID in Rscript command.")
}

counts <- diem::read_10x(dir10X[i])
quantile_pipe_out <- quantile_pipe(counts, dir_label=label, project=lab_ids[i])
