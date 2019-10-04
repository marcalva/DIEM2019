
setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "emptydrops"
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

# Make directories

	# Create directories
dp <- paste0("data/processed/", label, "/", method, "/")

emptydrops.l <- lapply(lab_ids, function(x) {
				 readRDS(paste0(dp, x, ".emptydrops_counts.rds"))})

countsl <- lapply(emptydrops.l, function(x) {
	return(x$counts)
	})

# slot 2 only has 1 sample
sizes <- sapply(countsl, function(x) ncol(x))

keep <- sizes >= 30

countsl <- countsl[keep]
lab_ids <- lab_ids[keep]

seurat_pipe_list(countsl, dir_label=label, project.l=lab_ids, project_out=label, method=method, scale.factor=1e3)

