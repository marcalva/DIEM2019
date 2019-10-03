
setwd("../../")

source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "quantile"
snames <- c("110686_121444", "54172_58356", "57289_64605", "62630_75306", "63099_76245", "63099_76246")
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

# Make directories

	# Create directories
dp <- paste0("data/processed/", label, "/", method, "/")

quantile.l <- lapply(lab_ids, function(x) {
				 readRDS(paste0(dp, x, ".quantile_counts.rds"))})

countsl <- lapply(quantile.l, function(x) {
	return(Matrix::Matrix(x$counts, sparse=TRUE))
	})

seurat_pipe_list(countsl, dir_label=label, project.l=lab_ids, project_out=label, method=method, scale.factor=1e3)
