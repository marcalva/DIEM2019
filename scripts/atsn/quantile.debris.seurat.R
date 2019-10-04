
setwd("../../")

source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

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
names(dir10X) <- lab_ids

#=========================================
#=========================================

# Make directories

	# Create directories
dp <- paste0("data/processed/", label, "/", method, "/")

quantile.l <- lapply(lab_ids, function(x) {
				 readRDS(paste0(dp, x, ".quantile_counts.rds"))})
names(quantile.l) <- lab_ids

countsl <- lapply(lab_ids, function(x) {
                  print(x)
                  counts <- diem::read_10x(dir10X[x])
                  ng <- Matrix::colSums(counts > 0)
                  nc <- Matrix::colSums(counts)
                  keep <- colnames(counts)[ ng >= 200 & nc < quantile.l[[x]]$thresh ]

        return(counts[,keep])
    })

method <- "quantile_debris"
dp <- paste0("data/processed/", label, "/", method, "/")

seurat_pipe_list(countsl, dir_label=label, project.l=lab_ids, project_out=label, method=method, scale.factor=1e3)
