
setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem"
snames <- c("110686_121444", "54172_58356", "57289_64605", "62630_75306", "63099_76245", "63099_76246")
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

# Make directories
dp <- paste0("data/processed/", label, "/", method, "/"); dir.create(dp, recursive=TRUE, showWarnings=FALSE)
dr <- paste0("results/", label, "/", method, "/"); dir.create(dr, recursive=TRUE, showWarnings=FALSE)
dir_plot <- paste0(dr, "plots/"); dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

#=========================================
#=========================================

scel <- lapply(lab_ids, function(x) {
			   readRDS(paste0(dp, x, ".diem_sce.rds"))})

countsl <- lapply(scel, function(x) {
	filtered <- get_clean_ids(x)
	counts <- x@counts[,filtered]
	return(counts)
	})
mdl <- lapply(scel, function(x) {
	filtered <- get_clean_ids(x)
	md <- x@droplet_data[filtered,]
	return(md)
	})

sizes <- sapply(countsl, function(x) ncol(x))

keep <- sizes >= 30

countsl <- countsl[keep]
lab_ids <- lab_ids[keep]

seurat_pipe_list(countsl, 
                 dir_label=label, 
                 project.l=lab_ids, 
                 meta.data.l = mdl, 
                 project_out=label, 
                 method=method, 
                 scale.factor=1e3)

