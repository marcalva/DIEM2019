
setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem"
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
	md <- x@test_data[filtered,]
    md[,"Sample"] <- x@name
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
                 project=label, 
                 method=method, 
                 scale.factor=1e3)

seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.0.8, mean)

