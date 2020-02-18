
# Integrate DIEM samples using Seurat v3 instead of raw merge

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
	md <- x@droplet_data[filtered,]
	return(md)
	})

sizes <- sapply(countsl, function(x) ncol(x))

keep <- sizes >= 30

countsl <- countsl[keep]
lab_ids <- lab_ids[keep]

seur.list <- lapply(1:length(countsl), function(i) {
                    CreateSeuratObject(countsl[[i]], 
                                       project = lab_ids[[i]], 
                                       min.features=200, 
                                       meta.data=mdl[[i]])
    })

seur.list <- lapply(seur.list, seurat_norm, scale.factor=1e3)
anchors <- FindIntegrationAnchors(object.list = seur.list, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
integrated <- seurat_cluster(integrated)

saveRDS(integrated, paste0(dp, "atsn.seur_intg.seur_obj.rds"))



