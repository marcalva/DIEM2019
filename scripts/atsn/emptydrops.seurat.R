
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

seurat_pipe_list(countsl, dir_label=label, project.l=lab_ids, project=label, method=method, scale.factor=1e3)

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))
seur$SpliceFrctn <- 100 * sf[colnames(seur),1]
saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))


# Rename cells
tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.0.8, mean)

from <- seq(0,13)
to <- as.character(from)
to <- c("Adipocyte-1", 
        "Stromal-1", 
        "Macrophage", 
        "Adipocyte-2", 
        "Adipocyte-3", 
        "Endothelial",  # 5
        "Debris-1", 
        "Debris-2", 
        "T-cell", 
        "Perivascular", 
        "Stromal-2", #10
        "Dendritic", 
        "Debris-3", 
        "Mast")

map <- to
names(map) <- from

seur@meta.data$CellType <- map[as.character(seur@meta.data$RNA_snn_res.0.8)]

saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

