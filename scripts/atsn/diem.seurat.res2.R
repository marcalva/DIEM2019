
setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/standard_seurat.R")

library(future)
plan("multiprocess", workers = 8)

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

seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.0.8, mean)

# Clustering resolution 2
seur <- FindClusters(seur, resolution = 2)
markers <- FindAllMarkers(seur, only.pos=TRUE, resolution = 2)
write.table(markers, "results/atsn/diem/atsn.seur_markers.res2.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
from <- seq(0,20)
to <- c("Stm-1", 
        "Myl-1", 
        "Adp-1", 
        "Adp-2", 
        "Adp-3", 
        "Adp-4",  # 5
        "Adp-5", 
        "Vsc-1", 
        "Stm-2", 
        "Adp-6", 
        "Adp-7", #10
        "T", 
        "Stm-3", 
        "Vsc-2", 
        "Myl-2", 
        "Dblt-1", #15
        "Vsc-3", 
        "Mast", 
        "Myl-3", 
        "Dblt-2", 
        "Vsc-4") #20
map <- to
names(map) <- from
seur@meta.data$CellType <- as.character(seur@meta.data$RNA_snn_res.2)
seur@meta.data$CellType <- map[seur@meta.data$CellType]
saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

markers <- read.table("results/atsn/diem/atsn.seur_markers.res2.txt", 
                      stringsAsFactors=FALSE, header = TRUE, sep = "\t") 
markers[,"cluster"] <- map[as.character(markers[,"cluster"])]
write.table(markers, "results/atsn/diem/atsn.seur_markers.res2.ct.txt", 
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.2, mean)

