
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

seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

# Clustering resolution 2
seur <- FindClusters(seur, resolution = 2)
#markers <- FindAllMarkers(seur, only.pos=TRUE, resolution = 2)
#write.table(markers, "results/atsn/emptydrops/atsn.seur_markers.res2.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
from <- seq(0,22)
to <- c("Myl-1", 
        "Adp-1", 
        "Stm-1", 
        "Adp-2", 
        "Vsc-1", 
        "Adp-3", # 5
        "Adp-4", 
        "Stm-2", 
        "T", 
        "Dbr-1", 
        "Dbr-2", #10
        "Myl-2", 
        "Vsc-2", 
        "Adp-5", 
        "Stm-3", 
        "Dbr-3", #15
        "Adp-6", 
        "Dblt-1", 
        "Adp-7", 
        "Mast", 
        "Vsc-3", #20
        "Adp-8", 
        "Adp-9")
map <- to
names(map) <- from
seur@meta.data$CellType <- as.character(seur@meta.data$RNA_snn_res.2)
seur@meta.data$CellType <- map[seur@meta.data$CellType]
saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

markers <- read.table("results/atsn/emptydrops/atsn.seur_markers.res2.txt", 
                      stringsAsFactors=FALSE, header = TRUE, sep = "\t") 
markers[,"cluster"] <- map[as.character(markers[,"cluster"])]
write.table(markers, "results/atsn/emptydrops/atsn.seur_markers.res2.ct.txt", 
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.2, mean)

