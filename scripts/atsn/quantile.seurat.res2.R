
setwd("../../")

source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "quantile"
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

# Make directories

dp <- paste0("data/processed/", label, "/", method, "/")

seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

# Clustering resolution 2
seur <- FindClusters(seur, resolution = 2)
#markers <- FindAllMarkers(seur, only.pos=TRUE, resolution = 2)
#write.table(markers, "results/atsn/quantile/atsn.seur_markers.res2.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
from <- seq(0,22)
to <- c("Adp-1", 
        "Stm-1", 
        "Adp-2", 
        "Adp-3", 
        "Adp-4", 
        "Adp-5", # 5
        "Adp-6", 
        "Dbr-1", 
        "Myl-1", 
        "Dbr-2", 
        "Dbr-3", #10
        "Vsc-1", 
        "Vsc-2", 
        "Stm-2", 
        "Myl-2", 
        "T", #15
        "Dbr-4", 
        "Adp-7", 
        "Adp-8", 
        "Adp-9", 
        "Dblt-1", #20
        "Vsc-3", 
        "Mast")
map <- to
names(map) <- from
seur@meta.data$CellType <- as.character(seur@meta.data$RNA_snn_res.2)
seur@meta.data$CellType <- map[seur@meta.data$CellType]
saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

markers <- read.table("results/atsn/quantile/atsn.seur_markers.res2.txt", 
                      stringsAsFactors=FALSE, header = TRUE, sep = "\t") 
markers[,"cluster"] <- map[as.character(markers[,"cluster"])]
write.table(markers, "results/atsn/quantile/atsn.seur_markers.res2.ct.txt", 
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.2, mean)
tapply(seur@meta.data$SpliceFrctn, seur@meta.data$CellType, mean)

