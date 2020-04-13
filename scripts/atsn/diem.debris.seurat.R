
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
                  debris <- get_removed_ids(x)
                  counts <- x@counts[,debris]
                  return(counts)
	})
mdl <- lapply(scel, function(x) {
              debris <- get_removed_ids(x)
              md <- x@droplet_data[debris,]
              return(md)
    })

method <- "diem_debris"
dp <- paste0("data/processed/", label, "/", method, "/"); dir.create(dp, recursive=TRUE, showWarnings=FALSE)
dr <- paste0("results/", label, "/", method, "/"); dir.create(dr, recursive=TRUE, showWarnings=FALSE)
dir_plot <- paste0(dr, "plots/"); dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

seurat_pipe_list(countsl, 
                 dir_label = label, 
                 project.l = lab_ids, 
                 meta.data.l = mdl, 
                 project = label, 
                 method = method, 
                 scale.factor = 1e3)

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

seur$SpliceFrctn <- 100 * sf[colnames(seur),1]

saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.0.8, mean)

# Find markers
markers <- FindAllMarkers(seur, only.pos = TRUE)

from <- seq(0,7)
to <- c("Dbr-1", 
        "Dbr-2", 
        "Dbr-3", 
        "Dbr-4", 
        "Dbr-5", 
        "Dbr-6", 
        "Dbr-7", 
        "Dbr-8") 

map <- to
names(map) <- from

seur@meta.data$CellType <- map[as.character(seur@meta.data$RNA_snn_res.0.8)]

saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

