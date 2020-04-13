
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

quantile.l <- lapply(lab_ids, function(x) {
				 readRDS(paste0(dp, x, ".quantile_counts.rds"))})

countsl <- lapply(quantile.l, function(x) {
	return(Matrix::Matrix(x$counts, sparse=TRUE))
	})

seurat_pipe_list(countsl, dir_label=label, project.l=lab_ids, project=label, method=method, scale.factor=1e3)


# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))
seur$SpliceFrctn <- 100 * sf[colnames(seur),1]
saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

seur$SpliceFrctn <- 100 * sf[colnames(seur),1]

saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

tapply(seur@meta.data$SpliceFrctn, seur@meta.data$RNA_snn_res.0.8, mean)

from <- seq(0,14)
to <- as.character(from)
to <- c("Adipocyte-1", 
        "Adipocyte-2", 
        "Stromal-1", 
        "Adipocyte-3", 
        "Macrophage", 
        "Adipocyte-4",  # 5
        "Debris-1", 
        "Endothelial-1", 
        "Debris-2", 
        "Perivascular", 
        "Dendritic", #10
        "T-cell", 
        "Endothelial-2", 
        "Mast", 
        "Adipocyte-5")

map <- to
names(map) <- from

seur@meta.data$CellType <- map[as.character(seur@meta.data$RNA_snn_res.0.8)]

saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))

