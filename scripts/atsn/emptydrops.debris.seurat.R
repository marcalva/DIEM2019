
setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/standard_seurat.R")
source("scripts/common/plotting.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "emptydrops"
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

emptydrops.l <- lapply(lab_ids, function(x) {
                       read.table(paste0("results/atsn/emptydrops/", x, ".emptydrops.out.txt", sep=""), stringsAsFactors=FALSE)
})
names(emptydrops.l) <- lab_ids

countsl <- lapply(lab_ids, function(x) {
                  print(x)
                  counts <- diem::read_10x(dir10X[x])
                  edo <- emptydrops.l[[x]]
                  edo <- edo[!is.na(edo$FDR),]
                  re <- rownames(edo)[edo$FDR > .05]
	    return(counts[,re])
	})

# slot 2 only has 1 sample

method <- "emptydrops_debris"
dp <- paste0("data/processed/", label, "/", method, "/")

seurat_pipe_list(countsl, dir_label=label, project.l=lab_ids, project=label, method=method, scale.factor=1e3)

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

seur <- readRDS(paste0(dp, label, ".seur_obj.rds"))

seur$SpliceFrctn <- 100 * sf[colnames(seur),1]

saveRDS(seur, paste0(dp, label, ".seur_obj.rds"))


