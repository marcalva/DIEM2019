#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 16
#$ -l h_data=4G,h_vmem=64G,h_rt=10:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -o diem.integrated.R.log

setwd("../../")

library(diem)
library(Seurat)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem_integrated"
dir10X <- c("data/raw/AT1/", 
            "data/raw/AT2/",
            "data/raw/AT3/",
            "data/raw/AT4/",
            "data/raw/AT5/", 
            "data/raw/AT6/")

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

counts <- lapply(1:6, function(i){
                 dt <- dir10X[i]
                 lb <- lab_ids[i]
                 ct <- diem::read_10x(dt)
                 colnames(ct) <- paste0(lb, "_", colnames(ct))
                 return(ct)
            })

counts <- do.call(cbind, counts)

k_init <- 20

sce <- diem_pipe1(counts, 
                  dir_label = label, 
                  project = label, 
                  method = method, 
                  k_init = k_init, 
                  threads = 16) 

sce <- diem_pipe2(sce, 
                  dir_label = label, 
                  project = label, 
                  method = method)

# Add percent of reads spliced
message("Adding percent spliced")
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

for (i in 1:6){
    keep <- grep(lab_ids[i], rownames(sf))
    sf <- sf[keep,,drop=FALSE]
    rownames(sf) <- sapply(rownames(sf), function(s) {
                           s <- strsplit(s, "_")
                           s <- s[[1]][length(s[[1]])]
                           return(s) })

    ts <- sce@test_set
    drops <- intersect(ts, rownames(sf))
    sce@droplet_data[drops,"SpliceFrctn"] <- 100*sf[drops,1]
    sce@test_data[drops,"SpliceFrctn"] <- 100*sf[drops,1]
}

dp <- paste0("data/processed/", label, "/", method, "/")
scefn <- paste0(dp, label, ".diem_sce.rds")
saveRDS(sce, scefn)

# Run Seurat
message("Running Seurat")

filtered <- get_clean_ids(sce)
counts <- sce@counts[,filtered]
md <- sce@test_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)

