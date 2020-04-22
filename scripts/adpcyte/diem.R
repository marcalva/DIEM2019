#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=2:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -o diem.R.log


setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "adpcyte"
method <- "diem"
dir10X <- "data/raw/adpcyte/"
lab_ids <- "ADCY"

#=========================================
#=========================================

dp <- paste0("data/processed/", label, "/", method, "/")
fn <- paste0(dp, label, ".diem_sce.rds")

counts <- diem::read_10x(dir10X)
sce <- diem_pipe1(counts, 
                  dir_label = label, 
                  project = label, 
                  k_init = 20, 
                  threads = 8) 

sce <- diem_pipe2(sce, 
                  dir_label = label, 
                  project = label)

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
keep <- grep("adpcyte", rownames(sf))
sf <- sf[keep,,drop=FALSE]
rownames(sf) <- sapply(rownames(sf), function(s) {
                       s <- strsplit(s, "_")
                       s <- s[[1]][length(s[[1]])]
                       return(s) })

ts <- sce@test_set
sce@droplet_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
sce@test_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
dp <- paste0("data/processed/", label, "/", method, "/")
saveRDS(sce, paste0(dp, label, ".diem_sce.rds"))

# Run Seurat
filtered <- get_clean_ids(sce)
counts <- sce@counts[,filtered]
md <- sce@test_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)

