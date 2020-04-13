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
#$ -t 1-6
#$ -o diem.R.log.$TASK_ID


setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem"
dir10X <- c("data/raw/AT1/", 
            "data/raw/AT2/",
            "data/raw/AT3/",
            "data/raw/AT4/",
            "data/raw/AT5/", 
            "data/raw/AT6/")

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

args = commandArgs(TRUE)
i = as.integer(args[1])

i = Sys.getenv(x = "SGE_TASK_ID")
if (i == ""){
    stop("Set SGE_TASK_ID environment variable")
}
i <- as.integer(i)
counts <- diem::read_10x(dir10X[i])

#=========================================
#=========================================

dp <- paste0("data/processed/", label, "/", method, "/")
sce_fn <- paste0(dp, lab_ids[i], ".diem_sce.rds")

k_init <- 20

sce <- diem_pipe1(counts, 
                  dir_label = label, 
                  project = lab_ids[i], 
                  k_init = k_init, 
                  model = "mltn", 
                  threads = 8) 

sce <- diem_pipe2(sce, 
                  dir_label = label, 
                  project = lab_ids[i])

# Add percent of reads spliced
message("Adding percent spliced")
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
keep <- grep(lab_ids[i], rownames(sf))
sf <- sf[keep,,drop=FALSE]
rownames(sf) <- sapply(rownames(sf), function(s) {
                       s <- strsplit(s, "_")
                       s <- s[[1]][length(s[[1]])]
                       return(s) })

ts <- sce@test_set
sce@droplet_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
sce@test_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
saveRDS(sce, sce_fn)

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

