#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=16G,h_rt=6:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -t 1-11
#$ -o diem.k.R.log.$TASK_ID

setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem.k"
dir10X <- c("data/raw/AT1/",
            "data/raw/AT2/",
            "data/raw/AT3/",
            "data/raw/AT4/",
            "data/raw/AT5/",
            "data/raw/AT6/")

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

infn <- "data/raw/splice_frctn/midpoints.txt"
midpoints <- read.table(infn, header = FALSE, row.names = 1)
midpoints[,1] <- 100 * midpoints[,1]

args = commandArgs(TRUE)
j = as.integer(args[1])

j = Sys.getenv(x = "SGE_TASK_ID")
if (j == ""){
    stop("Set SGE_TASK_ID environment variable")
}
j <- as.integer(j)

k_vals <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

k <- k_vals[j]
dir_label <- "atsn"
methd <- paste0("diem.k/k", as.character(k))

for (i in 1:6){
    counts <- diem::read_10x(dir10X[i])

    sce <- diem_pipe1(counts, 
                      dir_label = dir_label, 
                      project = lab_ids[i], 
                      method = methd, 
                      k_init = k, 
                      model = "mltn", 
                      threads = 8) 

    sce <- diem_pipe2(sce, 
                      dir_label = label, 
                      project = lab_ids[i], 
                      method = methd, 
                      thresh = 0.5)
    dp <- paste0("data/processed/", label, "/", methd, "/")
    scefn <- paste0(dp, lab_ids[i], ".diem_sce.rds")
    #sce <- readRDS(scefn)

    ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
    sf <- read.table(ifn, header = TRUE, row.names = 1)
    keep <- grep(lab_ids[i], rownames(sf))
    sf <- sf[keep,,drop=FALSE]
    rownames(sf) <- sapply(rownames(sf), function(s) {
                           s <- strsplit(s, "_")
                           s <- s[[1]][length(s[[1]])]
                           return(s) })
    ts <- rownames(sce@test_data)
    sce@droplet_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
    sce@test_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]

    sce@test_data[is.na(sce@test_data[ts,"SpliceFrctn"]), "SpliceFrctn"] <- 110

    sce@test_data[ts,"Truth"] <- rep("Nuclear", length(ts))
    keep <- sce@test_data[ts,"SpliceFrctn"] >= midpoints[lab_ids[i],1]
    sce@test_data[keep,"Truth"] <- "Background"
    sce@test_data[sce@test_data[ts,"SpliceFrctn"] == 110,"Truth"] <- NA
    sce@test_data[sce@test_data[ts,"SpliceFrctn"] == 110,"SpliceFrctn"] <- NA
    saveRDS(sce, scefn)
}

