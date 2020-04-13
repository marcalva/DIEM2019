#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 2
#$ -l h_data=8G,h_vmem=8G,h_rt=1:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -t 1-6
#$ -o diem.2.R.log.$TASK_ID


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

#=========================================
#=========================================

diem_pipe_out <- diem_pipe2(dir_label = label, 
                            project = lab_ids[i], 
                            thresh = 0.5)

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

ts <- diem_pipe_out@test_set
diem_pipe_out@droplet_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
diem_pipe_out@test_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
dp <- paste0("data/processed/", label, "/", method, "/")
saveRDS(diem_pipe_out, paste0(dp, lab_ids[i], ".diem_sce.rds"))

# Run Seurat
message("Running Seurat")

filtered <- get_clean_ids(diem_pipe_out)
counts <- diem_pipe_out@counts[,filtered]
md <- diem_pipe_out@test_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)

