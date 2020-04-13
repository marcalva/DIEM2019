#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -l h_data=4G,h_vmem=4G,h_rt=8:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -o diem.2.R.log


setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "mouse_nuclei_2k"
method <- "diem"
dir10X <- "data/raw/mouse_nuclei_2k/raw_gene_bc_matrices/mm10/"
lab_ids <- "mouse_brain"

#=========================================
#=========================================

diem_pipe_out <- diem_pipe2(dir_label = label, 
                            project = label, 
                            thresh = 0.5)

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
keep <- grep("mouse-nuclei_2k", rownames(sf))
sf <- sf[keep,,drop=FALSE]
rownames(sf) <- sapply(rownames(sf), function(s) {
                       s <- strsplit(s, "_")
                       s <- s[[1]][length(s[[1]])]
                       return(s) })

ts <- diem_pipe_out@test_set
diem_pipe_out@droplet_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
diem_pipe_out@test_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
dp <- paste0("data/processed/", label, "/", method, "/")
saveRDS(diem_pipe_out, paste0(dp, label, ".diem_sce.rds"))

# Run Seurat
filtered <- get_clean_ids(diem_pipe_out)
counts <- diem_pipe_out@counts[,filtered]
md <- diem_pipe_out@test_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)

