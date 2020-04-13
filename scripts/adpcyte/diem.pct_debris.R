
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
sce <- readRDS(diem_pipe_out, paste0(dp, label, ".diem_sce.rds"))

sce <- estim_prop_nnls(sce)
sce <- estimate_db(sce)
sce <- correct_counts(sce)


dir_res <- paste0("results/", label, "/", method, "/pct_debris/")
dir_plot <- paste0(dir_res, "/plots/")

saveRDS(sce, paste0(dp, label, ".diem_sce.rds"))

# Add percent of reads spliced
ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
keep <- grep("adpcyte", rownames(sf))
sf <- sf[keep,,drop=FALSE]
rownames(sf) <- sapply(rownames(sf), function(s) {
                       s <- strsplit(s, "_")
                       s <- s[[1]][length(s[[1]])]
                       return(s) })

ts <- diem_pipe_out@test_set
diem_pipe_out@droplet_data[ts,"SpliceFrctn"] <- 100*sf[ts,1]
dp <- paste0("data/processed/", label, "/", method, "/")
saveRDS(diem_pipe_out, paste0(dp, label, ".diem_sce.rds"))

# Run Seurat
filtered <- get_clean_ids(diem_pipe_out)
counts <- diem_pipe_out@counts[,filtered]
md <- diem_pipe_out@droplet_data[filtered,]

seurat_pipe_single(counts, 
                   dir_label = label, 
                   project = label, 
                   meta.data = md, 
                   method = method)

