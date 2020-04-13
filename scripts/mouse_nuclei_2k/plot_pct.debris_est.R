
setwd("../../")

library(diem)
library(ggplot2)
library(reshape2)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "mouse_nuclei_2k"
method <- "diem"
dir10X <- "data/raw/mouse_nuclei_2k/raw_gene_bc_matrices/mm10/"
lab_ids <- "mouse_brain"

ppre <- label
ptitle <- "Mouse Brain"

#=========================================
#=========================================

dp <- paste0("data/processed/", label, "/", method, "/")
dir_plot <- paste0("results/", label, "/diem/plots/")

sce <- readRDS(paste0(dp, label, ".diem_sce.rds"))

infn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(infn, header = TRUE, row.names = 1)
keep <- grep("mouse-nuclei_2k", rownames(sf))
sf <- sf[keep,,drop=FALSE]
rownames(sf) <- sapply(rownames(sf), function(s) {
                       s <- strsplit(s, "_")
                       s <- s[[1]][length(s[[1]])]
                       return(s) })

dropdat <- droplet_data(sce)
drops <- intersect(rownames(sf), rownames(dropdat))

#=========================================
# Plot percent debris against variables
#=========================================

datf <- data.frame("PercentDebris" = dropdat[drops,"pct.debris"], 
                   "SpliceFraction" = sf[drops,"SpliceFrctn"],
                   "PercentMT" = dropdat[drops,"pct.mt"], 
                   "LogUMI" = log10(dropdat[drops,"total_counts"]))

datfm <- melt(datf, id.vars = "PercentDebris")

p <- ggplot(datfm, aes(x = PercentDebris, y = value)) + 
geom_point(shape = 16, alpha = 0.2) + 
facet_wrap(~variable, scales = "free_y") + 
ggtitle(ptitle) + 
theme_minimal() + 
theme(axis.title.y = element_blank(), 
      text = element_text(size = 18), 
      plot.title = element_text(hjust = 0.5))

ggsave(paste0(dir_plot, ppre, ".pct_debris.vars.jpeg"), 
       height = 5, width = 10)

print(cor.test(datf[,"PercentDebris"], datf[,"SpliceFraction"]))
print(cor.test(datf[,"PercentDebris"], datf[,"PercentMT"]))
print(cor.test(datf[,"PercentDebris"], datf[,"LogUMI"]))

#=========================================
# Plot percent debris for clean vs. debris
#=========================================

datf <- data.frame("PercentDebris" = dropdat[drops,"pct.debris"], 
                   "Call" = factor(dropdat[drops,"Call"]))

p <- ggplot(datf, aes(x = Call, y = PercentDebris)) + 
geom_jitter(shape = 16, alpha = 0.2) + 
ggtitle(ptitle) + 
theme_minimal() + 
theme(text = element_text(size = 18), 
     plot.title = element_text(hjust = 0.5))

ggsave(paste0(dir_plot, ppre, ".pct_debris.call.jpeg"), 
       height = 5, width = 4)

t.test(datf[datf$Call == "Clean","PercentDebris"], 
       datf[datf$Call == "Debris","PercentDebris"])

