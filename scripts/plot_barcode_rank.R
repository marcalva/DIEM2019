
setwd("../")

library(diem)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
source("scripts/common/plotting.R")

#=========================================
# Plot functions
#=========================================

set_breaks_10 <- function(x){
    xmax <- x[2]
    bk <- 10
    brks <- c(bk)
    while (bk < xmax){
        bk <- bk * 10
        brks <- c(brks, bk)
    }
    return(brks)
}


ct <- theme(legend.position="none",
            text=element_text(size=16),
            plot.title=element_text(size=18, hjust=0.5, face="bold")
            )

yexp <- 1.1

xa_angle <- 30
hj <- 0.9
ts <- 18

#=========================================
#=========================================

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
expr <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- expr
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples

mthds <- c("DIEM", "EmptyDrops", "Quantile")

#=========================================
# Read in data
#=========================================

all_sf <- read.table("data/processed/meta_data.txt",
                     header=TRUE,
                     sep="\t")

all_sf <- all_sf[all_sf[,"n_genes"] >= 200,]

all_sf[,"Exprm"] <- factor(all_sf[,"Exprm"], levels=unique(expr))
all_sf[,"Sample"] <- factor(all_sf[,"Sample"], levels=samples)
all_sf[,"Label"] <- factor(all_sf[,"Label"], levels=labels)

# Adipocyte
sce_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.diem_sce.rds") # To get the barcode-rank plot
seur_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
markers_ad <- read.table("results/adpcyte/quantile/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
quant <- readRDS("data/processed/adpcyte/quantile/adpcyte.quantile_counts.rds")
thresh_ad <- quant$thresh

# Mouse brain
sce_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.diem_sce.rds")
seur_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
quant <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.quantile_counts.rds")
thresh_mb <- quant$thresh

# Adipose tissue
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

scel_at <- lapply(lab_ids, function(id) readRDS(paste0("data/processed/atsn/diem/", id, ".diem_sce.rds")) )
names(scel_at) <- lab_ids
seur_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")

thresh_at <- lapply(lab_ids, function(id) {
    quant <- readRDS(paste0("data/processed/atsn/quantile/", id, ".quantile_counts.rds"))
    return(quant$thresh)
})
names(thresh_at) <- lab_ids

#=========================================
# Read in data
#=========================================

scel <- c(sce_ad, sce_mb, scel_at)
names(scel) <- samples

datf <- lapply(1:8, function(i){
    x <- scel[[i]]
    lb <- labels[i]
    s <- samples[i]
    counts <- x@droplet_data[order(x@droplet_data[,"total_counts"], decreasing=TRUE), "total_counts"]
    counts <- sort(counts, decreasing=TRUE)
    ranks <- seq(length(counts))
    counts[duplicated(counts)] <- NA
    ranks[duplicated(counts)] <- NA
    df <- data.frame("Rank"=ranks, "Count"=counts, "Sample" = s, "Label" = lb)
    df <- df[!is.na(df[,2]),,drop=FALSE]
    return(df)
})
datf <- do.call(rbind, datf)
datf[,"Label"] <- factor(datf[,"Label"], levels = labels)

quants <- data.frame("Thresh" = c(thresh_ad, thresh_mb, unlist(thresh_at)), 
                     "Label" = labels)
quants[,"Label"] <- factor(quants[,"Label"], levels = labels)
quants[,"Thresh"] <- log10(quants[,"Thresh"])

p <- ggplot(datf, aes_string(x = "Rank", y = "Count")) + 
geom_point() + 
geom_abline(aes(intercept = Thresh, 
                 slope = 0), colour = "red", linetype = "dashed", data=quants) + 
facet_wrap(~Label, ncol = 2) + 
scale_x_continuous(trans='log10', breaks=set_breaks_10, labels=comma) +
scale_y_continuous(name="Droplet size", trans='log10', labels=comma)  +
theme_minimal() + 
theme(strip.text = element_text(size = 14),
      plot.title=element_text(hjust=0.5),
      axis.text.x=element_text(angle=45, hjust=1))
fn <- "results/plots/barcode_rank.pdf"
ggsave(fn, width = 4, height = 6)
fn <- "results/plots/barcode_rank.jpeg"
ggsave(fn, width = 4, height = 6)

