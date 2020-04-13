
setwd("../")

library(diem)
library(ggplot2)
library(scales)
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

plot_umap_pct <- function(x, names, colname="SpliceFrctn", 
                          legend_name="Percent\nreads\nspliced",
                          low = "lightgrey", high = "firebrick3", color_breaks=waiver(),
                          size=0.2, alpha=alpha, order_col = TRUE){
    dfl <- lapply(1:length(x), function(i) {
                  data.frame(Feat=x[[i]]@meta.data[,colname],
                             UMAP1=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_1"],
                             UMAP2=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_2"],
                             Method=names[[i]])
                          })
    df <- do.call(rbind, dfl)
    if (order_col) df <- df[order(df$Feat, decreasing=FALSE),,drop=FALSE]

    if (any(df[,"Feat"] > 1)){
        df[,"Feat"] <- df[,"Feat"] / 100
    }
    p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Feat)) +
    geom_point(size=size, alpha=0.8) + theme_bw() +
    facet_wrap(~Method, ncol=1, scale="free") +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), 
          strip.text = element_text(size = 14),
          strip.background = element_blank()) +
    scale_colour_gradient(labels = scales::percent, 
                          low=low, high=high,
                          name=legend_name,
                          breaks=color_breaks)
    return(p)
}


#=========================================
#=========================================

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
expr <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- expr
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples
eLabels <- c("DiffPA", "Mouse Brain", rep("Adipose Tissue", 6))
names(eLabels) <- samples

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
seur_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
seur_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
seur_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")

seur_l <- list(seur_ad, seur_mb, seur_at)

#=========================================
# Plot UMAP
#=========================================

dir_plot <- "results/plots/"
p <- plot_umap_pct(seur_l, names=c("DiffPA", "Mouse Brain", "Adipose Tissue"))

outfn <- paste0(dir_plot, "quantile.umap.sf.jpeg")
ggsave(outfn, width = 3, height = 6)
outfn <- paste0(dir_plot, "quantile.umap.sf.pdf")
ggsave(outfn, width = 3, height = 6)


