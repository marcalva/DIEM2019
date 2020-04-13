
setwd("../")

library(diem)
library(Seurat)
suppressMessages(library(DropletUtils))
library(ggplot2)
library(gplots)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
source("scripts/common/plotting.R")
source("scripts/common/overlap_graph.R")



#=========================================
# Functions
#=========================================

boxplot_ft_l <- function(x, 
                         names, 
                         ct_col="RNA_snn_res.0.8", 
                         colname="percent.mt", 
                         ylab="MT%", 
                         color_breaks=waiver(), 
                         scales = "free_x", 
                         ncol = 1,
                         size=1){
    datfl <- lapply(1:length(x), function(i) {
                    data.frame(Feat = x[[i]]@meta.data[,colname], 
                               Cluster = x[[i]]@meta.data[,ct_col], 
                             Method = names[[i]])
                           })
    datf <- do.call(rbind, datfl)
    datf <- reshape2::melt(datf)
    datf[,"value"] <- datf[,"value"]/100
    
    datf[,"Method"] <- factor(datf[,"Method"], 
                              levels = names)

    clust <- factor(datf[,"Cluster"])
    datf[,"Cluster"] <- factor(clust, levels = sort(levels(clust)))

    p <- ggplot(datf, aes(x=Cluster, y=value)) +
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    facet_wrap(~Method, scales=scales, ncol = ncol) + 
    scale_y_continuous(labels = scales::percent, 
                       breaks = c(0, .25, .5, .75, 1),
                       limits = c(0, 1)) +
    ylab(ylab) + 
    theme(legend.position = "none", 
          axis.title = element_text(size = 24),
          axis.text.x = element_text(size=12, angle=30, hjust=0.9), 
          strip.text = element_text(size = 18), 
          strip.background = element_blank(), 
          text = element_text(size=16))
    return(p)
}


#=========================================
# Read in data
#=========================================

dir_plot <- "results/plots/"

# Adipocyte
seur_quant_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")

# Mouse brain
seur_quant_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")

# Adipose tissue
seur_quant_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")

seur_l <- list("DiffPA" = seur_quant_ad, 
               "Mouse Brain" = seur_quant_mb, 
               "Adipose Tissue" = seur_quant_at)

#=========================================
# Plot
#=========================================

hnames <- c("DiffPA", "Mouse Brain", "Adipose Tissue")

# DiffPA
p <- boxplot_ft_l(seur_l, 
                  hnames, 
                  ncol = 3, 
                  ct_col = "RNA_snn_res.0.8",
                  colname = "SpliceFrctn", 
                  ylab = "Percent reads spliced")

outfn <- paste0(dir_plot, "quantile.boxplot.sf.jpeg")
ggsave(outfn, width = 10, height = 4)
outfn <- paste0(dir_plot, "quantile.boxplot.sf.pdf")
ggsave(outfn, width = 10, height = 4)


