
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
    
    names <- factor(names)
    datf[,"Method"] <- factor(datf[,"Method"], 
                              levels = levels(names))

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
          axis.title.x = element_blank(), 
          axis.title = element_text(size = 24),
          text = element_text(size=22))
    return(p)
}


#=========================================
# Read in data
#=========================================

dir_plot <- "results/plots/"

# Adipocyte
seur_diem_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.seur_obj.rds")
seur_quant_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
seur_ED_ad <- readRDS("data/processed/adpcyte/emptydrops/adpcyte.seur_obj.rds")

# Mouse brain
seur_diem_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.seur_obj.rds")
seur_quant_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
seur_ED_mb <- readRDS("data/processed/mouse_nuclei_2k/emptydrops/mouse_nuclei_2k.seur_obj.rds")

# Adipose tissue
seur_diem_at <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_quant_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
seur_ED_at <- readRDS("data/processed/atsn/emptydrops/atsn.seur_obj.rds")

# fresh 68k
seur_diem_b <- readRDS("data/processed/fresh_68k/diem/fresh_68k.seur_obj.rds")
seur_ED_b <- readRDS("data/processed/fresh_68k/emptydrops/fresh_68k.seur_obj.rds")

seur_ad <- list("DIEM" = seur_diem_ad, "EmptyDrops" = seur_ED_ad, "Quantile" = seur_quant_ad)
seur_mb <- list("DIEM" = seur_diem_mb, "EmptyDrops" = seur_ED_mb, "Quantile" = seur_quant_mb)
seur_at <- list("DIEM" = seur_diem_at, "EmptyDrops" = seur_ED_at, "Quantile" = seur_quant_at)
seur_b <- list("DIEM" = seur_diem_b, "EmptyDrops" = seur_ED_b)

#=========================================
# Plot
#=========================================

w <- 6
h <- 15

hnames <- c("DIEM clusters", "EmptyDrops clusters", "Quantile clusters")

# DiffPA
p <- boxplot_ft_l(seur_ad, 
                  hnames, 
                  ct_col = "RNA_snn_res.0.8",
                  colname = "SpliceFrctn", 
                  ylab = "Percent reads spliced")

outfn <- paste0(dir_plot, "sf.clust.method.pass.DiffPA.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "sf.clust.method.pass.DiffPA.pdf")
ggsave(outfn, width = w, height = h)


# Mouse brain
p <- boxplot_ft_l(seur_mb, 
                  hnames, 
                  ct_col = "RNA_snn_res.0.8",
                  colname = "SpliceFrctn", 
                  ylab = "Percent reads spliced")

outfn <- paste0(dir_plot, "sf.clust.method.pass.MouseBrain.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "sf.clust.method.pass.MouseBrain.pdf")
ggsave(outfn, width = w, height = h)


# Adipose
p <- boxplot_ft_l(seur_at, 
                  hnames, 
                  ct_col = "CellType",
                  colname = "SpliceFrctn", 
                  ylab = "Percent reads spliced")
p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 0.9))

outfn <- paste0(dir_plot, "sf.clust.method.pass.atsn.jpeg")
ggsave(outfn, width = 10, height = h)
outfn <- paste0(dir_plot, "sf.clust.method.pass.atsn.pdf")
ggsave(outfn, width = 10, height = h)


