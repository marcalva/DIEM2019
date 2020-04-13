
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
# Read in data
#=========================================

dir_plot <- "results/plots/"

seur_diem_at.clean <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_diem_at.debris <- readRDS("data/processed/atsn/diem_debris/atsn.seur_obj.rds")

seur_quant_at.clean <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
seur_quant_at.debris <- readRDS("data/processed/atsn/quantile_debris/atsn.seur_obj.rds")

seur_ED_at.clean <- readRDS("data/processed/atsn/emptydrops/atsn.seur_obj.rds")
seur_ED_at.debris <- readRDS("data/processed/atsn/emptydrops_debris/atsn.seur_obj.rds")


#=========================================
#=========================================

w = 6.5
h = 6

pu1 <- plot_umap_labels_r(seur_diem_at.clean, 
                          ct_col = "CellType", 
                          size = 0.4,
                          label_size = 6, 
                          labels = TRUE) + 
ggtitle(paste0("DIEM passing\n(n=", as.character(ncol(seur_diem_at.clean)), ")")) + 
theme(text=element_text(size=18), plot.title=element_text(hjust=0.5))

outfn <- paste0(dir_plot, "umap.pass.diem.clust.pdf")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "umap.pass.diem.clust.jpeg")
ggsave(outfn, width = w, height = h)




pu2 <- plot_umap_labels_r(seur_ED_at.clean, 
                          ct_col = "CellType", 
                          size = 0.4, 
                          label_size = 6,
                          labels = TRUE) + 
ggtitle(paste0("EmptyDrops passing\n(n=", as.character(ncol(seur_ED_at.clean)), ")")) + 
theme(text=element_text(size=18), plot.title=element_text(hjust=0.5))

outfn <- paste0(dir_plot, "umap.pass.emptydrops.clust.pdf")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "umap.pass.emptydrops.clust.jpeg")
ggsave(outfn, width = w, height = h)




pu3 <- plot_umap_labels_r(seur_quant_at.clean, 
                          ct_col = "CellType", 
                          label_size = 6,
                          size = 0.4,
                          labels = TRUE) + 
ggtitle(paste0("Quantile passing\n(n=", as.character(ncol(seur_quant_at.clean)), ")")) + 
theme(text=element_text(size=18), plot.title=element_text(hjust=0.5))

outfn <- paste0(dir_plot, "umap.pass.quantile.clust.pdf")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "umap.pass.quantile.clust.jpeg")
ggsave(outfn, width = w, height = h)

