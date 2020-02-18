
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

#=========================================
# Plot
#=========================================

l_ad <- list("DIEM" = colnames(seur_diem_ad), 
             "Quantile" = colnames(seur_quant_ad), 
             "EmptyDrops" = colnames(seur_ED_ad))

pdf("results/plots/drp_ovrlp.ad.pdf", width = 3.5, height = 4)
par(mar = c(0, 0, 0, 0))
venn(l_ad)
dev.off()

l_at <- list("DIEM" = colnames(seur_diem_at), 
             "Quantile" = colnames(seur_quant_at), 
             "EmptyDrops" = colnames(seur_ED_at))

pdf("results/plots/drp_ovrlp.at.pdf", width = 3.5, height = 4)
par(mar = c(0, 0, 0, 0))
venn(l_at)
dev.off()

l_mb <- list("DIEM" = colnames(seur_diem_mb), 
             "Quantile" = colnames(seur_quant_mb), 
             "EmptyDrops" = colnames(seur_ED_mb))

pdf("results/plots/drp_ovrlp.mb.pdf", width = 3.5, height = 4)
par(mar = c(0, 0, 0, 0))
venn(l_mb)
dev.off()

