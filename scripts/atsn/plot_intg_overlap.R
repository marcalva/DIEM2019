
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
library(viridis)
library(reshape2)
source("scripts/common/standard_seurat.R")
source("scripts/common/plotting.R")
source("scripts/common/overlap_graph.R")


#=========================================
# Functions
#=========================================
ct <- theme(text = element_text(size=16))

p_blank <- ggplot() + theme_void()

h_theme <- theme(panel.grid=element_blank(),
                 panel.border = element_blank(),
                 text=element_text(size=16))


plot_overlap <- function(datf, col1 = "DIEM Cluster", col2 = "Seurat Cluster", title = ""){
    datfm <- reshape2::melt(as.matrix(datf * 100))
    datfm[,1] <- factor(datfm[,1])
    datfm[,2] <- factor(datfm[,2])
    p <- ggplot(datfm, aes_string(x = "Var2", y = "Var1")) + 
    geom_tile(aes_string(fill="value"), colour="white") + 
    scale_fill_viridis(name = "Percent\nOverlap", limits = c(0,100)) + 
    theme_minimal() + 
    h_theme + 
    ggtitle(title) + 
    scale_x_discrete(position = "top") + 
    ylab(col1) + xlab(col2)
    return(p)
}


#=========================================
#=========================================

#=========================================
# Read in data
#=========================================


dp <- "data/processed/atsn/diem/"
seur_raw <- readRDS(paste0(dp, "atsn.seur_obj.rds"))
seur_intg <- readRDS(paste0(dp, "atsn.seur_intg.seur_obj.rds"))

md_raw <- seur_raw@meta.data
md_intg <- seur_intg@meta.data


md_all <- cbind(md_raw, md_intg[rownames(md_raw),])
tb <- table(md_all$integrated_snn_res.0.8, md_all$RNA_snn_res.0.8)
tb <- as.data.frame.matrix(tb)
pct <- sweep(tb, 2, colSums(tb), "/")

p1 <- plot_umap_labels(seur_raw, ct_col = "RNA_snn_res.0.8") + 
    ggtitle("Raw merge") + theme(plot.title = element_text(size = 18, hjust = 0.5))
p2 <- plot_umap_labels(seur_intg, ct_col = "integrated_snn_res.0.8") + 
    ggtitle("Seurat integration") + theme(plot.title = element_text(size = 18, hjust = 0.5))
p3 <- plot_overlap(pct, col2 = "Raw Merge", col1 = "Seurat Integration")

dir_plot <- "results/atsn/diem/plots/"
dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "raw_v_seur_intg.pdf")
jpgname <- paste0(dir_plot, "raw_v_seur_intg.jpeg")
pdf(pdfname, width=14, height=5)
ggarrange(p1, p2, p3, labels=c("a", "b", "c"), ncol = 3, font.label = list(size = 24, face="bold"))
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))

