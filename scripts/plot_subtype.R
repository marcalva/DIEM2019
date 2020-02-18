
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

ct <- theme(text=element_text(size=16),
			plot.title=element_text(size=18, hjust=0.5, face="bold")
			)
yexp <- 1.1


boxplot_ft_fct <- function(x, names, ct_col="RNA_snn_res.0.8",
                           colname="percent.mt", ylab="MT%",
                           color_breaks=waiver(),
                           size=1){
    dfl <- lapply(1:length(x), function(i) data.frame(Mito=x[[i]]@meta.data[,colname], Cluster=x[[i]]@meta.data[,ct_col], Method=names[[i]]))
    df <- do.call(rbind, dfl)
    df <- reshape2::melt(df)

    p <- ggplot(df, aes(x=Cluster, y=value)) +
    geom_boxplot(outlier.shape = NA) + theme_bw() +
    facet_wrap(~Method, scales="free_x") +
    ylab(ylab) +
    theme(legend.position="none",
          text=element_text(size=16),
          axis.text.x=element_text(size=10),
          plot.title=element_text(hjust=0.5))
    return(p)
}

plot_umap_fct <- function(x, names, colname="percent.mt", legend_name="MT%",
                          color_limits=NULL, color_breaks=waiver(),
                          size=1){
    dfl <- lapply(1:length(x), function(i) {
                  data.frame("Mito"=x[[i]]@meta.data[,colname],
                             "UMAP1"=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_1"],
                             "UMAP2"=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_2"],
                             "Method"=names[[i]])})
    df <- do.call(rbind, dfl)

    p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Mito)) +
    geom_point(size=size) + theme_bw() +
    facet_wrap(~Method, ncol=3, scale="free") +
    theme(text=element_text(size=16),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title=element_text(hjust=0.5)) +
    scale_color_distiller(palette = "Spectral", name=legend_name,
                          limits=color_limits, breaks=color_breaks)

    return(p)
}

p_blank <- ggplot() + theme_void()


#=========================================
#=========================================

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

library(ggplot2)

pu1 <- plot_umap_gene(seur_diem_at, "NACA2", size = 2)
pu2 <- plot_umap_gene(seur_diem_at, "NACA", size = 2)
pu3 <- plot_umap_gene(seur_diem_at, "CACNA2D1", size = 2)

pb1 <- boxplot_gene(seur_diem_at, "NACA2") + ct
pb2 <- boxplot_gene(seur_diem_at, "NACA") + ct
pb3 <- boxplot_gene(seur_diem_at, "CACNA2D1") + ct

pe1 <- boxplot_gene(seur_ED_at, "NACA2") + ct
pe2 <- boxplot_gene(seur_ED_at, "NACA") + ct
pe3 <- boxplot_gene(seur_ED_at, "CACNA2D1") + ct

pq1 <- boxplot_gene(seur_quant_at, "NACA2") + ct
pq2 <- boxplot_gene(seur_quant_at, "NACA") + ct
pq3 <- boxplot_gene(seur_quant_at, "CACNA2D1") + ct

g <- "GPAM"
pb4 <- boxplot_gene(seur_diem_at, g) + ct
pq4 <- boxplot_gene(seur_quant_at, g) + ct
pe4 <- boxplot_gene(seur_ED_at, g) + ct

fl <- list(size = 20, face="bold")
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "subtype.pdf")
jpgname <- paste0(dir_plot, "subtype.jpeg")
pdf(pdfname, width=15, height=15)
ggarrange(
    ggarrange(pq1, pq2, pq3, pq4, ncol = 1),
    ggarrange(pe1, pe2, pe3, pe4, ncol = 1),
    ggarrange(pb1, pb2, pb3, pb4, ncol = 1),
    nrow = 1)

dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))








