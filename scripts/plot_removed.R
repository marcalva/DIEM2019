
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

get_n_s_mt <- function(markers){
	markers <- markers[markers$p_val_adj < 0.05,]
	nmt <- length(grep("^mt-", markers[,"gene"], ignore.case=TRUE))
	return(nmt)
}

ct <- theme(text=element_text(size=16),
			plot.title=element_text(size=18, hjust=0.5, face="bold")
			)
yexp <- 1.1


lincs <- c("^malat1$", "^neat1$", "^xist$", "^meg3$")
get_nuc_linc <- function(seur, patterns=lincs, colname="percent.linc"){
	genes <- c()
	for (i in patterns){
		genes <- c(genes, grep(i, rownames(seur), ignore.case=TRUE, value=TRUE))
	}
	genes <- unique(genes)
	seur@meta.data[,colname] <- 100 * (Matrix::colSums(seur@assays$RNA@counts[genes,,drop=FALSE])/Matrix::colSums(seur@assays$RNA@counts))
	return(seur)
}

boxplot_ft_fct <- function(x, names, ct_col="RNA_snn_res.0.8", 
						   colname="percent.mt", ylab="MT%", 
						   color_breaks=waiver(), scales = "free_x", ncol = 1,
						   size=1){
	dfl <- lapply(1:length(x), function(i) data.frame(Mito=x[[i]]@meta.data[,colname], Cluster=x[[i]]@meta.data[,ct_col], Method=names[[i]]))
	df <- do.call(rbind, dfl)
	df <- reshape2::melt(df)
    names <- factor(names)
    df[,"Method"] <- factor(df[,"Method"], levels = levels(names))

	p <- ggplot(df, aes(x=Cluster, y=value)) +
	geom_boxplot(outlier.shape = NA) + theme_bw() + 
	facet_wrap(~Method, scales=scales, ncol = ncol) + 
	ylab(ylab) + 
	theme(legend.position="none", 
		  text=element_text(size=16),
		  axis.text.x=element_text(size=10), 
		  plot.title=element_text(hjust=0.5)) 
	return(p)
}

plot_umap_fct <- function(x, names, colname="percent.mt", legend_name="MT%", 
						  color_limits=NULL, color_breaks=waiver(),
						  size=1, scales = "free"){
	dfl <- lapply(1:length(x), function(i) {
                  data.frame("Mito"=x[[i]]@meta.data[,colname], 
                             "UMAP1"=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_1"], 
                             "UMAP2"=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_2"], 
                             "Method"=names[[i]])})
	df <- do.call(rbind, dfl)

	p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Mito)) + 
	geom_point(size=size) + theme_bw() + 
	facet_wrap(~Method, ncol=3, scale=scales) + 
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

# Variables
methd_names <- c("Quantile", "EmptyDrops", "DIEM")
methd <- factor(c(methd_names, methd_names, methd_names), levels=c("Quantile", "EmptyDrops", "DIEM"))


#=========================================
# Read in data
#=========================================

# Adipose tissue
labl <- "atsn"
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")



seur_diem_at.clean <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_diem_at.debris <- readRDS("data/processed/atsn/diem_debris/atsn.seur_obj.rds")

seur_quant_at.clean <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
seur_quant_at.debris <- readRDS("data/processed/atsn/quantile_debris/atsn.seur_obj.rds")

seur_ED_at.clean <- readRDS("data/processed/atsn/emptydrops/atsn.seur_obj.rds")
seur_ED_at.debris <- readRDS("data/processed/atsn/emptydrops_debris/atsn.seur_obj.rds")


#=========================================
#=========================================

methd_names <- c("Removed", "Passing")

pu1 <- plot_umap_labels(seur_diem_at.debris) + 
ggtitle(paste0("DIEM removed\n(n=", as.character(ncol(seur_diem_at.debris)), ")")) + 
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(hjust=0.5))

pu2 <- plot_umap_labels(seur_ED_at.debris) + 
ggtitle(paste0("EmptyDrops removed\n(n=", as.character(ncol(seur_ED_at.debris)), ")")) + 
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(hjust=0.5))

pu3 <- plot_umap_labels(seur_quant_at.debris) + 
ggtitle(paste0("Quantile removed\n(n=", as.character(ncol(seur_quant_at.debris)), ")")) + 
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(hjust=0.5))


pa1 <- boxplot_ft_fct(list(seur_diem_at.debris), methd_names, 
                      colname="SpliceFrctn", ylab = "% reads\nspliced") + ggtitle("DIEM")
pa2 <- boxplot_ft_fct(list(seur_ED_at.debris), methd_names, 
                      colname="SpliceFrctn", ylab = "% reads\nspliced") + ggtitle("EmptyDrops")
pa3 <- boxplot_ft_fct(list(seur_quant_at.debris), methd_names, 
                      colname="SpliceFrctn", ylab = "% reads\nspliced") + ggtitle("Quantile")

methd_names <- c("DIEM", "EmptyDrops", "Quantile")
methd_names <- factor(methd_names, levels = methd_names)
sl <- list("DIEM" = seur_diem_at.debris, "EmptyDrops" = seur_ED_at.debris, "Quantile" = seur_quant_at.debris)
pbox <- boxplot_ft_fct(sl, methd_names, scales = "free_x",  
                       colname="SpliceFrctn", ylab = "% reads\nspliced")


g <- "PLIN1"
pg1 <- plot_umap_gene(seur_diem_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(DIEM removed)")) + ct
pg2 <- plot_umap_gene(seur_ED_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(EmptyDrops removed)")) + ct
pg3 <- plot_umap_gene(seur_quant_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(Quantile removed)")) + ct

g <- "RBPJ"
pj1 <- plot_umap_gene(seur_diem_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(DIEM removed)")) + ct
pj2 <- plot_umap_gene(seur_ED_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(EmptyDrops removed)")) + ct
pj3 <- plot_umap_gene(seur_quant_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(Quantile removed)")) + ct


g <- "VWF"
pk1 <- plot_umap_gene(seur_diem_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(DIEM removed)")) + ct
pk2 <- plot_umap_gene(seur_ED_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(EmptyDrops removed)")) + ct
pk3 <- plot_umap_gene(seur_quant_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(Quantile removed)")) + ct

g <- "CFD"
pl1 <- plot_umap_gene(seur_diem_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(DIEM removed)")) + ct
pl2 <- plot_umap_gene(seur_ED_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(EmptyDrops removed)")) + ct
pl3 <- plot_umap_gene(seur_quant_at.debris, g, size = 1) + ggtitle(paste0(g, " expression\n(Quantile removed)")) + ct


# Arrange into plot
fl <- list(size = 20, face="bold")
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "removed.pdf")
jpgname <- paste0(dir_plot, "removed.jpeg")
pdf(pdfname, width=20, height=10)
ggarrange(ggarrange(pu1, pu2, pu3, ncol = 1), 
          ggarrange(pg1, pg2, pg3, ncol = 1), 
          ggarrange(pj1, pj2, pj3, ncol = 1),
          ggarrange(pk1, pk2, pk3, ncol = 1),
          ggarrange(pl1, pl2, pl3, ncol = 1), 
          pbox, 
          # ggarrange(pa1, pa2, pa3, ncol = 1),
          nrow = 1, labels=c("a", "b", "c", "d", "e", "f"), font.label = fl)
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


