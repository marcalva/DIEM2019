
setwd("../")

library(diem)
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
						   color_breaks=waiver(), 
						   size=1){
	dfl <- lapply(1:length(x), function(i) {
                  data.frame(Mito=x[[i]]@meta.data[,colname], 
                             Cluster=as.character(x[[i]]@meta.data[,ct_col]), 
                             Method=names[[i]])
                           })
	df <- do.call(rbind, dfl)
	df <- reshape2::melt(df)
    df$Method <- factor(df$Method, levels=names)

	p <- ggplot(df, aes(x=Cluster, y=value)) +
	geom_boxplot(outlier.shape = NA) + theme_bw() + 
	facet_wrap(~Method, scales = "free_x") + 
	ylab(ylab) + 
	theme(legend.position="none", 
		  text=element_text(size=12),
		  axis.text.x=element_text(size=8), 
		  plot.title=element_text(hjust=0.5)) 
	return(p)
}

plot_umap_fct <- function(x, names, colname="percent.mt", legend_name="MT%", 
						  color_limits=NULL, color_breaks=waiver(),
						  size=1){
	dfl <- lapply(1:length(x), function(i){
                  data.frame(Mito=x[[i]]@meta.data[,colname], 
                             UMAP1=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_1"], 
                             UMAP2=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_2"], 
                             Method=names[[i]])
                          })
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

# Variables
methd_names <- c("Quantile", "EmptyDrops", "DIEM")
methd <- factor(methd_names, levels=c("Quantile", "EmptyDrops", "DIEM"))

dir_plot <- "results/plots/"
w = 6
h = 3

#=========================================
# Read in data
#=========================================

# Read in DIEM SCE
sce_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.diem_sce.rds")
seur_diem_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.seur_obj.rds")
seur_quant_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
seur_ED_ad <- readRDS("data/processed/adpcyte/emptydrops/adpcyte.seur_obj.rds")
markers_diem_ad <- read.table("results/adpcyte/diem/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_quant_ad <- read.table("results/adpcyte/quantile/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_ED_ad <- read.table("results/adpcyte/emptydrops/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)

seur_quant_ad <- get_nuc_linc(seur_quant_ad)
seur_ED_ad <- get_nuc_linc(seur_ED_ad)
seur_diem_ad <- get_nuc_linc(seur_diem_ad)

pm <- boxplot_ft_fct(list(seur_quant_ad, seur_ED_ad, seur_diem_ad), 
                     names=methd, colname="MALAT1", ylab="MALAT1%") + ggtitle("DiffPA")
outfn <- paste0(dir_plot, "diffpa.cluster.malat1.pdf")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "diffpa.cluster.malat1.jpeg")
ggsave(outfn, width = w, height = h)

pu <-boxplot_ft_fct(list(seur_quant_ad, seur_ED_ad, seur_diem_ad), 
                    names=methd, colname="nCount_RNA", ylab="Total UMI counts") + coord_trans(y="log10")
outfn <- paste0(dir_plot, "diffpa.cluster.logUMI.pdf")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "diffpa.cluster.logUMI.jpeg")
ggsave(outfn, width = w, height = h)


