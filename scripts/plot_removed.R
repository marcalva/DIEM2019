
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
ggtitle(paste0("DIEM removed (n=", as.character(ncol(seur_diem_at.debris)), ")")) + 
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(hjust=0.5))

pu2 <- plot_umap_labels(seur_ED_at.debris) + 
ggtitle(paste0("EmptyDrops removed (n=", as.character(ncol(seur_ED_at.debris)), ")")) + 
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(hjust=0.5))

pu3 <- plot_umap_labels(seur_quant_at.debris) + 
ggtitle(paste0("Quantile removed (n=", as.character(ncol(seur_quant_at.debris)), ")")) + 
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(hjust=0.5))


pa1 <- boxplot_ft_fct(list(seur_diem_at.debris, seur_diem_at.clean), methd_names) + ylim(0,10) + ggtitle("DIEM")
pa2 <- boxplot_ft_fct(list(seur_ED_at.debris, seur_ED_at.clean), methd_names) + ylim(0,10) + ggtitle("EmptyDrops")
pa3 <- boxplot_ft_fct(list(seur_quant_at.debris, seur_quant_at.clean), methd_names) + ylim(0,10) + ggtitle("Quantile")

# Arrange into plot
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "removed.pdf")
jpgname <- paste0(dir_plot, "removed.jpeg")
pdf(pdfname, width=18, height=10)
ggarrange(pu1, pa1, pu2,  pa2, pu3,  pa3, nrow=3, 
          ncol=2, labels=c("a", "", "b", "", "c", ""), 
          font.label = list(size = 20, face="bold"), 
          widths=c(1,1.5))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


