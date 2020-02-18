
# Splice fractions

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

#=========================================
# Functions
#=========================================

ct <- theme(text=element_text(size=16),
			plot.title=element_text(size=18, hjust=0.5, face="bold")
			)


yexp <- 1.1

plot_umap_fct <- function(x, names, colname="SpliceFrctn", legend_name="Fraction Spliced", 
						  color_limits=NULL, color_breaks=waiver(),
						  size=1, alpha=alpha, order_col = TRUE){
	dfl <- lapply(1:length(x), function(i) {
                  data.frame(Mito=x[[i]]@meta.data[,colname], 
                             UMAP1=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_1"], 
                             UMAP2=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_2"], Method=names[[i]])
                          })
	df <- do.call(rbind, dfl)
    if (order_col) df <- df[order(df$Mito, decreasing=FALSE),,drop=FALSE]

	p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Mito)) + 
	geom_point(size=size, alpha=0.8) + theme_bw() + 
	facet_wrap(~Method, ncol=3, scale="free") + 
	theme(text=element_text(size=16), 
		  axis.text=element_blank(), 
		  axis.ticks=element_blank(), 
		  plot.title=element_text(hjust=0.5), 
          panel.grid=element_blank()) + 
    scale_colour_gradient(low="gray90", high="red3", 
                          name=legend_name, 
                          limits=color_limits, 
                          breaks=color_breaks)
	#scale_color_distiller(palette = "Spectral", name=legend_name, 
	#					  limits=color_limits, breaks=color_breaks) 

	return(p)
}

#=========================================
#=========================================

methd_names <- c("Quantile", "EmptyDrops", "DIEM")

#=========================================
# Read in data
#=========================================

post_merged <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
pre_merged <- readRDS("data/processed/atsn/diem_integrated.fltr.1/atsn.seur_obj.rds")



#=========================================
# UMAP plots
#=========================================

t1 <- paste0("Adipose Tissue merged\nafter DIEM (n = ", ncol(post_merged), ")")
t2 <- paste0("Adipose Tissue merged\nbefore DIEM (n = ", ncol(pre_merged), ")")
p1 <- plot_umap_labels(post_merged) + ggtitle(t1) + ct
p2 <- plot_umap_labels(pre_merged) + ggtitle(t2) + ct


dplot <- "results/plots/"
pdfname <- paste0(dplot, "umap.pre_post.pdf")
jpgname <- paste0(dplot, "umap.pre_post.jpeg")
pdf(pdfname, width=12,height=6)
ggarrange(p1, p2, nrow = 1, labels=c("a", "b"),
        font.label = list(size = 18, face="bold"))
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))





