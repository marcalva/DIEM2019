
# Figure 5

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

plot_umap_fct <- function(x, names, colname="percent.mt", legend_name="MT%", 
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

# Adipocyte
labl <- "adpcyte"
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

# Read in DIEM SCE
seur_diem_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.seur_obj.rds")
seur_quant_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
seur_ED_ad <- readRDS("data/processed/adpcyte/emptydrops/adpcyte.seur_obj.rds")

# Mouse brain
labl <- "mouse_nuclei_2k"
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

seur_diem_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.seur_obj.rds")
seur_quant_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
seur_ED_mb <- readRDS("data/processed/mouse_nuclei_2k/emptydrops/mouse_nuclei_2k.seur_obj.rds")

# Adipose tissue
labl <- "atsn"
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

seur_diem_at <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_quant_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
seur_ED_at <- readRDS("data/processed/atsn/emptydrops/atsn.seur_obj.rds")

ch <- seur_diem_at@meta.data$percent.mt > 10
seur_diem_at@meta.data$percent.mt[ch] <- 10
ch <- seur_quant_at@meta.data$percent.mt > 10
seur_quant_at@meta.data$percent.mt[ch] <- 10
ch <- seur_ED_at@meta.data$percent.mt > 10
seur_ED_at@meta.data$percent.mt[ch] <- 10

#=========================================
# Main Figure
#=========================================

pu1 <- plot_umap_fct(list(seur_quant_ad, seur_ED_ad, seur_diem_ad), methd_names, 
					 colname="percent.mt", legend_name="MT%", size=0.5) + ggtitle("DiffPA")
pu2 <- plot_umap_fct(list(seur_quant_mb, seur_ED_mb, seur_diem_mb), methd_names, 
					 colname="percent.mt", legend_name="MT%", size=0.5) + ggtitle("Mouse Brain")
pu3 <- plot_umap_fct(list(seur_quant_at, seur_ED_at, seur_diem_at), methd_names, 
					 colname="percent.mt", legend_name="MT%", size=0.5) + ggtitle("Adipose Tissue")

dir_plot <- paste0("results/plots/")
pdfname <- paste0(dir_plot, "umap_mt.pdf")
jpgname <- paste0(dir_plot, "umap_mt.jpeg")
pdf(pdfname, width=12, height=12)
ggarrange(pu1, pu2, pu3, ncol=1, labels=c("a", "b", "c"), 
          font.label = list(size = 20, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))

