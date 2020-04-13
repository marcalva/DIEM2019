
setwd("../")

library(diem)
suppressMessages(library(DropletUtils))
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(scales)
source("scripts/common/plotting.R")


# Common function and themes
# Boxplot of MT%
boxplot_mt <- function(x, mt_col="percent.mt", ct_col="RNA_snn_res.0.8", yl = "MT%"){
	require(reshape2)
	df <- data.frame(Mito=x@meta.data[,mt_col], Cluster=x@meta.data[,ct_col])
	df <- reshape2::melt(df)
    df[,"value"] <- df[,"value"] / 100
	p <- ggplot(df, aes(x=Cluster, y=value)) +
	geom_boxplot() +  
    theme_bw() + 
    scale_y_continuous(labels = scales::percent_format(accuracy=1)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 0.9)) + 
	ylab(yl)
	return(p)
}

ct <- theme(legend.position="none",
			text=element_text(size=14),
			plot.title=element_text(size=18, hjust=0.5)
			)

ctl <- theme(text=element_text(size=14),
			 plot.title=element_text(size=18, hjust=0.5)
			 )
yexp <- 1.1


#=========================================
# Read in data
#=========================================

# Adipocyte
seur_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
seur_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
seur_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")

# Adjust Seurat MT% for plotting
ch <- seur_at@meta.data$percent.mt > 10
seur_at@meta.data$percent.mt[ch] <- 10


#=========================================
# Call plot functions
#=========================================

low = "lightgrey"
high = "firebrick3"

p1 <- boxplot_mt(seur_ad) + ct + ggtitle("DiffPA")
p1u <- plot_umap_meta(seur_ad, 
                      size = 0.8, alpha = 0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low=low, high=high, name="MT%") +
ggtitle("DiffPA")



p2 <- boxplot_mt(seur_mb) + ct + ggtitle("Mouse Brain")
p2u <- plot_umap_meta(seur_mb, 
                      size = 0.8, alpha = 0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low=low, high=high, name="MT%") +
ggtitle("Mouse Brain")


p3 <- boxplot_mt(seur_at) + ct + ggtitle("Adipose Tissue") + 
scale_y_continuous(labels = scales::percent_format(accuracy=1), limits = c(0,.1))
p3u <- plot_umap_meta(seur_at, 
                      size = 0.8, alpha=0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low=low, high=high, name="MT%") + 
ggtitle("Adipose Tissue")  

#=========================================
# Put all in 1 plot
#=========================================

dir_plot <- "results/plots/"
dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)


pdfname <- paste0(dir_plot, "quant.MT.pdf")
jpgname <- paste0(dir_plot, "quant.MT.jpeg")
pdf(pdfname, width=10, height=6)
ggarrange(ggarrange(p1u, p2u, p3u, ncol = 3), 
          ggarrange(p1, p2, p3, ncol = 3), 
          labels=c("a", "b"), ncol = 1, font.label = list(size = 16, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))



