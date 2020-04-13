
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
                             UMAP2=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_2"], 
                             Method=names[[i]])
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

seur_ad = list("DIEM" = seur_diem_ad, "EmptyDrops" = seur_ED_ad, "Quantile" = seur_quant_ad)
seur_mb = list("DIEM" = seur_diem_mb, "EmptyDrops" = seur_ED_mb, "Quantile" = seur_quant_mb)
seur_at = list("DIEM" = seur_diem_at, "EmptyDrops" = seur_ED_at, "Quantile" = seur_quant_at)


#=========================================
# Main Figure
#=========================================

pu1 <- plot_umap_fct(list(seur_quant_ad, seur_ED_ad, seur_diem_ad), methd_names, 
					 colname="SpliceFrctn", legend_name="Fraction Spliced", size=0.5) + ggtitle("DiffPA")
pu2 <- plot_umap_fct(list(seur_quant_mb, seur_ED_mb, seur_diem_mb), methd_names, 
					 colname="SpliceFrctn", legend_name="Fraction Spliced", size=0.5) + ggtitle("Mouse Brain")
pu3 <- plot_umap_fct(list(seur_quant_at, seur_ED_at, seur_diem_at), methd_names, 
					 colname="SpliceFrctn", legend_name="Fraction Spliced", size=0.5) + ggtitle("Adipose Tissue")


dir_plot <- paste0("results/plots/")
pdfname <- paste0(dir_plot, "SpliceFrctn.umap.pdf")
jpgname <- paste0(dir_plot, "SpliceFrctn.umap.jpeg")
pdf(pdfname, width=10, height=12)
ggarrange(pu1, pu2, pu3, nrow=3)
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))



#=========================================
# Proportion of droplets with high spliced reads
#=========================================

above_sd <- function(seur, trait = "SpliceFrctn", thresh = NULL){
    ta <- seur@meta.data[,trait]
    if (is.null(thresh)){
        thresh <- 2 * sd(ta, na.rm = T) + mean(ta, na.rm = T) 
    }
    ab <- table(ta > thresh)
    return(ab[["TRUE"]] / length(ta))
}


a_ad <- sapply(seur_ad, above_sd, thresh = 50)
a_mb <- sapply(seur_mb, above_sd, thresh = 50)
a_at <- sapply(seur_at, above_sd, thresh = 50)

datf <- data.frame("DiffPA" = a_ad, "Mouse brain" = a_mb, "Adipose Tissue" = a_at)
rownames(datf) <- c("DIEM", "EmptyDrops", "Quantile")
datfm <- reshape2::melt(as.matrix(datf))
colnames(datfm) <- c("Method", "DataSet", "PD")
datfm$PD <- 100 * datfm$PD

# Barplot of total MT markers
labeler <- c("DiffPA", "Mouse brain", "Adipose Tissue")
names(labeler) <- c("DiffPA", "Mouse.brain", "Adipose.Tissue")
p <- ggplot(datfm, aes(x=Method, y=PD, fill=Method)) + 
geom_bar(stat="identity", color="black", position=position_dodge()) + 
facet_wrap(~DataSet, scale="free_y", ncol = 1, labeller = labeller(DataSet = labeler)) + 
theme_minimal() + ylab("Background droplets (percent of filtered)") + 
theme(legend.position = "none", 
      axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = .8),
      axis.title.x = element_blank(),
      text=element_text(size=18),
      plot.title=element_text(hjust=0.5), 
      panel.grid.major.x=element_blank())

pdf("results/plots/prop_sf_clusters.pdf", width=3, height=7)
p
dev.off()


