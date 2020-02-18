
setwd("../")

library(diem)
suppressMessages(library(DropletUtils))
library(ggplot2)
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
	p <- ggplot(df, aes(x=Cluster, y=value)) +
	geom_boxplot() + theme_minimal() + 
	ylab(yl)
	return(p)
}



ct <- theme(legend.position="none",
			text=element_text(size=16),
			plot.title=element_text(size=18, hjust=0.5, face="bold")
			)

ctl <- theme(text=element_text(size=16),
			 plot.title=element_text(size=18, hjust=0.5, face="bold")
			 )
yexp <- 1.1


#=========================================
# Read in data
#=========================================

# Adipocyte
labl <- "adpcyte"
dp <- paste0("data/processed/", labl, "/")
sce_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.diem_sce.rds") # To get the barcode-rank plot
seur_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
markers_ad <- read.table("results/adpcyte/quantile/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
markers_ad <- markers_ad[markers_ad$p_val_adj < 0.05,]
clusters_ad <- as.character(sort(as.numeric(unique(markers_ad$cluster))))
keep_ad <- sapply(clusters_ad, function(x) length(grep("^mt-", markers_ad[markers_ad$cluster == x,"gene"], ignore.case=TRUE)) != 0)

quant <- readRDS("data/processed/adpcyte/quantile/adpcyte.quantile_counts.rds")
thresh_ad <- quant$thresh

# Mouse brain
labl <- "mouse_nuclei_2k"
dp <- paste0("data/processed/", labl, "/")
sce_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.diem_sce.rds")
seur_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
markers_mb <- read.table("results/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_markers.txt",
					  header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
markers_mb <- markers_mb[markers_mb$p_val_adj < 0.05,]
clusters_mb <- as.character(sort(as.numeric(unique(markers_mb$cluster))))
keep_mb <- sapply(clusters_mb, function(x) length(grep("^mt-", markers_mb[markers_mb$cluster == x,"gene"], ignore.case=TRUE)) != 0)

quant <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.quantile_counts.rds")
thresh_mb <- quant$thresh

# Adipose tissue
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

labl <- "atsn"
dp <- paste0("data/processed/", labl, "/")
scel_at <- lapply(lab_ids, function(id) readRDS(paste0("data/processed/atsn/diem/", id, ".diem_sce.rds")) )
names(scel_at) <- lab_ids
seur_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
markers_at <- read.table("results/atsn/quantile/atsn.seur_markers.txt",
					  header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
markers_at <- markers_at[markers_at$p_val_adj < 0.05,]
clusters_at <- as.character(sort(as.numeric(unique(markers_at$cluster))))
keep_at <- sapply(clusters_at, function(x) length(grep("^mt-", markers_at[markers_at$cluster == x,"gene"], ignore.case=TRUE)) != 0)

thresh_at <- lapply(lab_ids, function(id) {
	quant <- readRDS(paste0("data/processed/atsn/quantile/", id, ".quantile_counts.rds"))
	return(quant$thresh)
})
names(thresh_at) <- lab_ids

# Adjust Seurat MT% for plotting
ch <- seur_at@meta.data$percent.mt > 10
seur_at@meta.data$percent.mt[ch] <- 10

#=========================================
# Get spliced reads fraction
#=========================================

props_rown <- c("Below", "Above")
prosp_coln <- c("Spliced", "Unspliced")

df_ad <- sce_ad@droplet_data[sce_ad@droplet_data$n_genes >= 200,]
props_ad <- table(df_ad$total_counts > thresh_ad, df_ad$SpliceFrctn < 50)
rownames(props_ad) <- props_rown
colnames(props_ad) <- prosp_coln
props_adm <- reshape2::melt(props_ad)
props_adm$DataSet <- "DiffPA"

df_mb <- sce_mb@droplet_data[sce_mb@droplet_data$n_genes >= 200,]
props_mb <- table(df_mb$total_counts > thresh_mb, df_mb$SpliceFrctn < 50)
rownames(props_mb) <- props_rown
colnames(props_mb) <- prosp_coln
props_mbm <- reshape2::melt(props_mb)
props_mbm$DataSet <- "Mouse Brain"

df_at <- lapply(1:6, function(i) {
                s <- scel_at[[i]]
                df_ats <- s@droplet_data[s@droplet_data$n_genes >= 200,] 
                ts <- table(df_ats$total_counts > thresh_at[i], df_ats$SpliceFrctn < 50)
        })

props_at <- df_at[[1]]
for (i in 2:6){
    props_at <- props_at + df_at[[i]]
}
rownames(props_at) <- props_rown
colnames(props_at) <- prosp_coln
props_atm <- reshape2::melt(props_at)
props_atm$DataSet <- "Adipose Tissue"

props_all <- do.call(rbind, list(props_atm, props_mbm, props_adm))
labell <- c("Background", "Nuclear"); names(labell) <- prosp_coln
props_all[,"Var2"] <- labell[props_all[,"Var2"]]
props_all$DataSet <- factor(props_all$DataSet, 
                            levels = c("DiffPA", "Mouse Brain", "Adipose Tissue"))

# Stacked bar plot

p_all <- ggplot(props_all, aes(x = Var1, y = value, fill = Var2)) + 
theme_classic() + 
ctl +
facet_wrap(~DataSet, scales = "free_y", ncol = 1) + 
geom_bar(position="stack", stat="identity") + 
labs(x = "Threshold", y = "Number of droplets")  + 
theme(axis.text.x = element_text(angle = 30, hjust = 1), 
      legend.title = element_blank(), 
      legend.position = "bottom")

#=========================================
# Call plot functions
#=========================================

# Barcode rank plot
pb1 <- barcode_rank_plot(sce_ad, title="DiffPA", ret=TRUE) + 
geom_hline(yintercept=thresh_ad, linetype="dashed", color = "red") + 
theme(plot.title=element_text(face="bold"))









p1 <- boxplot_mt(seur_ad) + ct + ggtitle("DiffPA")
#min_y <- layer_scales(p1)$y$range$range[1]
#max_y <- layer_scales(p1)$y$range$range[2]
#p1 <- p1 + ylim(min_y, max_y*yexp)
#for (i in clusters_ad[keep_ad]){
#	p1 <- p1 + annotate("text", x = i, y = max_y, label = "*", size = 12, fontface="bold", color="red", vjust=0.5)
#}

p1u <- plot_umap_meta(seur_ad, 
                      size = 0.8, alpha = 0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low="gray90", high="red3", name="MT%") +
ggtitle("DiffPA")

p1u_sf <- plot_umap_meta(seur_ad, 
                         colname = "SpliceFrctn", 
                      size = 0.8, alpha = 0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low="gray90", high="red3", name="% reads\nspliced") +
ggtitle("DiffPA")


p1bp_sf <- boxplot_mt(seur_ad, mt_col = "SpliceFrctn", yl = "% reads\nspliced") + 
    ct + ggtitle("DiffPA")











pb2 <- barcode_rank_plot(sce_mb, title="Mouse Brain", ret=TRUE) + 
geom_hline(yintercept=thresh_mb, linetype="dashed", color = "red") + 
theme(plot.title=element_text(face="bold"))

p2prs <- boxplot_mt(seur_mb,  mt_col="SpliceFrctn", ct_col="RNA_snn_res.0.8", yl="% reads spliced")


p2 <- boxplot_mt(seur_mb) + ct + ggtitle("Mouse Brain")
#min_y <- layer_scales(p2)$y$range$range[1]
#max_y <- layer_scales(p2)$y$range$range[2]
#p2 <- p2 + ylim(min_y, max_y*yexp)
#for (i in clusters_mb[keep_mb]){
#	p2 <- p2 + annotate("text", x = i, y = max_y, label = "*", size = 12, fontface="bold", color="red", vjust=0.5)
#}

p2u <- plot_umap_meta(seur_mb, 
                      size = 0.8, alpha = 0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low="gray90", high="red3", name="MT%") +
ggtitle("Mouse Brain")


p2u_sf <- plot_umap_meta(seur_mb, 
                         colname = "SpliceFrctn", 
                      size = 0.8, alpha = 0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low="gray90", high="red3", name="% reads\nspliced") +
ggtitle("Mouse Brain")

p2bp_sf <- boxplot_mt(seur_mb, mt_col = "SpliceFrctn", yl = "% reads\nspliced") + 
    ct + ggtitle("Mouse Brain")








pb3 <- lapply(lab_ids, function(x) {
	p <- barcode_rank_plot(scel_at[[x]], title=x, ret=TRUE)
	p <- p + geom_hline(yintercept=thresh_at[[x]], linetype="dashed", color = "red") + 
	theme(plot.title=element_text(face="bold"))
	return(p)
	})

p3prs <- boxplot_mt(seur_at,  mt_col="SpliceFrctn", ct_col="RNA_snn_res.0.8", yl="% reads spliced")


p3 <- boxplot_mt(seur_at) + ct + ggtitle("Adipose Tissue") + ylim(0,10)
#min_y <- layer_scales(p3)$y$range$range[1]
#max_y <- layer_scales(p3)$y$range$range[2]
#p3 <- p3 + ylim(min_y, max_y*yexp)
#for (i in clusters_at[keep_at]){
#	p3 <- p3 + annotate("text", x = i, y = max_y, label = "*", size = 12, fontface="bold", color="red", vjust=0.5)
#}

p3u <- plot_umap_meta(seur_at, 
                      size = 0.8, alpha=0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low="gray90", high="red3", name="MT%") + 
ggtitle("Adipose Tissue")  

p3u_sf <- plot_umap_meta(seur_at, 
                         colname = "SpliceFrctn", 
                      size = 0.8, alpha=0.8, 
                      order_col = TRUE) + 
ctl + 
scale_colour_gradient(low="gray90", high="red3", name="% reads\nspliced") + 
ggtitle("Adipose Tissue")  

p3bp_sf <- boxplot_mt(seur_at, mt_col = "SpliceFrctn", yl = "% reads\nspliced") + 
    ct + ggtitle("Adipose Tissue")




#=========================================
# Put all in 1 plot
#=========================================

dir_plot <- "results/plots/"
dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)

pb_all <- append(pb3, list("DiffPA"=pb1, "Mouse Brain"=pb2))
pdfname <- paste0(dir_plot, "quant.fig1.pdf")
jpgname <- paste0(dir_plot, "quant.fig1.jpeg")
pdf(pdfname, width=15,height=15)
ggarrange(ggarrange(ggarrange(plotlist=pb_all, ncol=2, nrow=4), 
                    p_all, 
                    ggarrange(p1u_sf, p2u_sf, p3u_sf, ncol=1), 
                    ncol=3, widths = c(0.35, .3, .35), 
                    labels=c("a", "b", "c"), font.label = list(size = 18, face="bold")),
          ggarrange(p1bp_sf, p2bp_sf, p3bp_sf, ncol = 3),
		  nrow=2, ncol=1, heights = c(.7, .3), labels=c("", "d"), 
          font.label = list(size = 18, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


pdfname <- paste0(dir_plot, "quant.MT.pdf")
jpgname <- paste0(dir_plot, "quant.MT.jpeg")
pdf(pdfname, width=13, height=10)
ggarrange(ggarrange(p1u, p2u, p3u, ncol = 3), 
          ggarrange(p1, p2, p3, ncol = 3), 
          labels=c("a", "b"), ncol = 1, font.label = list(size = 18, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))



