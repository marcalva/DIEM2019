
setwd("../")

library(diem)
suppressMessages(library(DropletUtils))
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
source("scripts/common/plotting.R")


# Common function and themes
# Boxplot of MT%
boxplot_mt <- function(x, mt_col="percent.mt", ct_col="RNA_snn_res.0.8"){
	require(reshape2)
	df <- data.frame(Mito=x@meta.data[,mt_col], Cluster=x@meta.data[,ct_col])
	df <- reshape2::melt(df)
	p <- ggplot(df, aes(x=Cluster, y=value, fill=Cluster)) +
	geom_boxplot() + theme_classic() + 
	ylab("MT%")
	return(p)
}

ct <- theme(legend.position="none",
			text=element_text(size=16),
			plot.title=element_text(size=18, hjust=0.5, face="bold")
			)

yexp <- 1.1


#=========================================
# Read in data
#=========================================

# Adipocyte

labl <- "adpcyte"
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

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
dir_plot <- paste0("results/", labl, "/plots/")

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
dir_plot <- paste0("results/", labl, "/plots/")

scel_at <- lapply(lab_ids, function(id) readRDS(paste0("data/processed/atsn/diem/", id, ".diem_sce.rds")) )
names(scel_at) <- lab_ids
seur_at <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
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

#=========================================
# Create plot-ready data frames
#=========================================


df_ad <- sce_ad@droplet_data
df_ad <- df_ad[df_ad$total_counts > 100,]
df_ad$inflection <- thresh_ad

df_ad$name <- rep("Below threshold", nrow(df_ad))
df_ad[df_ad$total_counts > thresh_ad, "name"] <- "Above threshold"



df_mb <- sce_mb@droplet_data
df_mb <- df_mb[df_mb$total_counts > 100,]
df_mb$inflection <- thresh_mb

df_mb$name <- rep("Below threshold", nrow(df_mb))
df_mb[df_mb$total_counts > thresh_mb, "name"] <- "Above threshold"



dfl <- lapply(1:length(scel_at), function(i) {
	sce <- scel_at[[i]]
	df <- sce@droplet_data
	df <- df[df$total_counts > 100,]
	df$Sample <- sce@name
	df$inflection <- thresh_at[[sce@name]]
	df$name <- rep("Below threshold", nrow(df))
	df[df$total_counts > thresh_at[[i]], "name"] <- "Above threshold"
	return(df)
})

df_at <- do.call(rbind, dfl)


#=========================================
# Create plots
#=========================================

# Adipocyte UMI-MT% plot
pu1 <- ggplot(df_ad, aes(x = total_counts, y = pct.mt)) + 
geom_point(alpha=0.1) + 
theme_minimal() + 
theme(text=element_text(size=16), 
	plot.title=element_text(hjust=0.5, face="bold"), 
	axis.text.x=element_text(angle=90, hjust=1), 
	strip.text=element_text(face="bold")) + 
xlab("Total UMI count") + 
ylab("MT%") + 
ggtitle("DiffPA") + 
ylim(c(0, 100)) + 
scale_x_continuous(trans='log10', labels = comma) + 
geom_vline(aes(xintercept = inflection), col="red", linetype="dashed")

# Mouse brain UMI-MT% plot
pu2 <- ggplot(df_mb, aes(x = total_counts, y = pct.mt)) + 
geom_point(alpha=0.1) + 
theme_minimal() + 
theme(text=element_text(size=16), 
	plot.title=element_text(hjust=0.5, face="bold"), 
	axis.text.x=element_text(angle=90, hjust=1), 
	strip.text=element_text(face="bold")) + 
xlab("Total UMI count") + 
ylab("MT%") + 
ggtitle("Mouse Brain") + 
ylim(c(0, 100)) + 
scale_x_continuous(trans='log10', labels = comma) + 
geom_vline(aes(xintercept = inflection), col="red", linetype="dashed")

# AT UMI-MT% plot
pu3 <- lapply(1:length(dfl), function(i){
	df <- dfl[[i]]
	p <- ggplot(df, aes(x = total_counts, y = pct.mt)) + 
	geom_point(alpha=0.1) + 
	theme_minimal() + 
	theme(text=element_text(size=16), 
		plot.title=element_text(hjust=0.5, face="bold"), 
		axis.text.x=element_text(angle=90, hjust=1), 
		strip.text=element_text(face="bold")) + 
	xlab("Total UMI count") + 
	ylab("MT%") + 
	ggtitle(lab_ids[[i]]) + 
	ylim(c(0, 10)) + 
	scale_x_continuous(trans='log10', labels = comma) + 
	geom_vline(aes(xintercept = inflection), col="red", linetype="dashed")
	return(p)
})


# Adipocyte MT% boxplot
pm1 <- ggplot(df_ad, aes(x = name, y = pct.mt)) + 
geom_boxplot() + 
theme_minimal() + 
theme(axis.title.x=element_blank(), text=element_text(size=18),
	plot.title=element_text(size=22, hjust=0.5, face="bold")) + 
ylab("MT%") + 
ggtitle("DiffPA")

# Mouse brain MT% boxplot
pm2 <- ggplot(df_mb, aes(x = name, y = pct.mt)) + 
geom_boxplot() + 
theme_minimal() + 
theme(axis.title.x=element_blank(), text=element_text(size=18),
	plot.title=element_text(size=22, hjust=0.5, face="bold")) + 
ylab("MT%") + 
ggtitle("Mouse Brain")

# AT MT% boxplot
pm3 <- lapply(1:length(dfl), function(i){
	df <- dfl[[i]]
	ggplot(df, aes(x = name, y = pct.mt)) + 
	geom_boxplot() + 
	theme_minimal() + 
	theme(axis.title.x=element_blank(), text=element_text(size=18),
		plot.title=element_text(size=22, hjust=0.5, face="bold")) + 
	ylab("MT%") + 
	ggtitle(lab_ids[[i]])
})



#=========================================
# Put all in 1 plot
#=========================================

dir_plot <- "results/plots/"; dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)

# For now just include correlation scatter plot
pdfname <- paste0(dir_plot, "mt_umi.pdf")
jpgname <- paste0(dir_plot, "mt_umi.jpeg")
pdf(pdfname, width=12,height=12)
ggarrange(ggarrange(pu1, pu2, nrow=1),
		  ggarrange(plotlist=pu3, nrow=2, ncol=3), 
		  nrow=2, ncol=1, heights=c(2,3), font.label = list(size = 18, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


#=========================================
# Print out quantifications
#=========================================

# Not part of the paper

get_pct_above <- function(df, sd_above=0){
	mt_mean_below <- mean(df[df$name == "Below threshold","pct.mt"])
	mt_sd_below <- sd(df[df$name == "Below threshold","pct.mt"])
	mt_mean_above <- mean(df[df$name == "Above threshold","pct.mt"])
	table_above <- table(df[df$name == "Above threshold","pct.mt"] > (mt_mean_below+(sd_above*mt_sd_below)))
	return((table_above/sum(table_above))[["TRUE"]])
}


# DiffPA
pct_above_ad <- get_pct_above(df_ad)
print(pct_above_ad)
# 0.6324153

# Mouse brain
pct_above_mb <- get_pct_above(df_mb)
print(pct_above_mb)
# 0.1625616

# AT
pct_above_l <- lapply(dfl, get_pct_above)
print(pct_above_l)
# 0.05770224

