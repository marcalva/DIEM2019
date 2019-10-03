
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

get_n_s_mt <- function(markers){
	markers <- markers[markers$p_val_adj < 0.05,]
	nmt <- length(grep("^mt-", markers[,"gene"], ignore.case=TRUE))
	return(nmt)
}


#=========================================
#=========================================

#=========================================
# Read in data
#=========================================


# Adipocyte
labl <- "adpcyte"
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

# Read in DIEM SCE
sce_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.diem_sce.rds")
seur_diem_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.seur_obj.rds")
seur_quant_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
seur_ED_ad <- readRDS("data/processed/adpcyte/emptydrops/adpcyte.seur_obj.rds")
markers_diem_ad <- read.table("results/adpcyte/diem/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_quant_ad <- read.table("results/adpcyte/quantile/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_ED_ad <- read.table("results/adpcyte/emptydrops/adpcyte.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)

# Mouse brain
labl <- "mouse_nuclei_2k"
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

sce_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.diem_sce.rds")
seur_diem_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.seur_obj.rds")
seur_quant_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
seur_ED_mb <- readRDS("data/processed/mouse_nuclei_2k/emptydrops/mouse_nuclei_2k.seur_obj.rds")
markers_diem_mb <- read.table("results/mouse_nuclei_2k/diem/mouse_nuclei_2k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_quant_mb <- read.table("results/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_ED_mb <- read.table("results/mouse_nuclei_2k/emptydrops/mouse_nuclei_2k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)


# Adipose tissue
labl <- "atsn"
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

scel <- lapply(lab_ids, function(id) readRDS(paste0("data/processed/atsn/diem/", id, ".diem_sce.rds")))
names(scel) <- lab_ids
seur_diem_at <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_quant_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
seur_ED_at <- readRDS("data/processed/atsn/emptydrops/atsn.seur_obj.rds")
markers_diem_at <- read.table("results/atsn/diem/atsn.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_quant_at <- read.table("results/atsn/quantile/atsn.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_ED_at <- read.table("results/atsn/emptydrops/atsn.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)


#=========================================
# Main Figure
#=========================================

#=========================================
# Other stuff
#=========================================


# Get number of nuclei
methd_names <- c("Quantile", "EmptyDrops", "DIEM")
methd <- factor(c(methd_names, methd_names, methd_names), levels=c("Quantile", "EmptyDrops", "DIEM"))
exper <- factor(c("DiffPA", "DiffPA", "DiffPA", 
		   "Mouse\nBrain", "Mouse\nBrain", "Mouse\nBrain", 
		   "Adipose\nTissue", "Adipose\nTissue", "Adipose\nTissue"),
                levels=c("DiffPA", "Mouse\nBrain", "Adipose\nTissue"))

# Get number of MT-enoded marker genes
nmt <- c(get_n_s_mt(markers_quant_ad), get_n_s_mt(markers_ED_ad), get_n_s_mt(markers_diem_ad),
		 get_n_s_mt(markers_quant_mb), get_n_s_mt(markers_ED_mb), get_n_s_mt(markers_diem_mb),
		 get_n_s_mt(markers_quant_at), get_n_s_mt(markers_ED_at), get_n_s_mt(markers_diem_at))
dat_mtn <- data.frame(Method=methd, Experiment=exper, MT=nmt)

# Barplot of total MT markers
p <- ggplot(dat_mtn, aes(x=Method, y=MT, fill=Method)) + 
geom_bar(stat="identity", color="black", position=position_dodge()) + 
facet_wrap(~Experiment, scale="free_y") + 
theme_minimal() + ylab("Total number of\nMT markers") + xlab("") + 
theme(axis.text.x=element_blank(),
      axis.title.x=element_blank(),
      text=element_text(size=18),
      plot.title=element_text(hjust=0.5), 
      panel.grid.major.x=element_blank())


#=========================================
# Supp
#=========================================

dir_plot <- paste0("results/plots/")
pdfname <- paste0(dir_plot, "barplot_mt.pdf")
jpgname <- paste0(dir_plot, "barplot_mt.jpeg")
pdf(pdfname, width=8, height=6)
p
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))

