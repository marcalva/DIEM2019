
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

plot_umi_genes <- function(x,
                           color="Call",
                           color_name="Call",
                           alpha=0.1,
                           ret=FALSE){

    df <- x@droplet_data[x@test_set,]

    p <- ggplot(df, aes(x=total_counts, y=n_genes)) + geom_point(alpha=alpha, aes_string(colour=color)) +
    xlab("UMI Counts") +
    ylab("Genes Detected") +
    scale_x_continuous(trans='log10', breaks=set_breaks_10, labels=scales::comma) + 
    scale_y_log10() +
    theme_minimal() + theme(text=element_text(size=22)) +
    scale_colour_discrete(name=color_name)
    if (ret) return(p)
    else print(p)
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

# Mouse brain
labl <- "mouse_nuclei_2k"
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

sce_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.diem_sce.rds")

# Adipose tissue
labl <- "atsn"
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

scel <- lapply(lab_ids, function(id) readRDS(paste0("data/processed/atsn/diem/", id, ".diem_sce.rds")))
names(scel) <- lab_ids

#=========================================
# Plot
#=========================================
sce_all <- rev(c(DiffPA=sce_ad, "Mouse Brain"=sce_mb, scel))

ts <- lapply(names(sce_all), function(s){
             datf <- sce_all[[s]]@droplet_data[sce_all[[s]]@test_set,]
             datf[,"Dataset"] <- s
             datf[datf$Call == "Removed","Call"] <- "Debris"
             datf <- datf[order(datf$Call),]
             return(datf)
                               })

ts_all <- do.call(rbind, ts)
ts_all$Call <- factor(ts_all$Call, levels=c("Debris", "Clean"))
ts_all$ProbDebris <- 1 - ts_all$CleanProb
ts_all$Dataset <- factor(ts_all$Dataset, levels=rev(names(sce_all)))
alpha <- 0.1; color_name <- color <-"ProbDebris"

p <- ggplot(ts_all, aes(x=total_counts, y=n_genes)) + 
geom_point(alpha=alpha, aes_string(colour=color)) +
facet_wrap(~Dataset, nrow=4, scales="free") + 
xlab("Total UMI count") +
ylab("Number of genes detected") +
scale_x_continuous(trans='log10', labels=scales::comma) +
scale_y_continuous(trans='log10', labels=scales::comma) +
theme_minimal() + theme(text=element_text(size=22)) +
scale_color_distiller(palette="RdBu", direction=-1, name="Probability\ndebris")
# scale_colour_discrete(name=color_name)

dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "n_umi.n_gene.call.diem.pdf")
jpgname <- paste0(dir_plot, "n_umi.n_gene.call.diem.jpeg")
pdf(pdfname, width=10, height=12)
p
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))


