
# Plot values of splice fraction in boxplots

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
# Get splice fractions
#=========================================

at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
at_path <- paste0("data/processed/atsn/diem/", at_labs, ".diem_sce.rds")
mpath <- "data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.diem_sce.rds"
adpath <- "data/processed/adpcyte/diem/adpcyte.diem_sce.rds"
atintpath <- "data/processed/atsn/diem_integrated.fltr.1/atsn.diem_sce.rds"
fresh68path<- "data/processed/fresh_68k/diem/fresh_68k.diem_sce.rds"

paths_all <- c(adpath, mpath, at_path, atintpath, fresh68path)
ids_all <- c("adpcyte", "mouse_nuclei_2k", at_labs, "AT_intg", "PBMC68K")
titles <- c("DiffPA", "Mouse Brain", at_labs, "AT integrated", "Fresh 68K PBMCs")

sce_all <- lapply(paths_all, function(p){ readRDS(p) })

plots <- lapply(1:(length(sce_all)), function(i){
    s <- sce_all[[i]]
    p <- plot_dist(s, ret=TRUE) + ggtitle(titles[i]) + 
    theme(axis.text.x = element_text(angle=45, hjust = 1)) 
    return(p)
})

dir_plot <- "results/plots/"; dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)

fl <- list(size = 30, face="bold")

pdfname <- paste0(dir_plot, "init_clust.filters.pdf")
jpgname <- paste0(dir_plot, "init_clust.filters.jpeg")
pdf(pdfname, width=12,height=18)
ggarrange(ggarrange(plotlist = plots[1:2], nrow = 1, 
                    labels=c("a", "b"), font.label = fl), 
          ggarrange(plotlist = plots[3:8], ncol = 3, nrow = 2), 
          ggarrange(plotlist = plots[9:10], nrow = 1, ncol = 2,
                    labels=c("d", "e"), font.label = fl), 
          nrow = 3, labels=c("", "c", ""), font.label = fl, heights = c(.3, .4, .3))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))

