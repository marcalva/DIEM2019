
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

#=========================================
# Plot
#=========================================

paths_all <- c(adpath, mpath, at_path)
titles <- c("DiffPA", "Mouse Brain", at_labs)
names(paths_all) <- titles
sce_all <- lapply(paths_all, function(p){ readRDS(p) })
names(sce_all) <- c("DiffPA", "Mouse Brain", at_labs)
sce_pbmc <- readRDS(fresh68path)

ts <- lapply(names(sce_all), function(s){
             datf <- sce_all[[s]]@test_data
             datf[,"Dataset"] <- s
             datf <- datf[order(datf$Call),]
             return(datf)
                               })

ts_pbmc <- sce_pbmc@test_data

ts_all <- as.data.frame(do.call(rbind, ts))
ts_all[,"Dataset"] <- factor(ts_all$Dataset, levels=names(sce_all))
alpha <- 0.1

p <- ggplot(ts_all, aes(y=score.debris, x=total_counts)) + 
geom_point(alpha=alpha, size = 0.5) +
facet_wrap(~Dataset, nrow=4, scales="free") + 
ylab("Debris score") +
xlab("Total UMI counts") +
scale_x_continuous(trans='log10', labels=scales::comma) +
theme_minimal() + 
theme(text=element_text(size=12), 
      axis.text.x = element_text(angle = 30, hjust = 0.9))

dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "n_umi.debris_score.pdf")
jpgname <- paste0(dir_plot, "n_umi.debris_score.jpeg")
ggsave(pdfname, width = 4, height = 7)
ggsave(jpgname, width = 4, height = 7)


p <- ggplot(ts_pbmc, aes(y=score.debris, x=total_counts)) + 
geom_point(alpha=alpha, size = 1.2) +
ylab("Debris score") +
xlab("Total UMI counts") +
ggtitle("Fresh 68K PBMCs") + 
scale_x_continuous(trans='log10', labels=scales::comma) +
theme_minimal() + 
theme(text=element_text(size=16), 
      axis.text.x = element_text(angle = 30, hjust = 0.9))

dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "n_umi.debris_score.pbmc.pdf")
jpgname <- paste0(dir_plot, "n_umi.debris_score.pbmc.jpeg")
ggsave(pdfname, width = 6, height = 7)
ggsave(jpgname, width = 6, height = 7)


