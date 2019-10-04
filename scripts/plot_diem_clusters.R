
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
library(viridis)
library(reshape2)
source("scripts/common/standard_seurat.R")
source("scripts/common/plotting.R")
source("scripts/common/overlap_graph.R")


#=========================================
# Functions
#=========================================
ct <- theme(text = element_text(size=16))

plot_pcs <- function(seur, x = "PC_1", y = "PC_2", color="Cluster", legend=TRUE){
    datf <- cbind(seur@meta.data, seur@reductions$pca@cell.embeddings)
    datf <- datf[sample(rownames(datf)),,drop=FALSE]
    p <- ggplot(datf, aes_string(x = x, y = y, color=color)) + 
    geom_point(size = 0.5) + 
    theme_minimal() + 
    scale_color_discrete(name="DIEM\nCluster") + 
    ct
    if (!legend){
        p <- p + theme(legend.position = "none")
    }
    return(p)
}

p_blank <- ggplot() + theme_void()

h_theme <- theme(panel.grid=element_blank(),
                 panel.border = element_blank(),
                 text=element_text(size=16))


plot_overlap <- function(datf, col1 = "DIEM Cluster", col2 = "Seurat Cluster"){
    datfm <- reshape2::melt(as.matrix(datf * 100))
    datfm[,1] <- factor(datfm[,1])
    datfm[,2] <- factor(datfm[,2])
    p <- ggplot(datfm, aes_string(x = "Var2", y = "Var1")) + 
    geom_tile(aes_string(fill="value"), colour="white") + 
    scale_fill_viridis(name = "Percent\nOverlap", limits = c(0,100)) + 
    theme_minimal() + 
    h_theme + 
    scale_x_discrete(position = "top") + 
    ylab(col1) + xlab(col2)
    return(p)
}


#=========================================
#=========================================

#=========================================
# Read in data
#=========================================

# Adipocyte
seur_diem_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.seur_obj.rds")

# Mouse brain
seur_diem_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.seur_obj.rds")

# Adipose tissue
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
seur_diem_at <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_at <- lapply(lab_ids, function(i){
                  cells_all <- rownames(seur_diem_at@meta.data)
                  keep <- seur_diem_at$orig.ident == i
                  cells <- cells_all[keep]
                  s <- subset(seur_diem_at, cells=cells)
                  s <- seurat_norm(s, scale.factor=1000)
                  s <- RunPCA(s, npcs = 30, verbose = TRUE)
})

seurl <- c(seur_diem_ad, seur_diem_mb, seur_at)
names(seurl) <- c("DiffPA", "Mouse Brain", lab_ids)

# Drop clusters
seurl <- lapply(seurl, function(s){
                s$Cluster <- as.integer(s$Cluster) - 2
                s$Cluster <- factor(s$Cluster)
                return(s)
})

#=========================================
# Cluster
#=========================================

seurl <- lapply(seurl, function(s){
                s <- FindNeighbors(s, dims = 1:30)
                s <- FindClusters(s)
                return(s)
})

#=========================================
# Overlap Seurat clusters and MN clusters
#=========================================

overlaps <- lapply(seurl, function(s){
                   tb <- table(s$Cluster,  s$RNA_snn_res.0.8)
                   tb <- as.data.frame.matrix(tb)
                   pct <- sweep(tb, 2, colSums(tb), "/")
                   return(pct)
})
                   

#=========================================
# Plot
#=========================================


pcs_to_plot <- list(c("PC_1", "PC_2"), c("PC_3", "PC_4"), c("PC_5", "PC_6"))

set.seed(1)
gplots <- sapply(1:length(seurl), function(i){
                 s <- seurl[[i]]
                 o <- overlaps[[i]]
                 p1 <- lapply(pcs_to_plot[1:2], function(p) plot_pcs(s, x = p[1], y = p[2], legend = FALSE))
                 p2 <- lapply(pcs_to_plot[3], function(p) plot_pcs(s, x = p[1], y = p[2]))
                 p3 <- list(plot_overlap(o))
                 return(c(p1, p2, p3))
})


# Arrange into plot
caps <- sapply(letters[1:length(seurl)], function(i) c(i, "", "", ""))
caps <- as.vector(caps)

dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "mncluster_pcs.pdf")
jpgname <- paste0(dir_plot, "mncluster_pcs.jpeg")
pdf(pdfname, width=20, height=24)
ggarrange(plotlist=gplots, ncol=4, nrow=8, widths=c(3,3,3.5,3.5), 
          labels=caps, 
          font.label = list(size = 24, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))

summary(unlist(lapply(overlaps, function(o) apply(o, 2, max)*100)))

