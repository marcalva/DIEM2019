
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
ct <- theme(text = element_text(size=16), 
            plot.title = element_text(size = 18, hjust = 0.5))

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


plot_overlap <- function(datf, col1 = "DIEM Cluster", col2 = "Seurat Cluster", title = ""){
    datfm <- reshape2::melt(as.matrix(datf * 100))
    datfm[,1] <- factor(datfm[,1])
    datfm[,2] <- factor(datfm[,2])
    p <- ggplot(datfm, aes_string(x = "Var2", y = "Var1")) + 
    geom_tile(aes_string(fill="value"), colour="white") + 
    scale_fill_viridis(name = "Percent\nOverlap", limits = c(0,100)) + 
    theme_minimal() + 
    h_theme + 
    ggtitle(title) + 
    scale_x_discrete(position = "top") + 
    ylab(col1) + xlab(col2) + ct
    return(p)
}


#=========================================
#=========================================

#=========================================
# Read in data
#=========================================

# DiffPA
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

all_names <- c("DiffPA", "Mouse Brain", lab_ids)
pcs_to_plot <- list(c("PC_1", "PC_2"), c("PC_3", "PC_4"), c("PC_5", "PC_6"))

set.seed(1)
#gplots <- sapply(1:length(seurl), function(i){
#                 s <- seurl[[i]]
#                 o <- overlaps[[i]]
#                 p1 <- lapply(pcs_to_plot[1:2], function(p) plot_pcs(s, x = p[1], y = p[2], legend = FALSE))
#                 p2 <- lapply(pcs_to_plot[3], function(p) plot_pcs(s, x = p[1], y = p[2]))
#                 p3 <- list(plot_overlap(o))
#                 return(c(p1, p2, p3))
#})
gplots <- sapply(1:length(seurl), function(i){
                 o <- overlaps[[i]]
                 list(plot_overlap(o, title = all_names[i]))})


# Arrange into plot
caps <- sapply(letters[1:length(seurl)], function(i) c(i, "", "", ""))
caps <- as.vector(caps)


summary(unlist(lapply(overlaps, function(o) apply(o, 2, max)*100)))




#=========================================
# Marker genes
#=========================================

marker_genes <- function(sce, only.pos =TRUE){
    zw <- colSums(sce@emo[[1]]$Z)
    zw <- zw[-1]
    pbar <- sce@emo[[1]]$params$Alpha
    pbar <- pbar[,-c(1)]
    pbar <- sweep(pbar, 2, colSums(pbar), "/")
    clusts <- 1:ncol(pbar)
    DE <- lapply(clusts, function(i){
                 zwo <- zw[-i] / sum(zw[-i])
                 pbarw <- pbar[,-c(i)] %*% zwo
                 diffs <- pbar[,i,drop=F] - pbarw[,drop=F]
                 if (only.pos){
                     diffs <- diffs[order(diffs[,1], decreasing=T),,drop=F]
                 } else {
                    diffs <- diffs[order(abs(diffs[,1]), decreasing=T),,drop=F]
                 }
                 return(diffs)
                 })
    return(DE)
}

sce_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.diem_sce.rds")
markers_ad <- marker_genes(sce_ad)
seur_diem_ad@meta.data$Cluster = factor(seur_diem_ad@meta.data$Cluster, 
                                        labels = 0:(length(unique(seur_diem_ad@meta.data$Cluster))-1))

sce_at <- readRDS("data/processed/atsn/diem/AT2.diem_sce.rds")
markers_at <- marker_genes(sce_at)
seur_diem_at2 <- seur_diem_at[,seur_diem_at$orig.ident == "AT2"]
seur_diem_at2@meta.data$Cluster = factor(seur_diem_at2@meta.data$Cluster, 
                                         labels = 0:(length(unique(seur_diem_at2@meta.data$Cluster))-1))
pumap1 <- plot_umap_labels(seur_diem_ad, "Cluster") + ggtitle("DiffPA DIEM clusters") + ct
pm1 <- plot_umap_gene(seur_diem_ad, "FN1", size = 1) + ggtitle("FN1") + ct
pm2 <- plot_umap_gene(seur_diem_ad, "CFD", size = 1) + ggtitle("CFD") + ct
pm3 <- plot_umap_gene(seur_diem_ad, "GPAM", size = 1) + ggtitle("GPAM") + ct

pumap2 <- plot_umap_labels(seur_diem_at2, "Cluster") + ggtitle("AT2 DIEM clusters") + ct
pn1 <- plot_umap_gene(seur_diem_at2, "VWF", size = 1) + ggtitle("VWF") + ct
pn2 <- plot_umap_gene(seur_diem_at2, "RBPJ", size = 1) + ggtitle("RBPJ") + ct
pn3 <- plot_umap_gene(seur_diem_at2, "GPAM", size = 1) + ggtitle("GPAM") + ct



fl <- list(size = 24, face="bold")

dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "mncluster.pdf")
jpgname <- paste0(dir_plot, "mncluster.jpeg")
pdf(pdfname, width=18, height=16)
ggarrange(
    ggarrange(pumap1, pm1, pm2, pm3, ncol=1), 
    ggarrange(pumap2, pn1, pn2, pn3, ncol=1),
    ggarrange(plotlist=gplots, ncol=2, nrow=4), 
    ncol = 3, labels = c("a", "b", "c"), font.label = fl, 
    widths = c(.25, .25, .5))
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))




