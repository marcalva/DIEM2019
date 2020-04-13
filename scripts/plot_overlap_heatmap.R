
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
source("scripts/common/plotting.R")
source("scripts/common/overlap_graph.R")


#=========================================
# Functions
#=========================================

ct <- theme(text = element_text(size = 20), 
            axis.title = element_text(size = 22), 
            plot.title = element_text(size = 26, hjust = 0.5)) 

plot_overlap <- function(datf, col1 = "DIEM Cluster", col2 = "Seurat Cluster", title = ""){
    require(viridis)
    datfm <- reshape2::melt(as.matrix(datf * 100))
    datfm[,1] <- factor(datfm[,1])
    datfm[,2] <- factor(datfm[,2])
    p <- ggplot(datfm, aes_string(x = "Var2", y = "Var1")) + 
    geom_tile(aes_string(fill="value"), colour="white") + 
    scale_fill_viridis(name = "Percent\nOverlap", limits = c(0,100)) + 
    theme_minimal() + 
    ggtitle(title) + 
    scale_x_discrete(position = "top") + 
    ylab(col1) + 
    xlab(col2) + 
    ct
    return(p)
}

# Compare clusters and return data frame for plotting
comp_clust <- function(seur1, seur2, ct_col1="seurat_clusters", ct_col2="seurat_clusters"){
	clusters1 <- as.character(seur1@meta.data[,ct_col1])
	clusters2 <- as.character(seur2@meta.data[,ct_col2])
	ct1 <- as.character(sort(unique(clusters1)))
	ct2 <- as.character(sort(unique(clusters2)))

	dfc <- as.data.frame(matrix(nrow=length(ct1), ncol=length(ct2)))
	rownames(dfc) <- ct1
	colnames(dfc) <- ct2
	for (c1 in ct1){
		md1 <- seur1@meta.data[clusters1 == c1,]
		names1 <- rownames(md1)
		n_all <- length(names1)
		for (c2 in ct2){
			md2 <- seur2@meta.data[clusters2 == c2,]
			names2 <- rownames(md2)
			n_intr <- length(intersect(names1, names2))
			dfc[c1, c2] <- n_intr/n_all
		}
	}
	return(dfc)
}

#=========================================
#=========================================

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

# fresh 68k
seur_diem_b <- readRDS("data/processed/fresh_68k/diem/fresh_68k.seur_obj.rds")
seur_ED_b <- readRDS("data/processed/fresh_68k/emptydrops/fresh_68k.seur_obj.rds")

seur_ad <- list("DIEM" = seur_diem_ad, "EmptyDrops" = seur_ED_ad, "Quantile" = seur_quant_ad)
seur_mb <- list("DIEM" = seur_diem_mb, "EmptyDrops" = seur_ED_mb, "Quantile" = seur_quant_mb)
seur_at <- list("DIEM" = seur_diem_at, "EmptyDrops" = seur_ED_at, "Quantile" = seur_quant_at)
seur_b <- list("DIEM" = seur_diem_b, "EmptyDrops" = seur_ED_b)

#=========================================
# Plot
#=========================================


dir_plot <- "results/plots/overlaps/"
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

w <- 6
h <- 5


# DiffPA
ed_v_diem <- comp_clust(seur_ED_ad, seur_diem_ad, ct_col1="RNA_snn_res.0.8", ct_col2="RNA_snn_res.0.8")

q_v_diem <- comp_clust(seur_quant_ad, seur_diem_ad, ct_col1="RNA_snn_res.0.8", ct_col2="RNA_snn_res.0.8")

q_v_ed <- comp_clust(seur_quant_ad, seur_ED_ad, ct_col1="RNA_snn_res.0.8", ct_col2="RNA_snn_res.0.8")

p <- plot_overlap(ed_v_diem, "EmptyDrops", "DIEM")
outfn <- paste0(dir_plot, "overlap.pass.DiffPA.ed_v_diem.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.DiffPA.ed_v_diem.pdf")
ggsave(outfn, width = w, height = h)

p <- plot_overlap(q_v_diem, "Quantile", "DIEM")
outfn <- paste0(dir_plot, "overlap.pass.DiffPA.q_v_diem.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.DiffPA.q_v_diem.pdf")
ggsave(outfn, width = w, height = h)

p <- plot_overlap(q_v_ed, "Quantile", "EmptyDrops")
outfn <- paste0(dir_plot, "overlap.pass.DiffPA.q_v_ed.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.DiffPA.q_v_ed.pdf")
ggsave(outfn, width = w, height = h)


# Mouse brain
ed_v_diem <- comp_clust(seur_ED_mb, seur_diem_mb, ct_col1="RNA_snn_res.0.8", ct_col2="RNA_snn_res.0.8")

q_v_diem <- comp_clust(seur_quant_mb, seur_diem_mb, ct_col1="RNA_snn_res.0.8", ct_col2="RNA_snn_res.0.8")

q_v_ed <- comp_clust(seur_quant_mb, seur_ED_mb, ct_col1="RNA_snn_res.0.8", ct_col2="RNA_snn_res.0.8")

p <- plot_overlap(ed_v_diem, "EmptyDrops", "DIEM")
outfn <- paste0(dir_plot, "overlap.pass.MouseBrain.ed_v_diem.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.MouseBrain.ed_v_diem.pdf")
ggsave(outfn, width = w, height = h)

p <- plot_overlap(q_v_diem, "Quantile", "DIEM")
outfn <- paste0(dir_plot, "overlap.pass.MouseBrain.q_v_diem.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.MouseBrain.q_v_diem.pdf")
ggsave(outfn, width = w, height = h)

p <- plot_overlap(q_v_ed, "Quantile", "EmptyDrops")
outfn <- paste0(dir_plot, "overlap.pass.MouseBrain.q_v_ed.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.MouseBrain.q_v_ed.pdf")
ggsave(outfn, width = w, height = h)


# Adipose tissue
ed_v_diem <- comp_clust(seur_ED_at, seur_diem_at, ct_col1="CellType", ct_col2="CellType")

q_v_diem <- comp_clust(seur_quant_at, seur_diem_at, ct_col1="CellType", ct_col2="CellType")

q_v_ed <- comp_clust(seur_quant_at, seur_ED_at, ct_col1="CellType", ct_col2="CellType")

p <- plot_overlap(ed_v_diem, "EmptyDrops", "DIEM") + 
theme(text = element_text(size = 16), 
      axis.text.x = element_text(angle = 90, hjust = 0), 
      axis.text.x.top = element_text(vjust = 0.5))
outfn <- paste0(dir_plot, "overlap.pass.atsn.ed_v_diem.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.atsn.ed_v_diem.pdf")
ggsave(outfn, width = w, height = h)

p <- plot_overlap(q_v_diem, "Quantile", "DIEM") + 
theme(text = element_text(size = 16),
      axis.text.x = element_text(angle = 90, hjust = 0), 
      axis.text.x.top = element_text(vjust = 0.5))
outfn <- paste0(dir_plot, "overlap.pass.atsn.q_v_diem.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.atsn.q_v_diem.pdf")
ggsave(outfn, width = w, height = h)

p <- plot_overlap(q_v_ed, "Quantile", "EmptyDrops") + 
theme(text = element_text(size = 16),
      axis.text.x = element_text(angle = 90, hjust = 0), 
      axis.text.x.top = element_text(vjust = 0.5))
outfn <- paste0(dir_plot, "overlap.pass.atsn.q_v_ed.jpeg")
ggsave(outfn, width = w, height = h)
outfn <- paste0(dir_plot, "overlap.pass.atsn.q_v_ed.pdf")
ggsave(outfn, width = w, height = h)

