
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

# Compare clusters and return data frame for plotting
comp_clust <- function(seur1, seur2, ct_col1="seurat_clusters", ct_col2="seurat_clusters"){
	clusters1 <- as.character(seur1@meta.data[,ct_col1])
	clusters2 <- as.character(seur2@meta.data[,ct_col2])

	ct1 <- as.character(sort(as.numeric(unique(clusters1))))
	ct2 <- as.character(sort(as.numeric(unique(clusters2))))

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

p_blank <- ggplot() + theme_void()

#=========================================
#=========================================

#=========================================
# Read in data
#=========================================

# fresh 68K
labl <- "fresh_68k"
dp <- paste0("data/processed/", labl, "/")
dir_plot <- paste0("results/", labl, "/plots/")

# Read in DIEM SCE
sce_ad <- readRDS("data/processed/fresh_68k/diem/fresh_68k.diem_sce.rds")
seur_diem <- readRDS("data/processed/fresh_68k/diem/fresh_68k.seur_obj.rds")
seur_quant <- readRDS("data/processed/fresh_68k/quantile/fresh_68k.seur_obj.rds")
seur_ED <- readRDS("data/processed/fresh_68k/emptydrops/fresh_68k.seur_obj.rds")

hb_genes <- c("HBA1", "HBA2", "HBB")
seur_diem <- PercentageFeatureSet(seur_diem, features=hb_genes, col.name="Hemoglobin")
seur_quant <- PercentageFeatureSet(seur_quant, features=hb_genes, col.name="Hemoglobin")
seur_ED <- PercentageFeatureSet(seur_ED, features=hb_genes, col.name="Hemoglobin")


markers_diem <- read.table("results/fresh_68k/diem/fresh_68k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_quant <- read.table("results/fresh_68k/quantile/fresh_68k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_ED <- read.table("results/fresh_68k/emptydrops/fresh_68k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)

#=========================================
# Plot
#=========================================

# Compare the characteristics of droplets removed by DIEM to those not
names_diff <- setdiff(colnames(seur_ED), colnames(seur_diem))
df_kept <- data.frame(seur_ED@meta.data[colnames(seur_diem),], Filter="DIEM &\nEmptyDrops")
df_rm <- data.frame(seur_ED@meta.data[names_diff,], Filter="EmptyDrops\nOnly")

datn <- data.frame(Method=c("DIEM &\nEmptyDrops", "EmptyDrops\nOnly"), 
                   y=c(18, 18), 
                   N=c(paste0("n=", as.character(nrow(df_kept))), paste0("n=", as.character(nrow(df_rm)))))

datf_kept_rm <- rbind(df_kept, df_rm)

datfm <- reshape2::melt(datf_kept_rm[,c("Filter", "percent.mt", "MALAT1")])

relabel <- c("percent.mt"="Mitochondria", "MALAT1"="MALAT1")

prmall <- ggplot(datfm, aes(x=Filter, y=value)) + 
geom_boxplot(outlier.shape = NA) + theme_bw() +
facet_wrap(~variable, labeller = labeller(variable=relabel)) + 
geom_text(data=datn, aes(x=Method, y=y, label=N), size=6) + 
ylim(0,20) + ylab("Percent") + 
theme(text=element_text(size=16),
      axis.title.x=element_blank(),
      panel.grid.major.x=element_blank(), 
      plot.title=element_text(hjust=0.5))




# Get number of nuclei
methd_names <- factor(c("Quantile", "EmptyDrops", "DIEM"), levels=c("Quantile", "EmptyDrops", "DIEM"))
ncel <- c(ncol(seur_quant), ncol(seur_ED), ncol(seur_diem))
datf <- data.frame(Method=methd_names, Nuclei=ncel)

pb1 <- ggplot(datf, aes(x=Method, y=Nuclei, fill=Method)) +
#scale_x_discrete(limits=c("DiffPA", "Mouse\nBrain", "Adipose\nTissue" )) + 
geom_bar(stat="identity", color="black", position=position_dodge()) +
geom_text(aes(label=Nuclei), size=8, position=position_dodge(width=0.9), vjust=-0.25) + 
theme_minimal() + ylab("Total number\nof Nuclei") + ylim(0,80e3) + 
ggtitle("Fresh 68K PBMCs") + 
theme(legend.position="none",
      axis.title.x=element_blank(),
      text=element_text(size=16),
      panel.grid.major.x=element_blank(), 
      plot.title=element_text(hjust=0.5))

# Cluster overlap between EmptyDrops and DIEM
ctype="RNA_snn_res.0.8"

overlapq <- comp_clust(seur_quant, seur_diem, ct_col1=ctype, ct_col2=ctype)
mt_quant <- sapply(rownames(overlapq), function(x){mean(seur_quant@meta.data[seur_quant@meta.data[,ctype] == x,"percent.mt"])})
mt_diemq <- sapply(colnames(overlapq), function(x){mean(seur_diem@meta.data[seur_diem@meta.data[,ctype] == x,"percent.mt"])})

overlape <- comp_clust(seur_ED, seur_diem, ct_col1=ctype, ct_col2=ctype)
mt_ED <- sapply(rownames(overlape), function(x){mean(seur_ED@meta.data[seur_ED@meta.data[,ctype] == x,"percent.mt"])})
mt_dieme <- sapply(colnames(overlape), function(x){mean(seur_diem@meta.data[seur_diem@meta.data[,ctype] == x,"percent.mt"])})

pe1 <- overlap_graph(overlapq, mt_quant, mt_diemq, labels1="Quantile", labels2="DIEM", 
                     main="Fresh 68K PBMCs", col_title="Total UMIs", 
                     scale_size = .8)
pe2 <- overlap_graph(overlape, mt_ED, mt_dieme, labels1="EmptyDrops", labels2="DIEM", 
                     main="Fresh 68K PBMCs", col_title="Total UMIs", 
                     scale_size = .8)

#
# UMAP
#

pu1 <- plot_umap_labels(seur_quant, legend_title="\nClusters") + 
ggtitle("Quantile") + 
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(size=22, hjust=0.5, face="bold"))

pu2 <- plot_umap_labels(seur_ED, legend_title="\nClusters") + 
ggtitle("EmptyDrops") +             
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(size=22, hjust=0.5, face="bold"))

pu3 <- plot_umap_labels(seur_diem, legend_title="\nClusters") + 
ggtitle("DIEM") +             
theme(legend.position="none", text=element_text(size=18), plot.title=element_text(size=22, hjust=0.5, face="bold"))

# Arrange in plot
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "fresh_68k.pdf")
jpgname <- paste0(dir_plot, "fresh_68k.jpeg")
pdf(pdfname, width=20, height=15)
ggarrange(ggarrange(prmall, pu1, ncol=1, nrow=2, labels=c("a", "b"), font.label = list(size = 26, face="bold")), 
          ggarrange(pu2, pu3, ncol=1, nrow=2, labels=c("c", "d"), font.label = list(size = 26, face="bold")), 
          ggarrange(pe1, pe2, nrow=1, ncol=2, labels=c("e", "f"), font.label = list(size = 26, face="bold")), nrow=1, ncol=3)
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))

