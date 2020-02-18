
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

get_mt_reduction <- function(x, seur1, seur2, thresh=0.25, ct_col1="seurat_clusters", ct_col2="seurat_clusters"){
	clusters1 <- as.character(seur1@meta.data[,ct_col1])
	clusters2 <- as.character(seur2@meta.data[,ct_col2])

	ct1 <- as.character(sort(as.numeric(unique(clusters1))))
	ct2 <- as.character(sort(as.numeric(unique(clusters2))))

	mt1 <- sapply(ct1, function(x){mean(seur1@meta.data[seur1$seurat_clusters == x,"SpliceFrctn"], na.rm=T)})
	mt2 <- sapply(ct2, function(x){mean(seur2@meta.data[seur2$seurat_clusters == x,"SpliceFrctn"], na.rm=T)})

	mt_diff <- c()
	for (c2 in ct2){
		c1 <- as.character(rownames(x)[which.max(x[,c2])])
		if (x[c1,c2] < thresh) mt_diff[c2] <- NA
		else mt_diff[c2] <- (mt2[c2]-mt1[c1])
	}
	return(mt_diff)
}

bp <- function(x, y="Change in background RNA"){
    ret <- ggplot(x, aes(x=0, y=x)) + 
    geom_boxplot() + 
    facet_wrap(~Experiment) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
    theme_classic() + 
    ylab(y) + 
    theme(axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          text=element_text(size=26))
    return(ret)
}

p_blank <- ggplot() + theme_void()

#=========================================
#=========================================

#=========================================
# Read in data
#=========================================

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

# Adipocyte
seur_diem_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.seur_obj.rds")
seur_diem_ad$SpliceFrctn <- 100*sf[rownames(seur_diem_ad@meta.data),"SpliceFrctn"]
seur_quant_ad <- readRDS("data/processed/adpcyte/quantile/adpcyte.seur_obj.rds")
seur_quant_ad$SpliceFrctn <- 100*sf[rownames(seur_quant_ad@meta.data),"SpliceFrctn"]
seur_ED_ad <- readRDS("data/processed/adpcyte/emptydrops/adpcyte.seur_obj.rds")
seur_ED_ad$SpliceFrctn <- 100*sf[rownames(seur_ED_ad@meta.data),"SpliceFrctn"]

# Mouse brain
seur_diem_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.seur_obj.rds")
seur_diem_mb$SpliceFrctn <- 100*sf[rownames(seur_diem_mb@meta.data),"SpliceFrctn"]
seur_quant_mb <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.seur_obj.rds")
seur_quant_mb$SpliceFrctn <- 100*sf[rownames(seur_quant_mb@meta.data),"SpliceFrctn"]
seur_ED_mb <- readRDS("data/processed/mouse_nuclei_2k/emptydrops/mouse_nuclei_2k.seur_obj.rds")
seur_ED_mb$SpliceFrctn <- 100*sf[rownames(seur_ED_mb@meta.data),"SpliceFrctn"]

# Adipose tissue
seur_diem_at <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_diem_at$SpliceFrctn <- 100*sf[rownames(seur_diem_at@meta.data),"SpliceFrctn"]
seur_quant_at <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
seur_quant_at$SpliceFrctn <- 100*sf[rownames(seur_quant_at@meta.data),"SpliceFrctn"]
seur_ED_at <- readRDS("data/processed/atsn/emptydrops/atsn.seur_obj.rds")
seur_ED_at$SpliceFrctn <- 100*sf[rownames(seur_ED_at@meta.data),"SpliceFrctn"]

seur_ad = list("DIEM" = seur_diem_ad, "EmptyDrops" = seur_ED_ad, "Quantile" = seur_quant_ad)
seur_mb = list("DIEM" = seur_diem_mb, "EmptyDrops" = seur_ED_mb, "Quantile" = seur_quant_mb)
seur_at = list("DIEM" = seur_diem_at, "EmptyDrops" = seur_ED_at, "Quantile" = seur_quant_at)

#=========================================
# Plot
#=========================================

coltitle <- "% reads\nspliced"

# Get number of nuclei
methd_names <- c("Quantile", "EmptyDrops", "DIEM")
methd <- factor(c(methd_names, methd_names, methd_names), levels=c("Quantile", "EmptyDrops", "DIEM"))
exper <- factor(c("DiffPA", "DiffPA", "DiffPA", 
		   "Mouse\nBrain", "Mouse\nBrain", "Mouse\nBrain", 
		   "Adipose\nTissue", "Adipose\nTissue", "Adipose\nTissue"),
                levels=c("DiffPA", "Mouse\nBrain", "Adipose\nTissue"))
ncel <- c(ncol(seur_quant_ad), ncol(seur_ED_ad), ncol(seur_diem_ad), 
		  ncol(seur_quant_mb), ncol(seur_ED_mb), ncol(seur_diem_mb), 
		  ncol(seur_quant_at), ncol(seur_ED_at), ncol(seur_diem_at))
datf <- data.frame(Method=methd, Experiment=exper, Nuclei=ncel)

# Cluster overlap between EmptyDrops and DIEM
ctype="RNA_snn_res.0.8"
ad_overlap <- comp_clust(seur_ED_ad, seur_diem_ad, ct_col1=ctype, ct_col2=ctype)
mt_ED_ad <- sapply(rownames(ad_overlap), function(x){
                   mean(seur_ED_ad@meta.data[seur_ED_ad@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})
mt_diem_ad <- sapply(colnames(ad_overlap), function(x){
                     mean(seur_diem_ad@meta.data[seur_diem_ad@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})

mb_overlap <- comp_clust(seur_ED_mb, seur_diem_mb, ct_col1=ctype, ct_col2=ctype)
mt_ED_mb <- sapply(rownames(mb_overlap), function(x){
                   mean(seur_ED_mb@meta.data[seur_ED_mb@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})
mt_diem_mb <- sapply(colnames(mb_overlap), function(x){
                     mean(seur_diem_mb@meta.data[seur_diem_mb@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})

at_overlap <- comp_clust(seur_ED_at, seur_diem_at, ct_col1=ctype, ct_col2=ctype)
mt_ED_at <- sapply(rownames(at_overlap), function(x){
                   mean(seur_ED_at@meta.data[seur_ED_at@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})
mt_diem_at <- sapply(colnames(at_overlap), function(x){
                     mean(seur_diem_at@meta.data[seur_diem_at@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})

pe1 <- overlap_graph(ad_overlap, mt_ED_ad, mt_diem_ad, labels1="EmptyDrops", labels2="DIEM", main="DiffPA", col_title=coltitle)
pe2 <- overlap_graph(mb_overlap, mt_ED_mb, mt_diem_mb, labels1="EmptyDrops", labels2="DIEM", main="Mouse Brain", col_title=coltitle)
pe3 <- overlap_graph(at_overlap, mt_ED_at, mt_diem_at, labels1="EmptyDrops", labels2="DIEM", main="Adipose Tissue", col_title=coltitle)

# Average MT% reduction
r_ed1 <- data.frame(x=get_mt_reduction(ad_overlap, seur_ED_ad, seur_diem_ad), Experiment="DiffPA")
r_ed2 <- data.frame(x=get_mt_reduction(mb_overlap, seur_ED_mb, seur_diem_mb), Experiment="Mouse Brain")
r_ed3 <- data.frame(x=get_mt_reduction(at_overlap, seur_ED_at, seur_diem_at), Experiment="Adipose Tissue")
r_ed_all <- do.call(rbind, list(r_ed1, r_ed2, r_ed3))
r_ed_allp <- bp(r_ed_all, y=expression("SF"[DIEM] - "SF"[EmptyDrops]))

# Cluster overlap between quantile and DIEM
ctype="RNA_snn_res.0.8"
ad_overlap <- comp_clust(seur_quant_ad, seur_diem_ad, ct_col1=ctype, ct_col2=ctype)
mt_quant_ad <- sapply(rownames(ad_overlap), function(x){mean(seur_quant_ad@meta.data[seur_quant_ad@meta.data[,ctype] == x,"SpliceFrctn"])})
mt_diem_ad <- sapply(colnames(ad_overlap), function(x){mean(seur_diem_ad@meta.data[seur_diem_ad@meta.data[,ctype] == x,"SpliceFrctn"])})

mb_overlap <- comp_clust(seur_quant_mb, seur_diem_mb, ct_col1=ctype, ct_col2=ctype)
mt_quant_mb <- sapply(rownames(mb_overlap), function(x){mean(seur_quant_mb@meta.data[seur_quant_mb@meta.data[,ctype] == x,"SpliceFrctn"])})
mt_diem_mb <- sapply(colnames(mb_overlap), function(x){mean(seur_diem_mb@meta.data[seur_diem_mb@meta.data[,ctype] == x,"SpliceFrctn"])})

at_overlap <- comp_clust(seur_quant_at, seur_diem_at, ct_col1=ctype, ct_col2=ctype)
mt_quant_at <- sapply(rownames(at_overlap), function(x){mean(seur_quant_at@meta.data[seur_quant_at@meta.data[,ctype] == x,"SpliceFrctn"])})
mt_diem_at <- sapply(colnames(at_overlap), function(x){mean(seur_diem_at@meta.data[seur_diem_at@meta.data[,ctype] == x,"SpliceFrctn"])})

pq1 <- overlap_graph(ad_overlap, mt_quant_ad, mt_diem_ad, labels1="Quantile", labels2="DIEM", main="DiffPA", col_title=coltitle)
pq2 <- overlap_graph(mb_overlap, mt_quant_mb, mt_diem_mb, labels1="Quantile", labels2="DIEM", main="Mouse Brain", col_title=coltitle)
pq3 <- overlap_graph(at_overlap, mt_quant_at, mt_diem_at, labels1="Quantile", labels2="DIEM", main="Adipose Tissue", col_title=coltitle)


# Average MT% reduction
r_q1 <- data.frame(x=get_mt_reduction(ad_overlap, seur_quant_ad, seur_diem_ad), Experiment="DiffPA")
r_q2 <- data.frame(x=get_mt_reduction(mb_overlap, seur_quant_mb, seur_diem_mb), Experiment="Mouse Brain")
r_q3 <- data.frame(x=get_mt_reduction(at_overlap, seur_quant_at, seur_diem_at), Experiment="Adipose Tissue")
r_q_all <- do.call(rbind, list(r_q1, r_q2, r_q3))
r_q_allp <- bp(r_q_all, y=expression("%Spliced"[DIEM] - "%Spliced"[Quantile]))

# Arrange into plot
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "clust_overlaps.quant.ED.pdf")
jpgname <- paste0(dir_plot, "clust_overlaps.quant.ED.jpeg")
pdf(pdfname, width=16, height=9)
ggarrange(p_blank, pq1, pq2, pq3, p_blank, pe1, pe2, pe3, 
          ncol=8, labels=c("a", "", "", "", "b", "", "", ""), 
          font.label = list(size = 28, face="bold"),
          widths=c(.25,1,1,1,.25,1,1,1))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


# Average MT% reduction
mean(get_mt_reduction(ad_overlap, seur_quant_ad, seur_diem_ad), na.rm=T)
mean(get_mt_reduction(mb_overlap, seur_quant_mb, seur_diem_mb), na.rm=T)
mean(get_mt_reduction(at_overlap, seur_quant_at, seur_diem_at), na.rm=T)


# Higher resolution in DIEM
seur_diem_at <- FindClusters(seur_diem_at, resolution=2)
seur_diem_mb <- FindClusters(seur_diem_mb, resolution=2)
#seur_ED_ad <- FindClusters(seur_ED_ad, resolution=2)
#seur_ED_mb <- FindClusters(seur_ED_mb, resolution=2)
#seur_ED_at <- FindClusters(seur_ED_at, resolution=2)

# Cluster overlap
ctype="RNA_snn_res.0.8"
ctype2="RNA_snn_res.2"

mbq_overlap <- comp_clust(seur_quant_mb, seur_diem_mb, ct_col1=ctype, ct_col2=ctype2)
mtq_quant_mb <- sapply(rownames(mbq_overlap), function(x){mean(seur_quant_mb@meta.data[seur_quant_mb@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})
mtq_diem_mb <- sapply(colnames(mbq_overlap), function(x){mean(seur_diem_mb@meta.data[seur_diem_mb@meta.data[,ctype2] == x,"SpliceFrctn"], na.rm=T)})

mb_overlap <- comp_clust(seur_ED_mb, seur_diem_mb, ct_col1=ctype, ct_col2=ctype2)
mt_quant_mb <- sapply(rownames(mb_overlap), function(x){mean(seur_ED_mb@meta.data[seur_ED_mb@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})
mt_diem_mb <- sapply(colnames(mb_overlap), function(x){mean(seur_diem_mb@meta.data[seur_diem_mb@meta.data[,ctype2] == x,"SpliceFrctn"], na.rm=T)})

atq_overlap <- comp_clust(seur_quant_at, seur_diem_at, ct_col1=ctype, ct_col2=ctype2)
atq_quant_at <- sapply(rownames(atq_overlap), function(x){mean(seur_quant_at@meta.data[seur_quant_at@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})
atq_diem_at <- sapply(colnames(atq_overlap), function(x){mean(seur_diem_at@meta.data[seur_diem_at@meta.data[,ctype2] == x,"SpliceFrctn"], na.rm=T)})

at_overlap <- comp_clust(seur_ED_at, seur_diem_at, ct_col1=ctype, ct_col2=ctype2)
mt_ED_at <- sapply(rownames(at_overlap), function(x){mean(seur_ED_at@meta.data[seur_ED_at@meta.data[,ctype] == x,"SpliceFrctn"], na.rm=T)})
mt_diem_at <- sapply(colnames(at_overlap), function(x){mean(seur_diem_at@meta.data[seur_diem_at@meta.data[,ctype2] == x,"SpliceFrctn"], na.rm=T)})


ps1 <- overlap_graph(mbq_overlap, mtq_quant_mb, mtq_diem_mb, labels1="Quantile", labels2="DIEM (Res. 2)", main="Mouse Brain", col_title=coltitle)
ps2 <- overlap_graph(mb_overlap, mt_ED_mb, mt_diem_mb, labels1="EmptyDrops", labels2="DIEM (Res. 2)", main="Mouse Brain", col_title=coltitle)
ps3 <- overlap_graph(atq_overlap, atq_quant_at, atq_diem_at, labels1="Quantile", labels2="DIEM (Res. 2)", main="Adipose Tissue", col_title=coltitle)
ps4 <- overlap_graph(at_overlap, mt_ED_at, mt_diem_at, labels1="EmptyDrops", labels2="DIEM (Res. 2)", main="Adipose Tissue", col_title=coltitle)

# Arrange into plot
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "clust_overlaps.res2.pdf")
jpgname <- paste0(dir_plot, "clust_overlaps.res2.jpeg")
pdf(pdfname, width=20, height=14)
ggarrange(ps1, ps2, ps3, ps4, nrow=1, 
          labels=c("a", "", "b", ""),
          font.label = list(size = 30, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))



