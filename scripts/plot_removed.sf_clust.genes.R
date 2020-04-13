
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
#=========================================

boxplot_ft_fct <- function(x, names, ct_col="RNA_snn_res.0.8", 
						   colname="percent.mt", ylab="MT%", 
						   color_breaks=waiver(), scales = "free_x", ncol = 1,
						   size=1){
	dfl <- lapply(1:length(x), function(i) {
                  data.frame(Mito=x[[i]]@meta.data[,colname], 
                             Cluster=x[[i]]@meta.data[,ct_col], 
                             Method=names[[i]])
                           })
	df <- do.call(rbind, dfl)
	df <- reshape2::melt(df)
    names <- factor(names)
    df[,"Method"] <- factor(df[,"Method"], levels = levels(names))

	p <- ggplot(df, aes(x=Cluster, y=value)) +
	geom_boxplot(outlier.shape = NA) + theme_bw() + 
	facet_wrap(~Method, scales=scales, ncol = ncol) + 
	ylab(ylab) + 
	theme(legend.position="none", 
		  text=element_text(size=10),
		  axis.text.x=element_text(size=8, angle = 30, hjust = 0.9), 
		  plot.title=element_text(hjust=0.5)) 
	return(p)
}

boxplot_gene_l <- function(x, names, gene, ct_col="CellType", 
                           scales = "free", 
                           ncol = 1,
                           size=1){
	dfl <- lapply(1:length(x), function(i) {
                  data.frame(Gene = x[[i]]@assays$RNA@data[gene,], 
                             Cluster = x[[i]]@meta.data[,ct_col], 
                             Method=names[[i]])
                           })
	df <- do.call(rbind, dfl)
	df <- reshape2::melt(df)
    names <- factor(names)
    df[,"Method"] <- factor(df[,"Method"], levels = levels(names))
    cts <- sort(unique(as.character(df[,"Cluster"])))
    df[,"Cluster"] <- factor(df[,"Cluster"], levels = cts)

	p <- ggplot(df, aes(x=Cluster, y=value)) +
	geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.2) + 
    theme_bw() + 
	facet_wrap(~Method, scales=scales, ncol = ncol) + 
	ylab("logUMI") + 
    ggtitle(gene) + 
	theme(legend.position="none", 
		  text=element_text(size=10),
          strip.background = element_blank(),
		  axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size=8), 
		  plot.title=element_text(hjust=0)) 
	return(p)
}

plot_umap_gene_l <- function(x, 
                             names, 
                           gene, 
                           legend_title="UMI", 
                           low = "lightgrey", high = "firebrick3", 
                           logt = FALSE,
                           assay=NULL, 
                           alpha=1, 
                           size=.5,
                           ncol = 1, 
                           shape = 16, 
                           rand=TRUE){
    require(ggplot2)
    require(RColorBrewer)

    datf <- lapply(1:length(x), function(i){
                   s <- x[[i]]
                   n <- names[i]
                   expr <- GetAssayData(s, slot="data", assay=assay)
                   datf <- data.frame(feat=expr[gene,], 
                                      UMAP1=s@reductions$umap@cell.embeddings[,"UMAP_1"],
                                      UMAP2=s@reductions$umap@cell.embeddings[,"UMAP_2"], 
                                      Method=n)
                           })

    datf <- do.call(rbind, datf)
    if (logt){
        datf[,"feat"] <- log1p(datf[,"feat"])
    }
    if (rand) datf <- datf[sample(1:nrow(datf)),]
    p <- ggplot(datf, aes(x=UMAP1, y=UMAP2, color=feat)) + 
    geom_point(alpha=alpha, shape = shape, size = size) + 
    theme_classic() + 
    ggtitle(gene) + 
    facet_wrap(~Method, scales = "free", ncol = ncol) + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(hjust = 0), 
          legend.text = element_text(size = 6), 
          legend.key.size = unit(.08, "in")) + 
    scale_color_gradient(low=low, high=high, name=legend_title)
    return(p)
}



# Variables
methd_names <- c("Quantile", "EmptyDrops", "DIEM")
methd <- factor(c(methd_names, methd_names, methd_names), levels=c("Quantile", "EmptyDrops", "DIEM"))


#=========================================
# Read in data
#=========================================

# Adipose tissue
labl <- "atsn"
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
dp <- paste0("data/processed/", labl, "/")



seur_diem_at.clean <- readRDS("data/processed/atsn/diem/atsn.seur_obj.rds")
seur_diem_at.debris <- readRDS("data/processed/atsn/diem_debris/atsn.seur_obj.rds")

seur_quant_at.clean <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")
seur_quant_at.debris <- readRDS("data/processed/atsn/quantile_debris/atsn.seur_obj.rds")

seur_ED_at.clean <- readRDS("data/processed/atsn/emptydrops/atsn.seur_obj.rds")
seur_ED_at.debris <- readRDS("data/processed/atsn/emptydrops_debris/atsn.seur_obj.rds")


#=========================================
#=========================================

# Box plot of % spliced reads
methd_names <- c("DIEM\nremoved", "EmptyDrops\nremoved", "Quantile\nremoved")
methd_names <- factor(methd_names, levels = methd_names)
sl <- list("DIEM" = seur_diem_at.debris, "EmptyDrops" = seur_ED_at.debris, "Quantile" = seur_quant_at.debris)
pbox <- boxplot_ft_fct(sl, methd_names, scales = "free_x", ct_col = "CellType", 
                       colname="SpliceFrctn", ylab = "Percent reads spliced", 
                       ncol = 3)

fn <- "results/plots/sf.celltype.method.box.removed.pdf"
ggsave(fn, width = 8, height = 3)
fn <- "results/plots/sf.celltype.method.box.removed.jpeg"
ggsave(fn, width = 8, height = 3)



seur_l <- list(seur_diem_at.debris, seur_ED_at.debris, seur_quant_at.debris)
m <- c("DIEM\nremoved", "EmptyDrops\nremoved", "Quantile\nremoved")

w <- 2
h <- 5

genes <- c("PLIN1", 
           "PDE3B", 
           "VWF", 
           "PTPRB", 
           "SKAP1", 
           "CD14", 
           "ALCAM", 
           "FBN1", 
           "COL1A2")

dir_plot <- "results/plots/genes_debris/"
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

for (g in genes){
    pg1 <- plot_umap_gene_l(seur_l, m, g)

    fn <- paste0(dir_plot, "umap.removed.", g, ".pdf")
    ggsave(fn, width = w, height = h)
    fn <- paste0(dir_plot, "umap.removed.", g, ".jpeg")
    ggsave(fn, width = w, height = h)
}


w <- 4
h <- 5

for (g in genes){
    pg1 <- boxplot_gene_l(seur_l, m, g)

    fn <- paste0(dir_plot, "box.removed.", g, ".pdf")
    ggsave(fn, width = w, height = h)
    fn <- paste0(dir_plot, "box.removed.", g, ".jpeg")
    ggsave(fn, width = w, height = h)

}

