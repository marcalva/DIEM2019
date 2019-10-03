
setwd("../")

library(diem)
library(reshape2)
library(ggplot2)
library(viridis)
library(gplots)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
source("scripts/common/plotting.R")

#=========================================
#=========================================
labl <- "atsn"
method <- "DE"

dp <- paste0("data/processed/", labl, "/")
dr <- paste0("results/", labl, "/", method, "/")

#=========================================
# Read in data
#=========================================

n_de <- read.table(paste0(dr, "edgeR.n_de.cell_types.txt"), 
                   row.names=1, sep="\t", stringsAsFactors=FALSE, 
                   header=TRUE, check.names=FALSE)
n_genes <- read.table(paste0(dr, "edgeR.n_genes.cell_types.txt"), 
                      row.names=1, sep="\t", stringsAsFactors=FALSE, 
                      header=TRUE, check.names=FALSE)
n_counts <- read.table(paste0(dr, "edgeR.n_counts.cell_types.txt"), 
                       row.names=1, sep="\t", stringsAsFactors=FALSE, 
                       header=TRUE, check.names=FALSE)

# Change column names
colnames(n_de) <- c("All other\ncell types", "Debris")

de <- read.table(paste0(dr, "edgeR.paired_de.diff.txt"), header=TRUE, 
                 row.names=1, sep="\t", stringsAsFactors=FALSE, 
                 check.names=FALSE)

# Volcano plot of nuclear vs debris
de$log_p<- -log10(de$p_adj)
de$DE <- "p_adj>0.05"
de[de[,"p_adj"] < 0.05, "DE"] <- "significant"

pv <- ggplot(de, aes(x=logFC, y=log_p, color=DE)) +
geom_point(shape=16) +
theme_minimal() +
xlim(c(-max(de$logFC), max(de$logFC))) +
scale_color_brewer(palette="Set1", breaks="significant", 
                   labels="significant", name="") +
ylab(expression(-log[10](P))) +
ggtitle("DE between debris\nand nuclei") + 
theme(text=element_text(size=18), 
      plot.title = element_text(hjust=0.5))

# Heatmap of number of DE genes between cell types and debris
n_de_m <- na.omit(melt(as.matrix(n_de*100)))
n_de_m$Var1 <- factor(n_de_m$Var1)
ph <- ggplot(n_de_m, aes(x=Var1, y=Var2)) + 
geom_tile(aes(fill=value), colour="white") + 
scale_fill_viridis(name="Percent of\ngenes DE") + 
theme_minimal() + 
ylab("Cell type") + 
theme(panel.grid=element_blank(), 
      panel.border = element_blank(), 
      text=element_text(size=18), 
      axis.title.x = element_blank())

# Boxplot number of DE genes between pairs
n_de_ct <- n_de[,1]
#n_de_ct <- n_de[1:((ncol(n_de)-1)), 1:((ncol(n_de)-1))]
#n_de_ct <- n_de_ct[upper.tri(n_de_ct)]
n_de_ct <- data.frame("DE" = 100*n_de_ct, "Pair" = "Cell type-cell type")
#n_de_dbr <- n_de[1:((ncol(n_de)-1)), ncol(n_de)]
n_de_dbr <- n_de[,2]
n_de_dbr <- data.frame("DE" = 100*n_de_dbr, "Pair" = "Cell type-debris")

pval <- t.test(n_de[,1], n_de[,2])$p.value

datf_de <- rbind(n_de_ct, n_de_dbr)
datf_de$Pair <- factor(datf_de$Pair, 
                       levels=c("Cell type-debris", "Cell type-cell type"))
pbde <- ggplot(datf_de, aes(x = Pair, y = DE)) + 
geom_boxplot() + geom_jitter(width=0.1) + 
annotate(geom="text", x = 1.5, y = 13, label = paste0("p=", as.character(round(pval, 2))), 
         size = 8) + 
theme_minimal() + 
ylab("Percent of genes DE") + 
ggtitle("DE between debris\n and cell types") + 
theme(text=element_text(size=18), 
      axis.title.x = element_blank(), 
      panel.grid.major.x = element_blank(), 
      plot.title = element_text(hjust=0.5))

# Boxplot number of number of counts
#n_c_ct <- n_counts[1:((ncol(n_counts)-1)), 1:((ncol(n_counts)-1))]
#n_c_ct <- n_c_ct[upper.tri(n_c_ct)]
#n_c_dbr <- n_counts[1:((ncol(n_counts)-1)), ncol(n_counts)]

dir_plot <- "results/plots/"; dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)

pdfname <- paste0(dir_plot, "de_genes.ct.pdf")
jpgname <- paste0(dir_plot, "de_genes.ct.jpeg")
pdf(pdfname, width=12,height=10)
ggarrange(
	ggarrange(pv, pbde, ncol=2, widths=c(0.6, 0.4), labels=c("a", "b"), font.label = list(size = 24, face="bold")), 
    ph, 
	nrow=2, ncol=1, heights=c(.7, .3), 
    labels=c("", "c"), font.label = list(size = 24, face="bold"))
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))

