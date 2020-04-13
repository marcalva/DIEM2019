
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
# sce_ad <- readRDS("data/processed/fresh_68k/diem/fresh_68k.diem_sce.rds")
seur_diem <- readRDS("data/processed/fresh_68k/diem.nc100/fresh_68k.seur_obj.rds")
seur_quant <- readRDS("data/processed/fresh_68k/quantile/fresh_68k.seur_obj.rds")
seur_ED <- readRDS("data/processed/fresh_68k/emptydrops.nc100/fresh_68k.seur_obj.rds")

hb_genes <- c("HBA1", "HBA2", "HBB")
seur_diem <- PercentageFeatureSet(seur_diem, features=hb_genes, col.name="Hemoglobin")
seur_quant <- PercentageFeatureSet(seur_quant, features=hb_genes, col.name="Hemoglobin")
seur_ED <- PercentageFeatureSet(seur_ED, features=hb_genes, col.name="Hemoglobin")


markers_diem <- read.table("results/fresh_68k/diem.nc100/fresh_68k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_quant <- read.table("results/fresh_68k/quantile/fresh_68k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)
markers_ED <- read.table("results/fresh_68k/emptydrops.nc100/fresh_68k.seur_markers.txt", header=TRUE, stringsAsFactors=FALSE)

#=========================================
# Plot
#=========================================


dplot <- "results/plots/fresh_68k/"
dir.create(dplot, showWarnings = FALSE);

seur_diem@meta.data[,"Method"] <- "DIEM"
seur_quant@meta.data[,"Method"] <- "Quantile"
seur_ED@meta.data[,"Method"] <- "EmptyDrops"

seur_diem@meta.data <- seur_diem@meta.data[,colnames(seur_ED@meta.data)]
seur_quant@meta.data <- seur_quant@meta.data[,colnames(seur_ED@meta.data)]

dfall <- do.call(rbind, list( seur_diem@meta.data, seur_ED@meta.data, seur_quant@meta.data))
dfall[,"percent.mt"] <- dfall[,"percent.mt"] / 100
dfall[,"MALAT1"] <- dfall[,"MALAT1"] / 100
lev <- sort(as.numeric(levels(dfall[,"RNA_snn_res.0.8"])))
dfall[,"RNA_snn_res.0.8"] <- factor(dfall[,"RNA_snn_res.0.8"], 
                                    levels = lev)
dfall[,"logUMI"] <- log10(dfall[,"nCount_RNA"])


# Plots
w <- 6
h <- 1.25

methds <- c("DIEM", "EmptyDrops", "Quantile")
feats <- c("percent.mt", "logUMI")
feat_names <- c("MT%", "logUMI")

# DIEM 
for (m in methds){
    dfs <- dfall[dfall$Method == m,]
    p <- ggplot(dfs, aes(x = RNA_snn_res.0.8, y = percent.mt)) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_minimal() + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,.3)) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA), 
          plot.title = element_text(hjust = 0.5)) + 
    ylab("MT%") + 
    xlab("Fresh 68K PBMC cluster") + 
    ggtitle(m)

    ggsave(paste0(dplot, "fresh_68k.", m, ".mt_pct.pdf"), width = w, height = h)
    ggsave(paste0(dplot, "fresh_68k.", m, ".mt_pct.jpeg"), width = w, height = h)
}

for (m in methds){
    dfs <- dfall[dfall$Method == m,]
    p <- ggplot(dfs, aes(x = RNA_snn_res.0.8, y = MALAT1)) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_minimal() + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,.3)) + 
    theme(axis.title.x = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_blank())  + 
    ylab("MALAT1%") + 
    xlab("Fresh 68K PBMC cluster") + 

    ggsave(paste0(dplot, "fresh_68k.", m, ".malat1.pdf"), width = w, height = 1)
    ggsave(paste0(dplot, "fresh_68k.", m, ".malat1.jpeg"), width = w, height = 1)
}

for (m in methds){
    dfs <- dfall[dfall$Method == m,]
    p <- ggplot(dfs, aes(x = RNA_snn_res.0.8, y = logUMI)) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_minimal() + 
    theme(panel.border = element_rect(colour = "black", fill = NA)) + 
    ylim(2,4) + 
    xlab("Fresh 68K PBMC cluster") + 

    ggsave(paste0(dplot, "fresh_68k.", m, ".logUMI.pdf"), width = w, height = 1.3)
    ggsave(paste0(dplot, "fresh_68k.", m, ".logUMI.jpeg"), width = w, height = 1.3)
}

# Compare the characteristics of droplets removed by DIEM to those not
names_diff <- setdiff(rownames(seur_ED@meta.data), rownames(seur_diem@meta.data))
diem_u <- setdiff(rownames(seur_diem@meta.data), rownames(seur_ED@meta.data))
df_kept <- data.frame(seur_ED@meta.data[colnames(seur_diem),], Filter="DIEM &\nEmptyDrops")
df_rm <- data.frame(seur_ED@meta.data[names_diff,], Filter="EmptyDrops\nOnly")
df_diem <- data.frame(seur_diem@meta.data[diem_u,colnames(seur_ED@meta.data)], Filter = "DIEM\nonly")

datf_kept_rm <- do.call(rbind, list(df_kept, df_diem, df_rm))
datf_kept_rm[,"percent.mt"] <- datf_kept_rm[,"percent.mt"] / 100
datf_kept_rm[,"MALAT1"] <- datf_kept_rm[,"MALAT1"] / 100

datfm <- reshape2::melt(datf_kept_rm[,c("Filter", "percent.mt", "MALAT1")])

relabel <- c("percent.mt"="Mitochondria", "MALAT1"="MALAT1")

prmall <- ggplot(datfm, aes(x=Filter, y=value)) + 
geom_boxplot(outlier.shape = NA) + theme_bw() +
facet_wrap(~variable, labeller = labeller(variable=relabel)) + 
geom_text(data = datn, aes(x=Method, y=y, label=N, angle = 90), size=4, hjust = 1, vjust = 0.5) + 
ylab("Percent") + 
ggtitle("Fresh 68K PBMCs") + 
scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,.6)) +
theme(text = element_text(size=16),
      strip.background = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1), 
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(), 
      plot.title = element_text(hjust=0.5))

ggsave(paste0(dplot, "fresh_68k.comp.mt.malat1.pdf"), width = 5, height = 7)
ggsave(paste0(dplot, "fresh_68k.comp.mt.malat1.jpeg"), width = 5, height = 7)



tb = table(seur_ED@meta.data[shared,"RNA_snn_res.0.8"] , seur_diem@meta.data[shared,"RNA_snn_res.0.8"])

tb <- sweep(tb, 2, colSums(tb), "/")
round(tb, 2)


