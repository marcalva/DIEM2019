
setwd("../")

library(diem)
library(Seurat)
library(ggplot2)
library(gplots)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(GGally)
library(RColorBrewer)
library(biomaRt)
source("scripts/common/plotting.R")


#=========================================
# Functions
#=========================================



norm_tmm <- function(counts, thresh=1){
    require(edgeR)
    y = edgeR::DGEList(counts=counts)
    keep <- rowSums(cpm(y) >= thresh) == 2
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)

    counts_norm <- sapply(1:2, function(i){ y$samples[i,"norm.factors"] * 1e6 * y$counts[,i] / sum(y$counts[,i]) })
    rownames(counts_norm) <- rownames(y$counts)
    colnames(counts_norm) <- colnames(y$counts)

    return(counts_norm)
}

get_lh <- function(sce, cutpoint=100, cpm_thresh=3, prefix=NULL){
    tc <- Matrix::colSums(sce@counts)

    low_counts <- Matrix::rowSums(sce@counts[,tc < cutpoint])
    high_counts <- Matrix::rowSums(sce@counts[,tc >= cutpoint])
    lh <- data.frame(low=low_counts, high=high_counts)
    lh <- norm_tmm(lh)

    if (!is.null(prefix)){
        colnames(lh) <- paste0(prefix, " ", colnames(lh))
    }
    return(lh)
}

get_log_fc <- function(sce, cutpoint=100, cpm_thresh=3){
    lh <- get_lh(sce)
    logfc <- log2(lh[,1]) - log2(lh[,2])
    return(logfc)
}

panel.cor <- function(x, y, digits = 2, cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  

  par(usr = c(0, 1, 0, 1)) 
  # correlation coefficient
  r <- cor(x, y, use="pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "") 
  text(0.5, 0.6, txt)

  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "") 
  if(is.na(p)){txt2 <- NA}
  else if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "") 
  text(0.5, 0.4, txt2)
}

convertMouseGeneList <- function(vec){
  x <- names(vec)
  require("biomaRt")
  human_ds = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse_ds = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  converted = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
    values = x , mart = mouse_ds, attributesL = c("hgnc_symbol"), martL = human_ds, uniqueRows=T)
    
  # Remove duplicated
  keep1 <- ! ( duplicated(converted[,1]) | duplicated(converted[,1], fromLast=TRUE) )
  keep2 <- ! ( duplicated(converted[,2]) | duplicated(converted[,2], fromLast=TRUE) )
  keep <- keep1 & keep2
  converted <- converted[keep,]
  rownames(converted) <- converted[,"MGI.symbol"]

  genes <- intersect(converted$MGI.symbol, x)
  vec <- vec[genes]
  names(vec) <- converted[genes, "HGNC.symbol"]
  return(vec)
}

#=========================================
# Read in data
#=========================================

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

# Adipocyte
labl <- "adpcyte"
sce_ad <- readRDS("data/processed/adpcyte/diem/adpcyte.diem_sce.rds")
quant <- readRDS("data/processed/adpcyte/quantile/adpcyte.quantile_counts.rds")
thresh_ad <- quant$thresh

log2fc_ad <- get_log_fc(sce_ad, 100)

# Mouse brain
labl <- "mouse_nuclei_2k"
sce_mb <- readRDS("data/processed/mouse_nuclei_2k/diem/mouse_nuclei_2k.diem_sce.rds")
quant <- readRDS("data/processed/mouse_nuclei_2k/quantile/mouse_nuclei_2k.quantile_counts.rds")
thresh_mb <- quant$thresh

log2fc_mb <- get_log_fc(sce_mb, 100)

# Adipose tissue
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

labl <- "atsn"
scel_at <- lapply(lab_ids, function(id) readRDS(paste0("data/processed/atsn/diem/", id, ".diem_sce.rds")) )
names(scel_at) <- lab_ids
thresh_at <- lapply(lab_ids, function(id) {
  quant <- readRDS(paste0("data/processed/atsn/quantile/", id, ".quantile_counts.rds"))
  return(quant$thresh)
})
names(thresh_at) <- lab_ids

log2fc_at <- lapply(lab_ids, function(id) {
  log2fc <- get_log_fc(scel_at[[id]], 100)
  return(log2fc)
})
names(log2fc_at) <- lab_ids

#=========================================
#=========================================

#=========================================
# Cor log2fc
#=========================================

log2fc_mb <- convertMouseGeneList(log2fc_mb) # Convert to human

log2fc <- c(DiffPA=list(log2fc_ad), MouseBrain=list(log2fc_mb), log2fc_at)

genes <- lapply(log2fc, function(x) names(x))
cg <- Reduce(intersect, genes)

log2fc <- lapply(log2fc, function(x) x[cg])

datf <- as.data.frame(do.call(cbind, log2fc))

colnames(datf)[1:2] <- c("DiffPA", "Mouse\nBrain")


p <- GGally::ggpairs(datf, axisLabels="internal") + 
theme_minimal()

pdfname <- paste0("results/plots/pairs.log2fc_cor.pdf")
jpgname <- paste0("results/plots/pairs.log2fc_cor.jpeg")
pdf("results/plots/pairs.log2fc_cor.pdf", width=10, height=10)
p
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


# What is the correlation between adipocyte and AT
rs <- cor(datf[,3:8], datf[,1])
summary(rs)

print(mean(cor(datf)))
