
# Plot results from varying filter in adipose tissue

setwd("../")

library(diem)
library(ggplot2)
library(ggpubr)

#====================================================
# Adipose
#====================================================

dp <- "data/processed/atsn/diem/"
dr <- "results/atsn/diem/"
dir.create(dr, showWarnings=FALSE, recursive=TRUE)

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
k_vals <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
sf[,1] <- 100 * sf[,1]

scel <- lapply(lab_ids, function(i){
    scefn <- paste0(dp, i, ".diem_sce.rds")
    sce <- readRDS(scefn)
    mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
    genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")
    md <- sce@droplet_data
    rownames(md) <- paste0(i, "_", rownames(md))
    return(md)
})
names(scel) <- lab_ids

kt <- lapply(k_vals, function(k){
    kc <- as.character(k)
    mdkl <- lapply(lab_ids, function(i){
        ifn <- paste0(dp, "filt/meta.data.", i, ".", kc, ".txt")
        mdk <- read.table(ifn, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
        rownames(mdk) <- paste0(i, "_", rownames(mdk))
        mdk[,"pct.mt"] <- scel[[i]][rownames(mdk),"pct.mt"]
        mdk[,"MALAT1"] <- scel[[i]][rownames(mdk),"MALAT1"]
        keep <- intersect(rownames(sf), rownames(mdk))
        mdk[keep,"SpliceFrctn"] <- sf[keep,1]
        mdk <- mdk[mdk$n_genes >= 200, ]
        mdk[,"ID"] <- i
        return(mdk)
    })
    mdkdf <- do.call(rbind, mdkl)
    return(mdkdf)
})
names(kt) <- k_vals

summ <- sapply(kt, function(datk){
    cl <- datk$Call == "Clean"
    mt <- mean(datk[cl,"pct.mt"])
    malat <- mean(datk[cl,"MALAT1"])
    ngene <- mean(datk[cl,"n_genes"])
    sfk <- mean(datk[cl,"SpliceFrctn"])
    nclust <- mean(tapply(datk[cl,"Cluster"], datk[cl, "ID"], function(i){length(unique(i))}))
    npass <- sum(cl)
    ntotal <- nrow(datk)
    nfilter <- 100 * npass / ntotal
    return(c(npass, nclust, nfilter, mt, malat, ngene, sfk))
})

summ <- t(summ)
summ <- as.data.frame(summ)
colnames(summ) <- c("NumPass", "Final K", "PercentPass", "AvgMT%", "AvgMALAT1%", "AvgNumGenes", "SpliceFraction")
rownames(summ) <- k_vals
summ[,"K"] <- k_vals

summ_at <- summ[,c("NumPass", "Final K", "PercentPass", "SpliceFraction", "K")]

summl_at <- reshape2::melt(summ_at, id.vars = "K")

p1 <- ggplot(summl_at, aes(x = K, y = value)) +
geom_line() +
facet_wrap(~variable, scales="free", nrow=1) +
ylab("") +
xlab("Filter") +
theme_minimal() +
ggtitle(expression(paste("Adipose Tissue ( ", K[init], " = 30)"))) +
theme(plot.title = element_text(hjust = 0.5))

# Keep Filt = 0.25

k_vals <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
kt <- lapply(k_vals, function(k){
    kc <- as.character(k)
    mdkl <- lapply(lab_ids, function(i){
        fn <- paste0("meta.data.", i, k, ".txt")
        ifn <- file.path("data/processed/atsn/diem/k_init/", fn)
        mdk <- read.table(ifn, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
        rownames(mdk) <- paste0(i, "_", rownames(mdk))
        mdk[,"pct.mt"] <- scel[[i]][rownames(mdk),"pct.mt"]
        mdk[,"MALAT1"] <- scel[[i]][rownames(mdk),"MALAT1"]
        keep <- intersect(rownames(sf), rownames(mdk))
        mdk[keep,"SpliceFrctn"] <- sf[keep,1]
        mdk <- mdk[mdk$n_genes >= 200, ]
        mdk[,"ID"] <- i
        return(mdk)
    })
    mdkdf <- do.call(rbind, mdkl)
    return(mdkdf)
})
names(kt) <- k_vals


summ <- sapply(kt, function(datk){
    cl <- datk$Call == "Clean"
    mt <- mean(datk[cl,"pct.mt"])
    malat <- mean(datk[cl,"MALAT1"])
    ngene <- mean(datk[cl,"n_genes"])
    sfk <- mean(datk[cl,"SpliceFrctn"])
    nclust <- mean(tapply(datk[cl,"Cluster"], datk[cl, "ID"], function(i){length(unique(i))}))
    npass <- sum(cl)
    ntotal <- nrow(datk)
    nfilter <- 100 * npass / ntotal
    return(c(npass, nclust, nfilter, mt, malat, ngene, sfk))
})

summ <- t(summ)
summ <- as.data.frame(summ)
colnames(summ) <- c("NumPass", "Final K", "PercentPass", "AvgMT%", "AvgMALAT1%", "AvgNumGenes", "SpliceFraction")
rownames(summ) <- k_vals
summ[,"K"] <- k_vals


summ_at <- summ[,c("NumPass", "Final K", "PercentPass", "SpliceFraction", "K")]

summl_at <- reshape2::melt(summ_at, id.vars = "K")

p4 <- ggplot(summl_at, aes(x = K, y = value)) +
geom_line() +
facet_wrap(~variable, scales="free", nrow=1) +
ylab("") +
xlab(expression(K[init])) +
theme_minimal() +
ggtitle("Adipose Tissue (Filter = 0.2)") +
theme(plot.title = element_text(hjust = 0.5))

#====================================================
# DiffPA
#====================================================


dp <- "data/processed/adpcyte/diem/"
dr <- "results/adpcyte/diem/"
dir.create(dr, showWarnings=FALSE, recursive=TRUE)

k_vals <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
sf[,1] <- 100 * sf[,1]


scefn <- paste0(dp, "adpcyte.diem_sce.rds")
sce <- readRDS(scefn)
mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")
sce <- sce@droplet_data
rownames(sce) <- paste0("adpcyte_", rownames(sce))

kt <- lapply(k_vals, function(k){
    kc <- as.character(k)
    ifn <- paste0(dp, "filt/meta.data.", kc, ".txt")
    mdk <- read.table(ifn, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    rownames(mdk) <- paste0("adpcyte_", rownames(mdk))
    mdk[,"pct.mt"] <- sce[rownames(mdk),"pct.mt"]
    mdk[,"MALAT1"] <- sce[rownames(mdk),"MALAT1"]
    keep <- intersect(rownames(sf), rownames(mdk))
    mdk[keep,"SpliceFrctn"] <- sf[keep,"SpliceFrctn"]
    mdk <- mdk[mdk$n_genes >= 200, ]
    return(mdk)
})
names(kt) <- k_vals

summ <- sapply(kt, function(datk){
    cl <- datk$Call == "Clean"
    mt <- mean(datk[cl,"pct.mt"])
    malat <- mean(datk[cl,"MALAT1"])
    ngene <- mean(datk[cl,"n_genes"])
    sfk <- mean(datk[cl,"SpliceFrctn"], na.rm = T)
    nclust <- length(unique(datk[cl,"Cluster"]))
    npass <- sum(cl)
    ntotal <- nrow(datk)
    nfilter <- 100 * npass / ntotal
    return(c(npass, nclust, nfilter, mt, malat, ngene, sfk))
})

summ <- t(summ)
summ <- as.data.frame(summ)
colnames(summ) <- c("NumPass", "Final K", "PercentPass", "AvgMT%", "AvgMALAT1%", "AvgNumGenes", "SpliceFraction")
rownames(summ) <- k_vals
summ[,"K"] <- k_vals

summ_ad <- summ[,c("NumPass", "Final K", "PercentPass", "SpliceFraction", "K")]

summl_ad <- reshape2::melt(summ_ad, id.vars = "K")

p2 <- ggplot(summl_ad, aes(x = K, y = value)) +
geom_line() +
facet_wrap(~variable, scales="free", nrow=1) +
ylab("") +
xlab("Filter") +
theme_minimal() +
ggtitle(expression(paste("DiffPA ( ", K[init], " = 30)"))) +
theme(plot.title = element_text(hjust = 0.5))

#====================================================
# Mouse brain
#====================================================

dp <- "data/processed/mouse_nuclei_2k/diem/"
dr <- "results/mouse_nuclei_2k/diem/"
dir.create(dr, showWarnings=FALSE, recursive=TRUE)

k_vals <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
sf[,1] <- 100 * sf[,1]

scefn <- paste0(dp, "mouse_nuclei_2k.diem_sce.rds")
sce <- readRDS(scefn)
mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")
sce <- sce@droplet_data
rownames(sce) <- paste0("mouse-nuclei_2k_", rownames(sce))



kt <- lapply(k_vals, function(k){
    kc <- as.character(k)
    ifn <- paste0(dp, "filt/meta.data.", kc, ".txt")
    mdk <- read.table(ifn, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    rownames(mdk) <- paste0("mouse-nuclei_2k_", rownames(mdk))
    mdk[,"pct.mt"] <- sce[rownames(mdk),"pct.mt"]
    mdk[,"MALAT1"] <- sce[rownames(mdk),"MALAT1"]
    keep <- intersect(rownames(sf), rownames(mdk))
    mdk[keep,"SpliceFrctn"] <- sf[keep,"SpliceFrctn"]
    mdk <- mdk[mdk$n_genes >= 200, ]
    return(mdk)
})
names(kt) <- k_vals

summ <- sapply(kt, function(datk){
    cl <- datk$Call == "Clean"
    mt <- mean(datk[cl,"pct.mt"])
    malat <- mean(datk[cl,"MALAT1"])
    ngene <- mean(datk[cl,"n_genes"])
    sf <- mean(datk[cl,"SpliceFrctn"])
    nclust <- length(unique(datk[cl,"Cluster"]))
    npass <- sum(cl)
    ntotal <- nrow(datk)
    nfilter <- 100 * npass / ntotal
    return(c(npass, nclust, nfilter, mt, malat, ngene, sf))
})

summ <- t(summ)
summ <- as.data.frame(summ)
colnames(summ) <- c("NumPass", "Final K", "PercentPass", "AvgMT%", "AvgMALAT1%", "AvgNumGenes", "SpliceFraction")
rownames(summ) <- k_vals
summ[,"K"] <- k_vals

summ_mb <- summ[,c("NumPass", "Final K", "PercentPass", "SpliceFraction", "K")]

summl_mb <- reshape2::melt(summ_mb, id.vars = "K")

p3 <- ggplot(summl_mb, aes(x = K, y = value)) +
geom_line() +
facet_wrap(~variable, scales="free", nrow=1) +
ylab("") +
xlab("Filter") +
theme_minimal() +
ggtitle(expression(paste("Mouse Brain ( ", K[init], " = 30)"))) +
theme(plot.title = element_text(hjust = 0.5))


#=========================================================
# Plots
#=========================================================

dplot <- "results/plots/"
pdfname <- paste0(dplot, "param.values.pdf")
jpgname <- paste0(dplot, "param.values.jpeg")
pdf(pdfname, width=8,height=11)
ggarrange(p2, p3, p1, p4, nrow = 4, labels=c("a", "b", "c", "d"),
        font.label = list(size = 18, face="bold"))
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))



