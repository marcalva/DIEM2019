
setwd("../")

library(diem)
library(ggplot2)
library(ggpubr)


set_breaks_10 <- function(x){
    xmax <- x[2]
    bk <- 10
    brks <- c(bk)
    while (bk < xmax){
        bk <- bk * 10
        brks <- c(brks, bk)
    }
    return(brks)
}

#====================================================
# Adipose
#====================================================


dp <- "data/processed/atsn/diem/"
dr <- "results/atsn/diem/"
dir.create(dr, showWarnings=FALSE, recursive=TRUE)

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

scel <- lapply(lab_ids, function(i){
    scefn <- paste0(dp, i, ".diem_sce.rds")
    sce <- readRDS(scefn)
    mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
    genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")
    md <- sce@droplet_data[,]
    md[md[,"pct.mt"] > 5, "pct.mt"] <- NA
    rownames(md) <- paste0(i, "_", rownames(md))
    keep <- intersect(rownames(md), rownames(sf))
    md[keep,"SpliceFrctn"] <- sf[keep,1]
    md <- md[keep,]
    md$DataSet <- i
    return(md[md$total_counts > 100,])
})
names(scel) <- lab_ids

scel_at <- scel


dp <- "data/processed/adpcyte/diem/"
dr <- "results/adpcyte/diem/"
dir.create(dr, showWarnings=FALSE, recursive=TRUE)

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)
sf <- as.matrix(sf)


scefn <- paste0(dp, "adpcyte.diem_sce.rds")
sce <- readRDS(scefn)
mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")
sce <- sce@droplet_data
rownames(sce) <- paste0("adpcyte_", rownames(sce))
keep <- intersect(rownames(sce), rownames(sf))
sce[keep,"SpliceFrctn"] <- sf[keep,1]
sce <- sce[keep,]
sce_ad <- sce[sce$total_counts > 100,]
sce_ad$DataSet <- "DiffPA"


dp <- "data/processed/mouse_nuclei_2k/diem/"
dr <- "results/mouse_nuclei_2k/diem/"
dir.create(dr, showWarnings=FALSE, recursive=TRUE)

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

scefn <- paste0(dp, "mouse_nuclei_2k.diem_sce.rds")
sce <- readRDS(scefn)
mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")
sce <- sce@droplet_data
rownames(sce) <- paste0("mouse-nuclei_2k_", rownames(sce))
keep <- intersect(rownames(sce), rownames(sf))
sce[keep,"SpliceFrctn"] <- sf[keep,1]
sce <- sce[keep,]
sce_mb <- sce[sce$total_counts > 100,]
sce_mb$DataSet <- "Mouse Brain"
sce_mb[sce_mb[,"MALAT1"] > 20,"MALAT1"] <- NA

sce_all <- c(scel_at, list(sce_mb, sce_ad))

datf <- do.call(rbind, sce_all)
datf$SpliceFrctn <- 100 * datf$SpliceFrctn


library(ggplot2)

p_tc <- p <- ggplot(datf, aes(x = SpliceFrctn, y = total_counts)) +
geom_point(shape=16, alpha = 0.01) +
facet_wrap(~DataSet, ncol=1) + 
xlab("Percent reads spliced") +
ylab("Total counts") +
scale_y_continuous(trans='log10', breaks=set_breaks_10, labels=scales::comma) +
theme_minimal() + theme(text=element_text(size=22))

p_malat <- ggplot(datf, aes(x = SpliceFrctn, y = MALAT1)) + 
geom_point(shape=16, alpha = 0.01) +  
facet_wrap(~DataSet, ncol=1, scales = "free_y") +
xlab("Percent reads spliced") + 
ylab("MALAT1%") + 
theme_minimal() + theme(text=element_text(size=22))


p_mt <- ggplot(datf, aes(x = SpliceFrctn, y = pct.mt)) + 
geom_point(shape=16, alpha = 0.01) +  
facet_wrap(~DataSet, ncol=1, scales = "free_y") + 
xlab("Percent reads spliced") + 
ylab("MT%") + 
theme_minimal() + theme(text=element_text(size=22))

p_dens <- ggplot(datf, aes(x = SpliceFrctn)) + 
facet_wrap(~DataSet, ncol=1, scales = "free_y") +
geom_density() + 
ylab("Density") + 
xlab("Percent reads spliced") + 
theme_minimal() + theme(text=element_text(size=22))

library(ggpubr)

# Arrange into plot
fl <- list(size = 20, face="bold")
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "sf_relat.pdf")
jpgname <- paste0(dir_plot, "sf_relat.jpeg")
pdf(pdfname, width=18, height=14)
ggarrange(p_tc, p_dens, p_mt, p_malat, labels=c("a", "b", "c", "d"), font.label = fl,ncol=4)
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


