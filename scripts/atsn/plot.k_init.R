
# Plot over different values of k_init for filter = 0.25

setwd("../../")

dp <- "data/processed/atsn/diem/k_init/"

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
k_vals <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

ifn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(ifn, header = TRUE, row.names = 1)

scel <- lapply(lab_ids, function(i){
    scefn <- paste0("data/processed/atsn/diem/", i, ".diem_sce.rds")
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
        fn <- paste0("meta.data.", i, k, ".txt")
        ifn <- file.path(dp, fn)
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
    sf <- mean(datk[cl,"SpliceFrctn"])
    nclust <- mean(tapply(datk[cl,"Cluster"], datk[cl, "ID"], function(i){length(unique(i))}))
    npass <- sum(cl)
    ntotal <- nrow(datk)
    nfilter <- 100 * npass / ntotal
    return(c(npass, nclust, nfilter, mt, malat, ngene))
})

summ <- t(summ)
summ <- as.data.frame(summ)
colnames(summ) <- c("NumPass", "Final K", "PercentPass", "AvgMT%", "AvgMALAT1%", "AvgNumGenes")
rownames(summ) <- k_vals
summ[,"K"] <- k_vals


summ_at <- summ[,c("NumPass", "Final K", "PercentPass", "AvgMT%", "K")]

summl_at <- reshape2::melt(summ_at, id.vars = "K")

p1 <- ggplot(summl_at, aes(x = K, y = value)) +
geom_line() +
facet_wrap(~variable, scales="free", nrow=1) +
ylab("") +
xlab("Filter") +
theme_minimal() +
ggtitle("Adipose Tissue") +
theme(plot.title = element_text(hjust = 0.5))

dir_o <- "results/atsn/diem/k_init/"
dir.create(dir_o, showWarnings = FALSE, recursive = TRUE)
pdfname <- paste0(dir_o, "k_init.filt25.pdf")
jpgname <- paste0(dir_o, "k_init.filt25.jpeg")
pdf(pdfname, width=6,height=3)
p1
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))


