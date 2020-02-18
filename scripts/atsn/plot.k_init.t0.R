# Show effect of increasing k_init keeping filt = 0

setwd("../../")

library(diem)

dp <- "data/processed/atsn/diem/"
dr <- "results/atsn/diem/k_init.t0/"
dir.create(dr, showWarnings=FALSE, recursive=TRUE)

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
k_vals <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

scel <- lapply(lab_ids, function(i){
    scefn <- paste0(dp, i, ".diem_sce.rds")
    sce <- readRDS(scefn)
    mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
    genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")
    md <- sce@droplet_data
    return(md)
})
names(scel) <- lab_ids

kt <- lapply(k_vals, function(k){
    kc <- as.character(k)
    mdkl <- lapply(lab_ids, function(i){
        ifn <- paste0(dp, "k_init.t0/meta.data.", i, ".", kc, ".txt")
        mdk <- read.table(ifn, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
        mdk[,"pct.mt"] <- scel[[i]][rownames(mdk),"pct.mt"]
        mdk[,"MALAT1"] <- scel[[i]][rownames(mdk),"MALAT1"]
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
    npass <- sum(cl)
    ntotal <- nrow(datk)
    nfilter <- 100 * npass / ntotal
    return(c(npass, nfilter, mt, malat, ngene))
})

summ <- t(summ)
summ <- as.data.frame(summ)
colnames(summ) <- c("NumPass", "PercentPass", "AvgMT%", "AvgMALAT1%", "AvgNumGenes")
rownames(summ) <- k_vals
summ[,"K"] <- k_vals

summ <- summ[,c("NumPass", "PercentPass", "AvgMT%", "K")]

summl <- reshape2::melt(summ, id.vars = "K")

#=========================================================
# Plots
#=========================================================

library(ggplot2)
p1 <- ggplot(summl, aes(x = K, y = value)) +
geom_line() + 
facet_wrap(~variable, scales="free", nrow=1) + 
ylab("") + 
xlab(expression(k[init])) +
theme_minimal() + 
ggtitle("Adipose Tissue") + 
theme(plot.title = element_text(hjust = 0.5))


pdfname <- paste0(dr, "k_init.t0.pdf")
jpgname <- paste0(dr, "k_init.t0.jpeg")
pdf(pdfname, width=10,height=3)
print(p1)
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))


