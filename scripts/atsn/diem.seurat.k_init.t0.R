
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

# Plot seurat UMAP
ct <- theme(text=element_text(size=16),
            plot.title=element_text(size=18, hjust=0.5, face="bold")
            )


yexp <- 1.1

plot_umap_fct <- function(x, names, colname="percent.mt", legend_name="MT%",
                          color_limits=NULL, color_breaks=waiver(),
                          size=1, alpha=alpha, order_col = TRUE){
    dfl <- lapply(1:length(x), function(i) {
                  data.frame(Mito=x[[i]]@meta.data[,colname],
                             UMAP1=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_1"],
                             UMAP2=x[[i]]@reductions$umap@cell.embeddings[,"UMAP_2"], Method=names[[i]])
                          })
    df <- do.call(rbind, dfl)
    if (order_col) df <- df[order(df$Mito, decreasing=FALSE),,drop=FALSE]

    p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Mito)) +
    geom_point(size=size, alpha=0.8) + theme_bw() +
    facet_wrap(~Method, ncol=3, scale="free") +
    theme(text=element_text(size=16),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title=element_text(hjust=0.5),
          panel.grid=element_blank()) +
    scale_colour_gradient(low="gray90", high="red3",
                          name=legend_name,
                          limits=color_limits,
                          breaks=color_breaks)
    #scale_color_distiller(palette = "Spectral", name=legend_name, 
    #                     limits=color_limits, breaks=color_breaks) 

    return(p)
}

seur_k1 <- readRDS("data/processed/atsn/diem/k_init.t0/seur.k1.rds")
seur_k30 <- readRDS("data/processed/atsn/diem/k_init.t0/seur.k30.rds")

ch <- seur_k1@meta.data$percent.mt > 10
seur_k1@meta.data$percent.mt[ch] <- 10
ch <- seur_k30@meta.data$percent.mt > 10
seur_k30@meta.data$percent.mt[ch] <- 10

ltitles <- c("Initial K = 1", 
             "Initial K = 30")
p2 <- plot_umap_fct(list(seur_k1, seur_k30), ltitles) + ggtitle("Adipose Tissue")

pdfname <- paste0(dr, "k_init.umap_mt.pdf")
jpgname <- paste0(dr, "k_init.umap_mt.jpeg")
pdf(pdfname, width=12,height=5)
print(p2)
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))



