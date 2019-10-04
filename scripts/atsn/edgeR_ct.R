
setwd("../../")

library(diem)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(Matrix)

labl <- "atsn"
method <- "DE"

paths <- c("data/raw/AT1/",
           "data/raw/AT2/",
           "data/raw/AT3/",
           "data/raw/AT4/",
           "data/raw/AT5/",
           "data/raw/AT6/")

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
names(paths) <- lab_ids

seur <- readRDS("data/processed/atsn/quantile/atsn.seur_obj.rds")

dr <- paste0("results/", labl, "/", method, "/")
dir.create(dr, recursive=TRUE, showWarnings=FALSE)

labl <- "atsn"
scel <- lapply(lab_ids, function(id){
    expr <- read_10x(paths[[id]])
    sce <- create_SCE(expr, name=id)
    })
names(scel) <- lab_ids


#=============================================================
# Correlation scatter plot of TMM-normalized CPM
#=============================================================
# Not very representative of the differences
# I deleted what was below

#=============================================================
# DE
#=============================================================

set.seed(1)
low_counts <- sapply(scel, function(sce) {
                     tc <- Matrix::colSums(sce@counts)
                     i <- which(tc < 100)
                     i <- sample(i, 9e3)
                     Matrix::rowSums(sce@counts[,i])
})

get_ncounts <- function(s, ct){
    keep <- s$seurat_clusters %in% ct
    return(sum(Matrix::rowSums(s@assays$RNA@counts[,keep, drop=FALSE])))
}

# Function to get DE between cell types
get_de <- function(s, ct1, ct2){
    subj <- sort(unique(s$orig.ident))
    ct1_counts <- sapply(subj, function(i){
                         keep <- s$orig.ident == i & s$seurat_clusters == ct1
                         Matrix::rowSums(s@assays$RNA@counts[,keep, drop=FALSE])
                         })
    colnames(ct1_counts) <- paste0(colnames(ct1_counts), "_CT1")
    ct2_counts <- sapply(subj, function(i){
                         keep <- s$orig.ident == i & s$seurat_clusters %in% ct2
                         Matrix::rowSums(s@assays$RNA@counts[,keep, drop=FALSE])
                         })
    colnames(ct2_counts) <- paste0(colnames(ct2_counts), "_CT2")
    counts <- cbind(ct1_counts, ct2_counts)

    pair <- c(subj, subj)
    group <- factor(c(rep("CT1", length(subj)), rep("CT2", length(subj))))
    design <- model.matrix(~pair + group)
    rownames(design) <- colnames(counts)
    
    y <- edgeR::DGEList(counts=counts)
    keep <- rowSums(cpm(y) > 0) >= ncol(y)/2
    y <- y[keep, , keep.lib.sizes=FALSE]
    geneNamesKeep = y$genes
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design=design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit)
    diffTable <- qlf$table
    diffTable <- diffTable[order(diffTable$PValue),]
    diffTable$p_adj <- p.adjust(diffTable$PValue, "bonferroni")
    return(list("percent"=sum(diffTable$p_adj < .05)/nrow(diffTable), 
                "n_genes"=nrow(diffTable)))
}

# Function to get DE between cell types and debris
get_de_debris <- function(s, ct1){
    subj <- sort(unique(s$orig.ident))
    ct1_counts <- sapply(subj, function(i){
                         keep <- s$orig.ident == i & s$seurat_clusters == ct1
                         Matrix::rowSums(s@assays$RNA@counts[,keep, drop=FALSE])
                         })
    colnames(ct1_counts) <- paste0(colnames(ct1_counts), "_CT1")
    ct2_counts <- low_counts
    colnames(ct2_counts) <- paste0(colnames(ct2_counts), "_CT2")
    counts <- cbind(ct1_counts, ct2_counts)

    pair <- c(subj, subj)
    group <- factor(c(rep("CT1", length(subj)), rep("CT2", length(subj))))
    design <- model.matrix(~pair + group)
    rownames(design) <- colnames(counts)
    
    y <- edgeR::DGEList(counts=counts)
    keep <- rowSums(cpm(y) > 0) >= ncol(y)/2
    y <- y[keep, , keep.lib.sizes=FALSE]
    geneNamesKeep = y$genes
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design=design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit)
    diffTable <- qlf$table
    diffTable <- diffTable[order(diffTable$PValue),]
    diffTable$p_adj <- p.adjust(diffTable$PValue, "bonferroni")
    return(list("percent"=sum(diffTable$p_adj < .05)/nrow(diffTable), 
                "n_genes"=nrow(diffTable)))
}

cell_types <- levels(seur@active.ident)
nct <- table(seur$orig.ident, seur$seurat_clusters)
cell_types <- colnames(nct)[apply(nct, 2, function(x) all(x > 0))]

#n_de <- matrix(NA, nrow=length(cell_types)+1, ncol=length(cell_types)+1)
n_de <- matrix(NA, nrow=length(cell_types), ncol=2)
n_genes <- n_de

# Get DE between the cell types
for (i in 1:(length(cell_types))){
    j <- setdiff(cell_types, i)
    ret <- get_de(seur, cell_types[i], j)
    n_de[i,1] <- ret$percent
    n_genes[i,1] <- ret$n_genes
    #for (j in (i+1):length(cell_types)){
    #    ret <- get_de(seur, cell_types[i], cell_types[j])
    #    n_de[i,j] <- ret$percent
    #    n_genes[i,j] <- ret$n_genes
    #}
}

# Get DE between debris and cell types
for (i in 1:length(cell_types)){
    ret <- get_de_debris(seur, cell_types[i])
    n_de[i,ncol(n_de)] <- ret$percent
    n_genes[i,ncol(n_de)] <- ret$n_genes
}

#rownames(n_de) <- c(cell_types, "Debris")
#rownames(n_genes) <- c(cell_types, "Debris")
rownames(n_de) <- rownames(n_genes) <- cell_types
colnames(n_de) <- colnames(n_genes) <- c("Other_CellTypes", "Debris")

# Write out results
write.table(n_de, paste0(dr, "edgeR.n_de.cell_types.txt"),
            row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

write.table(n_genes, paste0(dr, "edgeR.n_genes.cell_types.txt"),
            row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

# Get average depth between 2 paris
#n_counts <- matrix(NA, nrow=length(cell_types)+1, ncol=length(cell_types)+1)
n_counts <- matrix(NA, nrow=length(cell_types), ncol=2)
for (i in 1:(length(cell_types))){
    j <- setdiff(cell_types, i)
    nc1 <- get_ncounts(seur, cell_types[i])
    nc2 <- get_ncounts(seur, j)
    n_counts[i,1] <- (nc1 + nc2) / 2
    #for (j in (i+1):length(cell_types)){
    #    nc1 <- get_ncounts(seur, cell_types[i])
    #    nc2 <- get_ncounts(seur, cell_types[j])
    #    n_counts[i,j] <- (nc1 + nc2) / 2
    #}
}

nclow <- sum(colSums(low_counts))
for (i in 1:(length(cell_types))){
    n_ct1 <- get_ncounts(seur, cell_types[i])
    n_counts[i,ncol(n_counts)] <- (n_ct1 + nclow) / 2
}

rownames(n_counts) <- cell_types
colnames(n_counts) <- c("Other_CellTypes", "Debris")

write.table(n_counts, paste0(dr, "edgeR.n_counts.cell_types.txt"), 
            row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

