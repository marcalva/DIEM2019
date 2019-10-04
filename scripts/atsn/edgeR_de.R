
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


get_lh <- function(sce, cutpoint, cpm_thresh=3, prefix=NULL){
	tc <- Matrix::colSums(sce@counts)

	low_counts <- Matrix::rowSums(sce@counts[,tc < cutpoint])
	high_counts <- Matrix::rowSums(sce@counts[,tc >= cutpoint])

	low_cpm <- 1e6*low_counts/sum(low_counts)
	high_cpm <- 1e6*high_counts/sum(high_counts)

	keep <- (low_cpm >= cpm_thresh) & (high_cpm >= cpm_thresh)
	lh <- data.frame(Low=low_counts[keep], High=high_counts[keep])
	if (!is.null(prefix)){
		colnames(lh) <- paste0(prefix, ".", colnames(lh))
	}
	return(lh)
}

#=========================================
#=========================================

# Make directories

dp <- paste0("data/processed/", labl, "/")
dir.create(dp, recursive=TRUE, showWarnings=FALSE)
dr <- paste0("results/", labl, "/", method, "/")
dir.create(dr, recursive=TRUE, showWarnings=FALSE)

labl <- "atsn"
scel <- lapply(lab_ids, function(id){
	expr <- read_10x(paths[[id]])
	sce <- create_SCE(expr, name=id)
	})
names(scel) <- lab_ids

lh <- lapply(lab_ids, function(id) get_lh(scel[[id]], cutpoint=100, cpm_thresh=0, prefix=id))
counts <- do.call(cbind, lh)

# EdgeR
y = edgeR::DGEList(counts=counts)
keep <- rowSums(cpm(y) > 0) >= ncol(y)/2
y <- y[keep, , keep.lib.sizes=FALSE]
geneNamesKeep = y$genes
y <- calcNormFactors(y)

# Calculate design
pair <- factor(sapply(lab_ids, function(id) rep(id, 2)))
group <- factor(rep(c("low", "high"), length(lab_ids)))
design <- model.matrix(~pair + group)
rownames(design) <- rownames(y$samples)
#y$samples$group <- design[,ncol(design)]

y <- estimateDisp(y, design=design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
diffTable <- qlf$table
# et <- exactTest(y)
# diffTable <- et$table
diffTable$FDR <- p.adjust(diffTable$PValue, "fdr")
diffTable$p_adj <- p.adjust(diffTable$PValue, "bonferroni")
diffTable <- diffTable[order(diffTable$FDR),]

print(table(diffTable$p_adj < .05))

write.table(diffTable, paste0(dr, "edgeR.paired_de.diff.txt"),
	row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

saveRDS(y, paste0(dr, "edgeR.out.rds"))
# 6,562 genes FDR < 0.01 after testing 13,565 genes

counts_norm <- cpm(y)

#=========================================
# Plots
#=========================================

dt <- diffTable[,c("logFC", "p_adj")]
dt$DE <- ""
dt[dt$p_adj < .05,"DE"] <- "DE"
dt$FDR <- -10 * log10(dt$p_adj)

# Volcano plot

pdf(paste0(dr, "de.volcano.pdf"), width=5, height=4)
p <- ggplot(dt, aes(x=logFC, y=FDR, color=DE)) + 
geom_point(shape=16) + 
theme_minimal() + 
scale_colour_brewer(palette="Set1", breaks="DE", labels="Bonferroni\np < 0.05", name="") + 
ylab(expression(-log[10](p)))

# Annotate
xlim <- c(min(dt$logFC), max(dt$logFC))
ylim <- c(min(dt$FDR), max(dt$FDR))
p <- p + annotate("text", xlim[2], ylim[2], label="Debris", hjust=1, vjust=1)
p <- p + annotate("text", xlim[1], ylim[2], label="Nuclear", hjust=0, vjust=1)
print(p)
dev.off()



