
setwd("../../")

library(diem)
library(future)
library(future.apply)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem"
dir10X <- c("data/raw/AT1/", 
            "data/raw/AT2/",
            "data/raw/AT3/",
            "data/raw/AT4/",
            "data/raw/AT5/", 
            "data/raw/AT6/")

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

args = commandArgs(TRUE)
i = as.integer(args[1])


dr <- paste0("results/", label, "/", method, "/")


testk <- lapply(1:6, function(i){
                fn <- paste0(dr, "test_k.", lab_ids[i], ".txt")
                read.table(fn, header = TRUE, row.names = 1, sep = "\t")
            })


testk <- lapply(1:6, function(i){
                datf <- testk[[i]]
                datf[,"Subject"] <- lab_ids[i]
                datf[,"K"] <- rownames(datf)
                return(datf)
            })

testk_all <- do.call(rbind, testk)
testk_all$K <- as.numeric(testk_all$K)

dplot <- paste0(dr, "plots/")
library(ggplot2)
pdf(paste0(dplot, "test_k.pdf"))
ggplot(testk_all, aes(x = K, y = Pass)) + 
geom_line() + 
facet_wrap( ~ Subject, scales = "free") + 
theme_bw()
dev.off()

sce <- create_SCE(counts, name=lab_ids[i])

mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")

sce@pp_thresh <- 0.95
sce <- set_debris_test_set(sce, top_n = 1e4, min_counts = 100, min_genes = 100)
sce <- filter_genes(sce, cpm_thresh = 0)

sce <- test_k(sce, K = 1:15)

dat <- sce@test_k

dr <- paste0("results/", label, "/", method, "/")
dir.create(dr, recursive=TRUE, showWarnings=FALSE)

fn <- paste0(dr, "test_k.", lab_ids[i], ".txt")
write.table(dat, 
            fn, 
            row.names = TRUE, 
            col.names = NA, 
            quote = FALSE, 
            sep = "\t")

