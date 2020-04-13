
setwd("../")

library(diem)
library(ggplot2)
library(scales)
library(reshape2)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# scatter plot
#=========================================

bp <- function(datf, x, y, fac, ncol = 1, ylab){
    p <- ggplot(datf, aes_string(x = x, y = y)) +
    geom_point(shape=3) +
    facet_wrap(c(fac), ncol = ncol, scales = "free") +
    theme_minimal() +
    theme(axis.title = element_text(size = 18)) +
    ylab(ylab)
    return(p)
}


#=========================================
# Set variables
#=========================================

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
exprm <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- exprm
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples
eLabels <- c("DiffPA", "Mouse Brain", rep("Adipose Tissue", 6))
names(eLabels) <- samples

mthds <- c("DIEM", "EmptyDrops", "Quantile")

#=========================================
#=========================================

dir_plot <- "results/plots/k/"
dir.create(dir_plot, showWarnings = FALSE)

#exprm <- c("adpcyte", "mouse_nuclei_2k", "atsn")
#eLabels <- c("DiffPA", "Mouse Brain", "Adipose Tissue")

all_sf <- read.table("data/processed/meta_data.txt",
                     header=TRUE,
                     sep="\t")

all_sf <- all_sf[all_sf[,"n_genes"] >= 200,]

all_sf[,"Exprm"] <- factor(all_sf[,"Exprm"], levels=unique(exprm))
all_sf[,"Sample"] <- factor(all_sf[,"Sample"], levels=samples)
all_sf[,"Label"] <- factor(all_sf[,"Label"], levels=labels)

k_vals <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
p <- 0.5

datfl <- list()

for (i in at_labs){
    ret <- datf <- as.data.frame(matrix(nrow = length(k_vals), ncol = 6))
    rownames(ret) <- k_vals
    colnames(ret) <- c("Label", "k", "N", "SF", "Nuc", "BG")

    for (k in k_vals){
        kc <- as.character(k)
        fn <- paste0("data/processed/atsn/diem.k/k", kc, "/", i, ".diem_sce.rds")
        sce <- readRDS(fn)
        test_dat <- sce@test_data
        nona <- ! is.na(test_dat[,"SpliceFrctn"])
        test_dat <- test_dat[nona,]
        keep <- test_dat[,"score.debris"] < p
        test_dat_s <- test_dat[keep,]
        l <- test_dat[,"score.debris"] >= p
        test_dat_bg <- test_dat[l,,drop=FALSE]
        
        ret[kc,1] <- i
        ret[kc,2] <- kc
        ret[kc,3] <- nrow(test_dat_s)
        ret[kc,4] <- mean(test_dat_s[,"SpliceFrctn"], na.rm = TRUE)
        ret[kc,5] <- sum(test_dat_s[,"Truth"] == "Nuclear") / nrow(test_dat_s)
        ret[kc,6] <- sum(test_dat_bg[,"Truth"] == "Nuclear") / nrow(test_dat_bg)
    }
    datfl[[i]] <- ret
}

datf <- do.call(rbind, datfl)
datf[,"Label"] <- factor(datf[,"Label"], levels = at_labs)
datf[,"SF"] <- datf[,"SF"] / 100
datf[,"k"] <- as.numeric(datf[,"k"])


w = 3
h = 8

# plot number of droplets

p1 <- bp(datf, x = "k", y = "N", fac = "Label",
         ylab = "Number of droplets passing")

fn <- paste0(dir_plot, "k.num_pass.exprm.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "k.num_pass.exprm.jpeg")
ggsave(fn, width = w, height = h)

# Average splice fraction
p2 <- bp(datf, x = "k", y = "SF", fac = "Label",
         ylab = "Average percent reads spliced")
p2 <- p2 + scale_y_continuous(labels = scales::percent_format(accuracy=.1))

fn <- paste0(dir_plot, "k.sf.exprm.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "k.sf.exprm.jpeg")
ggsave(fn, width = w, height = h)

# Proprtion of passing droplets that are nuclear
p3 <- bp(datf, x = "k", y = "Nuc", fac = "Label",
         ylab = "Percent of passing droplets that are nuclear")
p3 <- p3 + scale_y_continuous(labels = scales::percent_format(accuracy=.1))

fn <- paste0(dir_plot, "k.sf.pass_nuc.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "k.sf.pass_nuc.jpeg")
ggsave(fn, width = w, height = h)

# Proprtion of removed droplets that are nuclear
p3 <- bp(datf, x = "k", y = "BG", fac = "Label",
         ylab = "Percent of removed droplets that are nuclear")
p3 <- p3 + scale_y_continuous(labels = scales::percent_format(accuracy=.1))

fn <- paste0(dir_plot, "k.sf.rm_nuc.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "k.sf.rm_nuc.jpeg")
ggsave(fn, width = w, height = h)


