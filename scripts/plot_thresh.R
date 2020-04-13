
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

dir_plot <- "results/plots/thresh/"
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

p_vals <- seq(0, 1, .05)

datfl <- list()
for (i in 1:8){
    e <- samples[i]
    k <- all_sf$n_genes >= 200 & all_sf$Sample == e & all_sf$Method == "DIEM"
    sfs <- all_sf[k,]
    ret <- as.data.frame(matrix(nrow = length(p_vals), ncol = 6))
    rownames(ret) <- p_vals
    colnames(ret) <- c("Label", "Threshold", "N", "SF", "Nuc", "BG")
    for (p in p_vals){
        pc <- as.character(p)
        k <- sfs[,"score.debris"] < p
        sfss <- sfs[k,]
        l <- sfs[,"score.debris"] >= p
        sfsbg <- sfs[l,,drop=FALSE]
        ret[pc,1] <- labels[i]
        ret[pc,2] <- p
        ret[pc,3] <- nrow(sfss)
        ret[pc,4] <- mean(sfss[,"SpliceFrctn"], na.rm = TRUE)
        ret[pc,5] <- sum(sfss[,"Truth"] == "Nuclear") / nrow(sfss)
        ret[pc,6] <- sum(sfsbg[,"Truth"] == "Background") / nrow(sfsbg)
    }
    datfl[[e]] <- ret
}

datf <- do.call(rbind, datfl)
datf[,"Label"] <- factor(datf[,"Label"], levels = labels)
datf[,"SF"] <- datf[,"SF"] / 100

w = 3
h = 10

# plot number of droplets

p1 <- bp(datf, x = "Threshold", y = "N", fac = "Label",
         ylab = "Number of droplets passing")


fn <- paste0(dir_plot, "thresh.num_pass.exprm.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "thresh.num_pass.exprm.jpeg")
ggsave(fn, width = w, height = h)

# Average splice fraction
p2 <- bp(datf, x = "Threshold", y = "SF", fac = "Label",
         ylab = "Average percent reads spliced")
p2 <- p2 + scale_y_continuous(labels = scales::percent_format(accuracy=1))

fn <- paste0(dir_plot, "thresh.sf.exprm.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "thresh.sf.exprm.jpeg")
ggsave(fn, width = w, height = h)

# Proprtion of passing droplets that are nuclear
p3 <- bp(datf, x = "Threshold", y = "Nuc", fac = "Label",
         ylab = "Percent of passing droplets that are nuclear")
p3 <- p3 + scale_y_continuous(labels = scales::percent_format(accuracy=1)) 

fn <- paste0(dir_plot, "thresh.sf.pass_nuc.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "thresh.sf.pass_nuc.jpeg")
ggsave(fn, width = w, height = h)

# Proprtion of removed droplets that are nuclear
p3 <- bp(datf, x = "Threshold", y = "BG", fac = "Label",
         ylab = "Percent of removed droplets that are background")
p3 <- p3 + scale_y_continuous(labels = scales::percent_format(accuracy=1)) 

fn <- paste0(dir_plot, "thresh.sf.rm_nuc.pdf")
ggsave(fn, width = w, height = h)
fn <- paste0(dir_plot, "thresh.sf.rm_nuc.jpeg")
ggsave(fn, width = w, height = h)

