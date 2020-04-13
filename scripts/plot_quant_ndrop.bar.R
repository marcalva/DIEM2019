
setwd("../")

library(diem)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
source("scripts/common/plotting.R")

#=========================================
#=========================================

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
expr <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- expr
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples
eLabels <- c("DiffPA", "Mouse Brain", rep("Adipose Tissue", 6))
names(eLabels) <- samples

mthds <- c("DIEM", "EmptyDrops", "Quantile")

#=========================================
# Read in data
#=========================================

all_sf <- read.table("data/processed/meta_data.txt",
                     header=TRUE,
                     sep="\t")

all_sf <- all_sf[all_sf[,"n_genes"] >= 200,]

all_sf[,"Exprm"] <- factor(all_sf[,"Exprm"], levels=unique(expr))
all_sf[,"Sample"] <- factor(all_sf[,"Sample"], levels=samples)
all_sf[,"Label"] <- factor(all_sf[,"Label"], levels=labels)
all_sf[,"eLabel"] <- factor(eLabels[all_sf[,"Sample"]], levels=unique(eLabels))

mapcall <- c("Above", "Below")
names(mapcall) <- c("Clean", "Debris")
all_sf[,"Call"] <- factor(all_sf[,"Call"], 
                          labels = mapcall[levels(all_sf[,"Call"])])
all_sf[,"Call"] <- factor(all_sf[,"Call"], levels <- rev(mapcall))

#=========================================
# Plot barplot
#=========================================

datf <- subset(all_sf, Method == "Quantile")

datf_g <- dplyr::group_by(datf, Call, Exprm, eLabel)
pct_debris <- dplyr::summarise(datf_g,
    Nuclear = sum(Truth == "Nuclear"),
    Background = sum(Truth == "Background"))

pct_debris_m <- reshape2::melt(pct_debris, id.vars = c("Call", "Exprm", "eLabel"))
pct_debris_m[,"variable"] <- factor(pct_debris_m[,"variable"], 
                    levels = c("Background", "Nuclear"))
#pct_debris_m$variable = factor(pct_debris_m$variable, levels = c("b", "n"))

p <- ggplot(pct_debris_m, aes(x = Call, y = value)) +
geom_bar(aes(fill = variable), stat = "identity", position=position_stack()) +
facet_wrap(~eLabel, ncol = 1, scales = "free") +
ylab("Droplets\npassing") +
theme_bw() +
theme(strip.background = element_blank(),
      strip.text = element_text(size = 14),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(size = 14))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "quantile.call.n_drop.barplot.jpeg")
ggsave(outfn, width = 3, height = 6)
outfn <- paste0(dir_plot, "quantile.call.n_drop.barplot.pdf")
ggsave(outfn, width = 3, height = 6)


