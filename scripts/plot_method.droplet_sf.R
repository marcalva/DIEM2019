
# Compare fraction spliced in each of the experiments
# and methods

setwd("../")

library(diem)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
source("scripts/common/plotting.R")

#=========================================
# Plot functions
#=========================================

ct <- theme(legend.position="none",
            text=element_text(size=16),
            plot.title=element_text(size=18, hjust=0.5, face="bold")
            )

yexp <- 1.1

xa_angle <- 30
hj <- 0.9
ts <- 18

#=========================================
#=========================================

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
expr <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- expr
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples
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

all_sf[,"SpliceFrctn"] <- all_sf[,"SpliceFrctn"] / 100

midpoints <- read.table("data/raw/splice_frctn/midpoints.txt", 
                        stringsAsFactors=FALSE)
rownames(midpoints) <- midpoints[,1]
midpoints[,1] <- labels[midpoints[,1]]
colnames(midpoints) <- c("Label", "h")
midpoints[,"Label"] <- factor(midpoints[,"Label"], levels = labels)


#=========================================
# Plot method vs fraction spliced in passing droplets
#=========================================

datf <- subset(all_sf, Call == "Clean")

p <- ggplot(datf, aes(x = Method, y = SpliceFrctn)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(alpha = 0.1, width = 0.1, shape=16) +
facet_wrap(~Label, nrow = 1) +
geom_abline(data = midpoints, aes(intercept = h, slope = 0, col="red")) + 
ylab("Percent spliced reads\nin passingd droplets") +
theme_bw() +
ct +
scale_y_continuous(labels = scales::percent, 
                   breaks = c(0, .25, .5, .75, 1),
                   limits = c(0, 1.25)) +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      axis.title = element_text(size = 20),
      axis.title.y = element_text(size = 22),
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "method.sample.pct_sf_droplet.clean.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "method.sample.pct_sf_droplet.clean.pdf")
ggsave(outfn, width = 20, height = 4)

#=========================================
# Plot method vs fraction spliced in debris droplets
#=========================================

datf <- subset(all_sf, Call == "Debris")

p <- ggplot(datf, aes(x = Method, y = SpliceFrctn)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(alpha = 0.1, width = 0.1, shape=16) +
facet_wrap(~Label, nrow = 1) +
geom_abline(data = midpoints, aes(intercept = h, slope = 0, col="red")) + 
ylab("Percent spliced reads\nin removed droplets") +
scale_y_continuous(labels = scales::percent, 
                   breaks = c(0, .25, .5, .75, 1),
                   limits = c(0, 1.25)) +
theme_bw() +
ct +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      axis.title = element_text(size = 20),
      axis.title.y = element_text(size = 22),
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "method.sample.pct_sf_droplet.debris.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "method.sample.pct_sf_droplet.debris.pdf")
ggsave(outfn, width = 20, height = 4)


#=========================================
# Run stats
#=========================================

sf_col <- "SpliceFrctn"

# Wilcox test of droplets
datf_s <- all_sf[all_sf[,"Call"] == "Clean",]
for (s in samples){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    qcl <- datf$Method == "Quantile"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[qcl,sf_col], datf[dcl,sf_col], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value * 8)
}
for (s in samples){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    ecl <- datf$Method == "EmptyDrops"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[ecl,sf_col], datf[dcl,sf_col], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value) * 8
}

datf_s <- all_sf[all_sf[,"Call"] == "Debris",]
for (s in samples){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    qcl <- datf$Method == "Quantile"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[qcl,sf_col], datf[dcl,sf_col], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value * 8)
}
for (s in samples){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    ecl <- datf$Method == "EmptyDrops"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[ecl,sf_col], datf[dcl,sf_col], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value) * 8
}


# Quantifications
datf_s <- all_sf[all_sf$Call == "Clean",]
s <- "adpcyte"
for (m in mthds){
    k <- datf_s$Method == m & datf_s$Sample == s
    tb <- table(datf_s[k,"Truth"])
    tbb <- tb["Background"]
    tbc <- tb["Nuclear"]
    tall <- sum(tb)
    message(m, "\t: all:", tall, "; BG: ", tbb, "(", 100*tbb/tall, ")", 
            "; N: ", tbc, "(", 100*tbc/tall, ")")
}

s <- "mouse_nuclei_2k"
for (m in mthds){
    k <- datf_s$Method == m & datf_s$Sample == s
    tb <- table(datf_s[k,"Truth"])
    tbb <- tb["Background"]
    tbc <- tb["Nuclear"]
    tall <- sum(tb)
    message(m, "\t: all:", tall, "; BG: ", tbb, "(", 100*tbb/tall, ")", 
            "; N: ", tbc, "(", 100*tbc/tall, ")")
}

e <- "atsn"
for (m in mthds){
    k <- datf_s$Method == m & datf_s$Exprm == e
    tb <- table(datf_s[k,"Truth"])
    tbb <- tb["Background"]
    tbc <- tb["Nuclear"]
    tall <- sum(tb)
    message(m, "\t: all:", tall, "; BG: ", tbb, "(", 100*tbb/tall, ")", 
            "; N: ", tbc, "(", 100*tbc/tall, ")")
}


