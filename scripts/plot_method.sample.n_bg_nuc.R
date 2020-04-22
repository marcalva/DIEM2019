
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


#=========================================
# Plot number passing droplets along with
# number of nuclear droplets passing
#=========================================

datf <- subset(all_sf, Call == "Clean")

datf_g <- dplyr::group_by(datf, Method, Label)
pct_debris <- dplyr::summarise(datf_g,
    n = sum(Truth == "Nuclear"),
    b = sum(Truth == "Background"))

pct_debris_m <- reshape2::melt(pct_debris, id.vars = c("Method", "Label"))
pct_debris_m$variable = factor(pct_debris_m$variable, levels = c("b", "n"))

p <- ggplot(pct_debris_m, aes(x = Method, y = value)) +
geom_bar(aes(fill = variable), stat = "identity", position=position_stack()) +
scale_fill_discrete(labels = c("Background", "Nuclear")) +
facet_wrap(~Label, nrow = 1, scales = "free") +
ylab("Droplets\npassing") + 
theme_bw() +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      legend.position = "top", 
      legend.title = element_blank(), 
      legend.text=element_text(size=24, 
                               margin = margin(r = 12, unit = "pt")),
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      axis.title = element_text(size = 20),
      axis.title.y = element_text(size = 22),
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "method.sample.n_droplets.clean.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "method.sample.n_droplets.clean.pdf")
ggsave(outfn, width = 20, height = 4)

#=========================================
# Plot number removed droplets along with
# number of nuclear droplets passing
#=========================================

datf <- subset(all_sf, Call == "Debris")

datf_g <- dplyr::group_by(datf, Method, Label)
pct_debris <- dplyr::summarise(datf_g,
    n = sum(Truth == "Nuclear"),
    b = sum(Truth == "Background"))

pct_debris_m <- reshape2::melt(pct_debris, id.vars = c("Method", "Label"))
pct_debris_m$variable = factor(pct_debris_m$variable, levels = c("b", "n"))

p <- ggplot(pct_debris_m, aes(x = Method, y = value)) +
geom_bar(aes(fill = variable), stat = "identity", position=position_stack()) +
scale_fill_discrete(labels = c("Background", "Nuclear")) +
facet_wrap(~Label, nrow = 1, scales = "free") +
ylab("Droplets\nremoved") +
theme_bw() +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      legend.position = "top", 
      legend.title = element_blank(), 
      legend.text=element_text(size=24, 
                               margin = margin(r = 12, unit = "pt")),
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      axis.title = element_text(size = 20),
      axis.title.y = element_text(size = 22),
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "method.sample.n_droplets.debris.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "method.sample.n_droplets.debris.pdf")
ggsave(outfn, width = 20, height = 4)



# Compare n debris paired in passing
datf_s <- all_sf[all_sf$Call == "Clean",]
at_ndebris <- sapply(at_labs, function(s){
                     datf <- datf_s[datf_s[,"Sample"] == s,]
                     qcl <- datf$Method == "Quantile"
                     ecl <- datf$Method == "EmptyDrops"
                     dcl <- datf$Method == "DIEM"
                     mq <- sum(datf[qcl,"Truth"] == "Background") / sum(qcl)
                     me <- sum(datf[ecl,"Truth"] == "Background") / sum(ecl)
                     md <- sum(datf[dcl,"Truth"] == "Background") / sum(dcl)
                     return(c("EmptyDrops" = me, "DIEM" = md, "Quantile" = mq))
})
at_ndebris <- t(at_ndebris)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"EmptyDrops"], paired=T)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"Quantile"], paired=T)

# Compare n nuclear paired in passing
datf_s <- all_sf[all_sf$Call == "Clean",]
at_ndebris <- sapply(at_labs, function(s){
                     datf <- datf_s[datf_s[,"Sample"] == s,]
                     qcl <- datf$Method == "Quantile"
                     ecl <- datf$Method == "EmptyDrops"
                     dcl <- datf$Method == "DIEM"
                     mq <- sum(datf[qcl,"Truth"] == "Nuclear") / sum(qcl)
                     me <- sum(datf[ecl,"Truth"] == "Nuclear") / sum(ecl)
                     md <- sum(datf[dcl,"Truth"] == "Nuclear") / sum(dcl)
                     return(c("EmptyDrops" = me, "DIEM" = md, "Quantile" = mq))
})
at_ndebris <- t(at_ndebris)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"EmptyDrops"], paired=T)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"Quantile"], paired=T)

# Compare n nuclear paired in removed
datf_s <- all_sf[all_sf$Call == "Debris",]
at_ndebris <- sapply(at_labs, function(s){
                     datf <- datf_s[datf_s[,"Sample"] == s,]
                     qcl <- datf$Method == "Quantile"
                     ecl <- datf$Method == "EmptyDrops"
                     dcl <- datf$Method == "DIEM"
                     mq <- sum(datf[qcl,"Truth"] == "Nuclear") / sum(qcl)
                     me <- sum(datf[ecl,"Truth"] == "Nuclear") / sum(ecl)
                     md <- sum(datf[dcl,"Truth"] == "Nuclear") / sum(dcl)
                     return(c("EmptyDrops" = me, "DIEM" = md, "Quantile" = mq))
})
at_ndebris <- t(at_ndebris)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"EmptyDrops"], paired=T)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"Quantile"], paired=T)

# Compare n nuclear paired in removed
datf_s <- all_sf[all_sf$Call == "Debris",]
at_ndebris <- sapply(at_labs, function(s){
                     datf <- datf_s[datf_s[,"Sample"] == s,]
                     qcl <- datf$Method == "Quantile"
                     ecl <- datf$Method == "EmptyDrops"
                     dcl <- datf$Method == "DIEM"
                     mq <- sum(datf[qcl,"Truth"] == "Background") / sum(qcl)
                     me <- sum(datf[ecl,"Truth"] == "Background") / sum(ecl)
                     md <- sum(datf[dcl,"Truth"] == "Background") / sum(dcl)
                     return(c("EmptyDrops" = me, "DIEM" = md, "Quantile" = mq))
})
at_ndebris <- t(at_ndebris)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"EmptyDrops"], paired=T)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"Quantile"], paired=T)


# Quantifications

# Passing droplets
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


# Removed droplets
datf_s <- all_sf[all_sf$Call == "Debris",]
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

for (s in samples){
    for (m in mthds){
        k <- datf_s$Method == m & datf_s$Sample == s
        tb <- table(datf_s[k,"Truth"])
        tbb <- tb["Background"]
        tbc <- tb["Nuclear"]
        tall <- sum(tb)
        message(m, "\t", s, "\t: all:", tall, "; BG: ", tbb, "(", 100*tbb/tall, ")", 
                "; N: ", tbc, "(", 100*tbc/tall, ")")
    }
}
