
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

#=========================================
# Read in data
#=========================================

infn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(infn, header = TRUE, row.names = 1)

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
expr <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- expr
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples

# Droplet prefixes
dpre <- c("adpcyte", "mouse-nuclei_2k", at_labs)
names(dpre) <- samples

sf_calls <- list()

for (i in 1:length(expr)){
    s <- samples[i]
    e <- names(samples)[i]

    # Read in results
    diem_fn <- paste0("data/processed/", e, "/diem/", s, ".diem_sce.rds")
    e_diem <- readRDS(diem_fn)
    ed_fn <- paste0("data/processed/", e, "/emptydrops/", s, ".emptydrops_counts.rds")
    e_ed <- readRDS(ed_fn)
    quant_fn <- paste0("data/processed/", e, "/quantile/", s, ".quantile_counts.rds")
    e_quant <- readRDS(quant_fn)

    test_data <- e_diem@test_data

    # Add droplet ID prefixes  
    rownames(e_diem@test_data) <- paste(dpre[s], rownames(e_diem@test_data), sep="_")
    diem_si <- intersect(rownames(e_diem@test_data), rownames(sf))
    e_diem@test_data <- e_diem@test_data[diem_si,]

    colnames(e_ed$counts) <- paste(dpre[s], colnames(e_ed$counts), sep="_")
    ed_si <- intersect(colnames(e_ed$counts), rownames(sf))
    e_ed$counts <- e_ed$counts[,ed_si]

    colnames(e_quant$counts) <- paste(dpre[s], colnames(e_quant$counts), sep="_")
    q_si <- intersect(colnames(e_quant$counts), rownames(sf))
    e_quant$counts <- e_quant$counts[,q_si]

    test_drop <- rownames(e_diem@test_data)

    diem_call <- list("Clean" = rownames(e_diem@test_data)[e_diem@test_data$Call == "Clean"],
                      "Debris" = rownames(e_diem@test_data)[e_diem@test_data$Call == "Debris"])
    ed_call <- list("Clean" = colnames(e_ed$counts),
                    "Debris" = setdiff(test_drop, colnames(e_ed$counts)))
    quant_call <- list("Clean" = colnames(e_quant$counts),
                       "Debris" = setdiff(test_drop, colnames(e_quant$counts)))

    diem_sf <- list("Clean" = data.frame("Droplet" = diem_call$Clean, 
                                         "SpliceFraction" = sf[diem_call$Clean,],
                                         "Call" = "Clean",
                                         "Method" = "DIEM"),
                    "Debris" = data.frame("Droplet" = diem_call$Debris,
                                          "SpliceFraction" = sf[diem_call$Debris,],
                                          "Call" = "Debris",
                                          "Method" = "DIEM"))
    ed_sf <- list("Clean" = data.frame("Droplet" = ed_call$Clean, 
                                       "SpliceFraction" = sf[ed_call$Clean,],
                                       "Call" = "Clean",
                                       "Method" = "EmptyDrops"),
                  "Debris" = data.frame("Droplet" = ed_call$Debris, 
                                        "SpliceFraction" = sf[ed_call$Debris,],
                                        "Call" = "Debris",
                                        "Method" = "EmptyDrops"))
    quant_sf <- list("Clean" = data.frame("Droplet" = quant_call$Clean,
                                          "SpliceFraction" = sf[quant_call$Clean,],
                                          "Call" = "Clean",
                                          "Method" = "Quantile"),
                     "Debris" = data.frame("Droplet" = quant_call$Debris, 
                                           "SpliceFraction" = sf[quant_call$Debris,],
                                           "Call" = "Debris",
                                           "Method" = "Quantile"))

    diem_df <- do.call(rbind, diem_sf) 
    diem_df[,"Droplet"] <- as.character(diem_df[,"Droplet"])
    ed_df <- do.call(rbind, ed_sf)
    ed_df[,"Droplet"] <- as.character(ed_df[,"Droplet"])
    quant_df <- do.call(rbind, quant_sf)
    quant_df[,"Droplet"] <- as.character(quant_df[,"Droplet"])

    diem_df[,"n_genes"] <- e_diem@test_data[diem_df$Droplet,"n_genes"]
    ed_df[,"n_genes"] <- e_diem@test_data[ed_df$Droplet,"n_genes"]
    quant_df[,"n_genes"] <- e_diem@test_data[quant_df$Droplet,"n_genes"]

    all_df <- as.data.frame(do.call(rbind, list(diem_df,ed_df,quant_df)))
    all_df <- all_df[!is.na(all_df$SpliceFraction),]
    all_df[,"Exprm"] <- expr[i]
    all_df[,"Sample"] <- labels[i]
    sf_calls[[s]] <- all_df
}

all_sf <- do.call(rbind, sf_calls)
all_sf[,"Sample"] <- factor(all_sf[,"Sample"], labels)
all_sf[,"Exprm"] <- factor(all_sf[,"Exprm"], c("adpcyte", "mouse_nuclei_2k", "atsn"))

all_sf <- all_sf[all_sf[,"n_genes"] >= 200,]


#=========================================
# Plot method vs fraction spliced in debris droplets
#=========================================

datf <- subset(all_sf, Call == "Debris")

p <- ggplot(datf, aes(x = Method, y = SpliceFraction)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(alpha = 0.1, width = 0.1, shape=16) +
facet_wrap(~Sample, nrow = 1) +
ylab("Fraction spliced reads\nin removed droplets") +
ylim(0,1.25) + 
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
outfn <- paste0(dir_plot, "sf.method.debris.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "sf.method.debris.pdf")
ggsave(outfn, width = 20, height = 4)

#=========================================
# Plot method vs fraction spliced in clean droplets
#=========================================

datf <- subset(all_sf, Call == "Clean")

p <- ggplot(datf, aes(x = Method, y = SpliceFraction)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(alpha = 0.1, width = 0.1, shape=16) +
facet_wrap(~Sample, nrow = 1) +
ylab("Fraction spliced reads\nin passing droplets") +
ylim(0,1.25) + 
theme_bw() +
ct +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      legend.text=element_text(size=24),
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      axis.title = element_text(size = 20),
      axis.title.y = element_text(size = 22), 
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "sf.method.clean.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "sf.method.clean.pdf")
ggsave(outfn, width = 20, height = 4)

#=========================================
# Plot number passing droplets along with
# number of nuclear droplets passing
#=========================================

datf <- subset(all_sf, Call == "Clean")

datf_g <- dplyr::group_by(datf, Method, Sample)
pct_debris <- dplyr::summarise(datf_g,
    n = sum(SpliceFraction < 0.5),
    b = sum(SpliceFraction >= 0.5))

pct_debris_m <- reshape2::melt(pct_debris, id.vars = c("Method", "Sample"))
pct_debris_m$variable = factor(pct_debris_m$variable, levels = c("b", "n"))

p <- ggplot(pct_debris_m, aes(x = Method, y = value)) +
geom_bar(aes(fill = variable), stat = "identity", position=position_stack()) +
scale_fill_discrete(labels = c("Background", "Nuclear")) +
facet_wrap(~Sample, nrow = 1, scales = "free") +
ylab("Droplets\npassing") + 
theme_bw() +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      legend.position = "top", 
      legend.title = element_blank(), 
      legend.text=element_text(size=24),
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

datf_g <- dplyr::group_by(datf, Method, Sample)
pct_debris <- dplyr::summarise(datf_g,
    n = sum(SpliceFraction < 0.5),
    b = sum(SpliceFraction >= 0.5))

pct_debris_m <- reshape2::melt(pct_debris, id.vars = c("Method", "Sample"))
pct_debris_m$variable = factor(pct_debris_m$variable, levels = c("b", "n"))

p <- ggplot(pct_debris_m, aes(x = Method, y = value)) +
geom_bar(aes(fill = variable), stat = "identity", position=position_stack()) +
scale_fill_discrete(labels = c("Background", "Nuclear")) +
facet_wrap(~Sample, nrow = 1, scales = "free") +
ylab("Droplets\nremoved") +
theme_bw() +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      legend.position = "top", 
      legend.title = element_blank(), 
      legend.text=element_text(size=24),
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      axis.title = element_text(size = 20),
      axis.title.y = element_text(size = 22),
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "method.sample.n_droplets.debris.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "method.sample.n_droplets.debris.pdf")
ggsave(outfn, width = 20, height = 4)

#=========================================
# Plot number passing droplets that are debris (type I error)
#=========================================

datf <- subset(all_sf, Call == "Clean")

datf_g <- dplyr::group_by(datf, Method, Sample, Exprm)
pct_debris <- dplyr::summarise(datf_g,
    pct_debris = sum(SpliceFraction >= 0.5)/length(SpliceFraction))

p <- ggplot(pct_debris, aes(x = Method, y = pct_debris, fill = Method)) +
geom_bar(stat = "identity") +
facet_wrap(~Sample, ncol = 1, scales = "free") +
ylab("Percent of passing droplets that are debris\n(type I error)") +
scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
theme_bw() +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      legend.position = "none", 
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "method.sample.type1.jpeg")
ggsave(outfn, width = 3.5, height = 20)
outfn <- paste0(dir_plot, "method.sample.type1.pdf")
ggsave(outfn, width = 3.5, height = 20)

# Paired t-test
d_pd <- subset(pct_debris, Exprm == "atsn" & Method == "DIEM")[,"pct_debris"]
e_pd <- subset(pct_debris, Exprm == "atsn" & Method == "EmptyDrops")[,"pct_debris"]
q_pd <- subset(pct_debris, Exprm == "atsn" & Method == "Quantile")[,"pct_debris"]

#=========================================
# Plot number removed droplets that are nuclear (type II error)
#=========================================

datf <- subset(all_sf, Call == "Debris")

datf_g <- dplyr::group_by(datf, Method, Sample)
pct_debris <- dplyr::summarise(datf_g,
    pct_debris = sum(SpliceFraction < 0.5)/length(SpliceFraction))

p <- ggplot(pct_debris, aes(x = Method, y = pct_debris, fill = Method)) +
geom_bar(stat = "identity") +
facet_wrap(~Sample, ncol = 1, scales = "free") +
ylab("Percent of removed droplets that are nuclear\n(type II error)") +
scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
theme_bw() +
theme(strip.background = element_blank(),
      strip.placement = "outside",
      text=element_text(size=ts),
      legend.position = "none", 
      legend.title = element_blank(),       
      axis.text.x = element_text(angle=xa_angle,hjust=hj),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20))

dir_plot <- "results/plots/"
outfn <- paste0(dir_plot, "method.sample.type2.jpeg")
ggsave(outfn, width = 3.5, height = 20)
outfn <- paste0(dir_plot, "method.sample.type2.pdf")
ggsave(outfn, width = 3.5, height = 20)



# Compare means paired
at_means <- sapply(at_labs, function(s){
                   datf <- all_sf[all_sf[,"Sample"] == s,]
                   qcl <- datf$Method == "Quantile" & datf$Call == "Clean"
                   ecl <- datf$Method == "EmptyDrops" & datf$Call == "Clean"
                   dcl <- datf$Method == "DIEM" & datf$Call == "Clean"
                   mq <- mean(datf[qcl,"SpliceFraction"])
                   me <- mean(datf[ecl,"SpliceFraction"])
                   md <- mean(datf[dcl,"SpliceFraction"])
                   return(c("EmptyDrops" = me,"DIEM" = md, "Quantile" = mq))
})
at_means <- t(at_means)
t.test(at_means[,"DIEM"],at_means[,"EmptyDrops"], paired=T)
t.test(at_means[,"DIEM"],at_means[,"Quantile"], paired=T)


# Compare n debris paired
at_ndebris <- sapply(at_labs, function(s){
                     datf <- all_sf[all_sf[,"Sample"] == s,]
                     qcl <- datf$Method == "Quantile" & datf$Call == "Clean"
                     ecl <- datf$Method == "EmptyDrops" & datf$Call == "Clean"
                     dcl <- datf$Method == "DIEM" & datf$Call == "Clean"
                     mq <- sum(datf[qcl,"SpliceFraction"] >= .5)
                     me <- sum(datf[ecl,"SpliceFraction"] >= .5)
                     md <- sum(datf[dcl,"SpliceFraction"] >= .5)
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
                     mq <- sum(datf[qcl,"SpliceFraction"] < .5)
                     me <- sum(datf[ecl,"SpliceFraction"] < .5)
                     md <- sum(datf[dcl,"SpliceFraction"] < .5)
                     return(c("EmptyDrops" = me, "DIEM" = md, "Quantile" = mq))
})
at_ndebris <- t(at_ndebris)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"EmptyDrops"], paired=T)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"Quantile"], paired=T)

# Compare n nuclear paired in debris
datf_s <- all_sf[all_sf$Call == "Debris",]
at_ndebris <- sapply(at_labs, function(s){
                     datf <- datf_s[datf_s[,"Sample"] == s,]
                     qcl <- datf$Method == "Quantile"
                     ecl <- datf$Method == "EmptyDrops"
                     dcl <- datf$Method == "DIEM"
                     mq <- sum(datf[qcl,"SpliceFraction"] < .5)
                     me <- sum(datf[ecl,"SpliceFraction"] < .5)
                     md <- sum(datf[dcl,"SpliceFraction"] < .5)
                     return(c("EmptyDrops" = me, "DIEM" = md, "Quantile" = mq))
})
at_ndebris <- t(at_ndebris)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"EmptyDrops"], paired=T)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"Quantile"], paired=T)

# Compare % debris droplets
at_ndebris <- sapply(at_labs, function(s){
                     datf <- all_sf[all_sf[,"Sample"] == s,]
                     qcl <- datf$Method == "Quantile" & datf$Call == "Clean"
                     ecl <- datf$Method == "EmptyDrops" & datf$Call == "Clean"
                     dcl <- datf$Method == "DIEM" & datf$Call == "Clean"
                     mq <- sum(datf[qcl,"SpliceFraction"] >= .5) / sum(qcl)
                     me <- sum(datf[ecl,"SpliceFraction"] >= .5) / sum(ecl)
                     md <- sum(datf[dcl,"SpliceFraction"] >= .5) / sum(dcl)
                     return(c("EmptyDrops" = me, "DIEM" = md, "Quantile" = mq))
})
at_ndebris <- t(at_ndebris)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"EmptyDrops"], paired=T)
wilcox.test(at_ndebris[,"DIEM"],at_ndebris[,"Quantile"], paired=T)

# Wilcox test of droplets
datf_s <- all_sf[all_sf[,"Call"] == "Clean",]
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    qcl <- datf$Method == "Quantile"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[qcl,"SpliceFraction"], datf[dcl,"SpliceFraction"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value * 8)
}
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    ecl <- datf$Method == "EmptyDrops"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[ecl,"SpliceFraction"], datf[dcl,"SpliceFraction"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value) * 8
}

datf_s <- all_sf[all_sf[,"Call"] == "Debris",]
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    qcl <- datf$Method == "Quantile"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[qcl,"SpliceFraction"], datf[dcl,"SpliceFraction"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value * 8)
}
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    ecl <- datf$Method == "EmptyDrops"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[ecl,"SpliceFraction"], datf[dcl,"SpliceFraction"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value) * 8
}



