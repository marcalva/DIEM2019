
# Compare fraction spliced in each of the experiments
# and methods

setwd("../../")

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

mt_calls <- list()

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
    rownames(e_diem@droplet_data) <- paste(dpre[s], rownames(e_diem@droplet_data), sep="_")
    colnames(e_ed$counts) <- paste(dpre[s], colnames(e_ed$counts), sep="_")
    colnames(e_quant$counts) <- paste(dpre[s], colnames(e_quant$counts), sep="_")

    test_drop <- rownames(e_diem@test_data)

    diem_call <- list("Clean" = rownames(e_diem@test_data)[e_diem@test_data$Call == "Clean"],
                      "Debris" = rownames(e_diem@test_data)[e_diem@test_data$Call == "Debris"])
    ed_call <- list("Clean" = colnames(e_ed$counts),
                    "Debris" = setdiff(test_drop, colnames(e_ed$counts)))
    quant_call <- list("Clean" = colnames(e_quant$counts),
                       "Debris" = setdiff(test_drop, colnames(e_quant$counts)))

    diem_sf <- list("Clean" = data.frame("Droplet" = diem_call$Clean, 
                                         "MT" = e_diem@test_data[diem_call$Clean,"pct.mt"], 
                                         "Call" = "Clean",
                                         "Method" = "DIEM"),
                    "Debris" = data.frame("Droplet" = diem_call$Debris,
                                          "MT" = e_diem@test_data[diem_call$Debris,"pct.mt"],
                                          "Call" = "Debris",
                                          "Method" = "DIEM"))
    ed_sf <- list("Clean" = data.frame("Droplet" = ed_call$Clean, 
                                       "MT" = e_diem@test_data[ed_call$Clean,"pct.mt"],
                                       "Call" = "Clean",
                                       "Method" = "EmptyDrops"),
                  "Debris" = data.frame("Droplet" = ed_call$Debris, 
                                        "MT" = e_diem@test_data[ed_call$Debris,"pct.mt"],
                                        "Call" = "Debris",
                                        "Method" = "EmptyDrops"))
    quant_sf <- list("Clean" = data.frame("Droplet" = quant_call$Clean,
                                          "MT" = e_diem@test_data[quant_call$Clean,"pct.mt"], 
                                          "Call" = "Clean",
                                          "Method" = "Quantile"),
                     "Debris" = data.frame("Droplet" = quant_call$Debris, 
                                           "MT" = e_diem@test_data[quant_call$Debris,"pct.mt"], 
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
    all_df <- all_df[!is.na(all_df$MT),]
    all_df[,"Exprm"] <- expr[i]
    all_df[,"Sample"] <- labels[i]
    mt_calls[[s]] <- all_df
}

all_mt <- do.call(rbind, mt_calls)
all_mt[,"Sample"] <- factor(all_mt[,"Sample"], labels)
all_mt[,"Exprm"] <- factor(all_mt[,"Exprm"], c("adpcyte", "mouse_nuclei_2k", "atsn"))

all_mt <- all_mt[all_mt[,"n_genes"] >= 200,]


#=========================================
# Plot method vs fraction spliced in debris droplets
#=========================================

datf <- subset(all_mt, Call == "Debris")

p <- ggplot(datf, aes(x = Method, y = MT)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(alpha = 0.1, width = 0.1, shape=16) +
facet_wrap(~Sample, nrow = 1, scales = "free") +
ylab("Fraction spliced reads\nin removed droplets") +
#ylim(0,125) + 
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
outfn <- paste0(dir_plot, "mt.method.debris.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "mt.method.debris.pdf")
ggsave(outfn, width = 20, height = 4)

# Wilcox test of droplets
# Negative estimate means DIEM has higher SF
for (s in labels){
    datf <- all_mt[all_mt[,"Sample"] == s,]
    qcl <- datf$Method == "Quantile" & datf$Call == "Debris"
    dcl <- datf$Method == "DIEM" & datf$Call == "Debris"
    print(s)
    tr <- wilcox.test(datf[qcl,"MT"], datf[dcl,"MT"], conf.int=T)
    message("\t", tr$estimate * 8)
    message("\t", tr$p.value * 8)
}
for (s in labels){
    datf <- all_mt[all_mt[,"Sample"] == s,]
    ecl <- datf$Method == "EmptyDrops" & datf$Call == "Debris"
    dcl <- datf$Method == "DIEM" & datf$Call == "Debris"
    print(s)
    tr <- wilcox.test(datf[ecl,"MT"], datf[dcl,"MT"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value * 8)
}

#=========================================
# Plot method vs fraction spliced in clean droplets
#=========================================

datf <- subset(all_mt, Call == "Clean")

p <- ggplot(datf, aes(x = Method, y = MT)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(alpha = 0.1, width = 0.1, shape=16) +
facet_wrap(~Sample, nrow = 1, scales = "free") +
ylab("Fraction spliced reads\nin passing droplets") +
#ylim(0,1.25) + 
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
outfn <- paste0(dir_plot, "mt.method.clean.jpeg")
ggsave(outfn, width = 20, height = 4)
outfn <- paste0(dir_plot, "mt.method.clean.pdf")
ggsave(outfn, width = 20, height = 4)


# Compare means paired
at_means <- sapply(at_labs, function(s){
                   datf <- all_mt[all_mt[,"Sample"] == s,]
                   qcl <- datf$Method == "Quantile" & datf$Call == "Clean"
                   ecl <- datf$Method == "EmptyDrops" & datf$Call == "Clean"
                   dcl <- datf$Method == "DIEM" & datf$Call == "Clean"
                   mq <- mean(datf[qcl,"MT"])
                   me <- mean(datf[ecl,"MT"])
                   md <- mean(datf[dcl,"MT"])
                   return(c("EmptyDrops" = me,"DIEM" = md, "Quantile" = mq))
})
at_means <- t(at_means)
t.test(at_means[,"DIEM"],at_means[,"EmptyDrops"], paired=T)
t.test(at_means[,"DIEM"],at_means[,"Quantile"], paired=T)

# Wilcox test of droplets
datf_s <- all_mt[all_mt[,"Call"] == "Clean",]
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    qcl <- datf$Method == "Quantile"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[qcl,"MT"], datf[dcl,"MT"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value * 8)
}
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    ecl <- datf$Method == "EmptyDrops"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[ecl,"MT"], datf[dcl,"MT"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value) * 8
}

datf_s <- all_mt[all_mt[,"Call"] == "Debris",]
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    qcl <- datf$Method == "Quantile"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[qcl,"MT"], datf[dcl,"MT"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value * 8)
}
for (s in labels){
    datf <- datf_s[datf_s[,"Sample"] == s,]
    ecl <- datf$Method == "EmptyDrops"
    dcl <- datf$Method == "DIEM"
    print(s)
    tr <- wilcox.test(datf[ecl,"MT"], datf[dcl,"MT"], conf.int=T)
    message("\t", tr$estimate)
    message("\t", tr$p.value) * 8
}



