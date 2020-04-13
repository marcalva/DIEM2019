
setwd("../")

library(diem)
library(Seurat)
library(reshape2)
library(ggplot2)
library(gplots)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(viridis)
source("scripts/common/plotting.R")
source("scripts/common/overlap_graph.R")


#=========================================
# Functions
#=========================================

ct <- theme(text = element_text(size = 10), 
            axis.title = element_text(size = 14), 
            plot.title = element_text(size = 16, hjust = 0.5))

plot_overlap <- function(datf, col1 = "DIEM Cluster", col2 = "Seurat Cluster", title = ""){
    require(viridis)
    datfm <- reshape2::melt(as.matrix(datf * 100))
    datfm[,1] <- factor(datfm[,1])
    datfm[,2] <- factor(datfm[,2])
    p <- ggplot(datfm, aes_string(x = "Var2", y = "Var1")) + 
    geom_tile(aes_string(fill="value"), colour="white") + 
    scale_fill_viridis(name = "Percent\nOverlap", limits = c(0,100)) + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) + 
    ggtitle(title) + 
    scale_x_discrete(position = "top") + 
    ylab(col1) + 
    xlab(col2) + 
    ct
    return(p)
}

#=========================================
#=========================================

dir_plot <- "results/plots/rgrs.lk.clust/"
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)
dir_plot2 <- "results/plots/rgrs.sf/"
dir.create(dir_plot2, showWarnings = FALSE, recursive = TRUE)
dir_plot3 <- "results/plots/score.sf/"
dir.create(dir_plot3, showWarnings = FALSE, recursive = TRUE)


w <- 4
h <- 3.5

#=========================================
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

pct_ovrlp <- c()
pct_cor <- c()
score_cor <- c()

for (i in 1:length(expr)){
    s <- samples[i]
    e <- names(samples)[i]
    l <- labels[i]

    # Read in results
    diem_fn <- paste0("data/processed/", e, "/diem/", s, ".diem_sce.rds")
    sce <- readRDS(diem_fn)

    sce@test_data[,"SpliceFrctn"] <- sce@droplet_data[rownames(sce@test_data),"SpliceFrctn"]
    testdat <- sce@test_data
    #testdat <- testdat[testdat$Call == "Clean",]
    testdat <- testdat[testdat$n_genes >= 200,]

    score_means <- tapply(testdat$score.debris, testdat$Cluster, mean)
    deb_clust <- names(score_means)[score_means >= score_means[1]]
    sce <- estimate_db_pct(sce, deb_clust=deb_clust)
    testdat <- sce@test_data
    testdat <- testdat[testdat$n_genes >= 200,]
    props <- sce@prop[,rownames(testdat)]


    regres.clust <- apply(props, 2, which.max)

    cmp <- table(testdat$Cluster, regres.clust)
    cmp <- cmp %*% diag(1/colSums(cmp))

    pct_ovrlp[i] <- round(mean(diag(cmp)), 5)
    message(format(Sys.time(), "%a %b %d %H:%M:%S %Y"), 
            "\tSample ", s, " avg pct overlap: ", 
            pct_ovrlp[i]) 

    dimnames(cmp) <- NULL
    # cmpm <- melt(cmp)

    p <- plot_overlap(cmp, 
                      "Likelihood-based\ncluster estimate", 
                      "Regression-based\ncluster estimate", 
                      title = l)
    outfn <- paste0(dir_plot, "overlap.", s, ".jpeg")
    ggsave(outfn, width = w, height = h)
    outfn <- paste0(dir_plot, "overlap.", s, ".pdf")
    ggsave(outfn, width = w, height = h)

    # Plot correlation of percent debris with %SF
    pct_cor[i] <- cor(testdat[,"pct.debris"], 
                     sce@droplet_data[rownames(testdat),"SpliceFrctn"], 
                     use = "p")
    message(format(Sys.time(), "%a %b %d %H:%M:%S %Y"), 
            "\tSample ", s, " SF R2: ", 
            pct_cor[i]) 


    datf2 <- data.frame("PercentDebris" = testdat[,"pct.debris"],
                        "SpliceFrctn" = sce@droplet_data[rownames(testdat),"SpliceFrctn"])
    ggplot(datf2, aes(x = SpliceFrctn, y = PercentDebris)) + 
    geom_point(shape = 16, size = 1, alpha = 0.5) + 
    xlab("Percent reads spliced") + 
    ylab("Percent debris") + 
    xlim(0,100) + ylim(0,100) + 
    theme_bw()
    outfn <- paste0(dir_plot2, "PercentDebris.SpliceFrctn.", s, ".jpeg")
    ggsave(outfn, width = w, height = h)
    outfn <- paste0(dir_plot2, "PercentDebris.SpliceFrctn.", s, ".pdf")
    ggsave(outfn, width = w, height = h)

    # Plot correlation of debris score with %SF
    score_cor[i] <- cor(testdat[,"score.debris"], 
                     sce@droplet_data[rownames(testdat),"SpliceFrctn"], 
                     use = "p")
    message(format(Sys.time(), "%a %b %d %H:%M:%S %Y"), 
            "\tSample ", s, " Score R2: ", 
            score_cor[i]) 


    datf2 <- data.frame("PercentDebris" = testdat[,"score.debris"],
                        "SpliceFrctn" = sce@droplet_data[rownames(testdat),"SpliceFrctn"])
    ggplot(datf2, aes(x = SpliceFrctn, y = PercentDebris)) + 
    geom_point(shape = 16, size = 1, alpha = 0.5) + 
    xlab("Percent reads spliced") + 
    ylab("Debris score") + 
    xlim(0,100) + 
    theme_bw()
    outfn <- paste0(dir_plot3, "DebrisScore.SpliceFrctn.", s, ".jpeg")
    ggsave(outfn, width = w, height = h)
    outfn <- paste0(dir_plot3, "DebrisScore.SpliceFrctn.", s, ".pdf")
    ggsave(outfn, width = w, height = h)
}

