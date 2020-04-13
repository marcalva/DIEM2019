
setwd("../")

library(diem)
library(ggplot2)
library(scales)
library(ggrepel)
library(dplyr)

#=========================================
#=========================================

infn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(infn, header = TRUE, row.names = 1)

infn <- "data/raw/splice_frctn/midpoints.txt"
types <- read.table(infn, header = FALSE, row.names = 1)
types[,1] <- 100 * types[,1]

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
expr <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- expr
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples
mthds <- c("DIEM", "EmptyDrops", "Quantile")

#=========================================
#=========================================

midpoints <- read.table("data/raw/splice_frctn/midpoints.txt",
                        stringsAsFactors=FALSE)
rownames(midpoints) <- midpoints[,1]
midpoints[,1] <- labels[midpoints[,1]]
colnames(midpoints) <- c("Label", "h")
midpoints[,"Label"] <- factor(midpoints[,"Label"], levels = labels)


datfl <- list()

for (i in 1:length(expr)){
    s <- samples[i]
    e <- names(samples)[i]

    # Read in results
    diem_fn <- paste0("data/processed/", e, "/diem/", s, ".diem_sce.rds")
    sce <- readRDS(diem_fn)
    test_data <- sce@test_data

    clusters <- sort(unique(test_data$Cluster))
    clustersc <- as.character(clusters)
    avg_sf <- tapply(test_data$SpliceFrctn, test_data$Cluster, mean, na.rm=TRUE)
    avg_ng <- tapply(test_data$n_genes, test_data$Cluster, mean, na.rm=TRUE)

    tg <- dplyr::group_by(test_data, Cluster)
    avgs <- tg %>% summarise(Label = labels[i], 
                             ngenes = mean(n_genes, na.rm = TRUE), 
                             SpliceFrctn = mean(SpliceFrctn, na.rm = TRUE), 
                             Midpoint = midpoints[s,"h"])
    datfl[[s]] <- avgs                             

}

datf <- do.call(rbind, datfl)
datf <- datf %>% mutate(Label = factor(Label, levels = labels))
datf <- datf %>% mutate(SpliceFrctn = SpliceFrctn/100)



p <- ggplot(datf, aes(x = ngenes, y = SpliceFrctn, label = Cluster)) + 
geom_abline(aes(intercept = Midpoint, slope = 0), col = "red") + 
geom_point(shape = 16) + 
geom_text_repel() + 
scale_y_continuous(labels = scales::percent_format(accuracy=1),
                   breaks = c(0.2, 0.4, 0.6, 0.8, 1),
                   limits = c(.2, 1)) + 
scale_x_log10() + 
facet_wrap(~Label, ncol = 4, scales = "free") + 
theme_minimal() + 
xlab("Average number of genes detected") + 
ylab("Average percent reads spliced") + 
theme(text = element_text(size = 16), 
      axis.text.x = element_text(angle = 30, hjust = 0.9))


fn <- "results/plots/ngenes.pct_splice.init_clust.pdf"
ggsave(fn, width = 10, height = 6)
fn <- "results/plots/ngenes.pct_splice.init_clust.jpeg"
ggsave(fn, width = 10, height = 6)


