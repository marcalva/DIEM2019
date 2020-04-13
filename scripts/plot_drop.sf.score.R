
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

all_sf <- read.table("data/processed/meta_data.txt",
                     header=TRUE,
                     sep="\t")

# all_sf <- all_sf[all_sf[,"n_genes"] >= 200,]

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


p <- ggplot(all_sf, aes(x = score.debris, y = SpliceFrctn)) + 
geom_point(shape = 16, alpha = 0.05) + 
scale_y_continuous(labels = scales::percent_format(accuracy=1),
                   breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                   limits = c(0, 1)) + 
facet_wrap(~Label, ncol = 2, scales = "free") + 
geom_abline(data = midpoints, aes(intercept = h, slope = 0), col = "red") + 
geom_vline(aes(xintercept = 0.5), col = "blue") + 
theme_minimal() + 
ylab("Percent reads spliced") + 
xlab("Debris score") + 
theme(text = element_text(size = 16), 
      axis.text.x = element_text(angle = 30, hjust = 0.9))


fn <- "results/plots/drops.pct_splice.debris_score.pdf"
ggsave(fn, width = 6, height = 10)
fn <- "results/plots/drops.pct_splice.debris_score.jpeg"
ggsave(fn, width = 6, height = 10)

cors <- sapply(labels, function(l){
               datf <- subset(all_sf, Label == l)
               cor(datf[,"SpliceFrctn"], datf[,"score.debris"])
      })

summary(cors)


