
setwd("../")

library(diem)
library(ggplot2)
library(ggpubr)


set_breaks_10 <- function(x){
    xmax <- x[2]
    bk <- 10
    brks <- c(bk)
    while (bk < xmax){
        bk <- bk * 10
        brks <- c(brks, bk)
    }
    return(brks)
}

#====================================================
# Adipose
#====================================================

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

all_sf[,"Exprm"] <- factor(all_sf[,"Exprm"], levels=unique(expr))
all_sf[,"Sample"] <- factor(all_sf[,"Sample"], levels=samples)
all_sf[,"Label"] <- factor(all_sf[,"Label"], levels=labels)

all_sf[,"SpliceFrctn"] <- all_sf[,"SpliceFrctn"] / 100
all_sf[,"pct.mt"] <- all_sf[,"pct.mt"] / 100
all_sf[,"MALAT1"] <- all_sf[,"MALAT1"] / 100


midpoints <- read.table("data/raw/splice_frctn/midpoints.txt", 
                        stringsAsFactors=FALSE)
rownames(midpoints) <- midpoints[,1]
midpoints[,1] <- labels[midpoints[,1]]
colnames(midpoints) <- c("Label", "h")
midpoints[,"Label"] <- factor(midpoints[,"Label"], levels = labels)

#=========================================
#=========================================

k <- all_sf[,"Exprm"] == "atsn" & all_sf[,"pct.mt"] > 5
all_sf[k,"pct.mt"] <- NA

k <- all_sf[,"Exprm"] == "mouse_nuclei_2k" & all_sf[,"MALAT1"] >20
all_sf[k,"MALAT1"] <- NA

lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

library(ggplot2)

#=========================================
# Total counts
#=========================================

p_tc <- ggplot(all_sf, aes(x = SpliceFrctn, y = total_counts)) +
geom_point(shape=16, alpha = 0.01) +
facet_wrap(~Label, ncol=1) + 
geom_vline(data = midpoints, aes(xintercept = h), col="red") + 
scale_x_continuous(labels = scales::percent_format(accuracy = 1), 
                   breaks = c(0, .25, .5, .75, 1)) +
scale_y_continuous(trans='log10', 
                   breaks=set_breaks_10, 
                   labels=scales::comma) +
xlab("Percent reads spliced") +
ylab("Total counts") +
theme_minimal() + theme(text=element_text(size=22))

#=========================================
# MALAT1
#=========================================

p_malat <- ggplot(all_sf, aes(x = SpliceFrctn, y = MALAT1)) + 
geom_point(shape=16, alpha = 0.01) +  
facet_wrap(~Label, ncol=1, scales = "free_y") +
geom_vline(data = midpoints, aes(xintercept = h), col="red") + 
scale_x_continuous(labels = scales::percent_format(accuracy = 1), 
                   breaks = c(0, .25, .5, .75, 1)) +
scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
xlab("Percent reads spliced") + 
ylab("MALAT1%") + 
theme_minimal() + theme(text=element_text(size=22))

#=========================================
# MT%
#=========================================

p_mt <- ggplot(all_sf, aes(x = SpliceFrctn, y = pct.mt)) + 
geom_point(shape=16, alpha = 0.01) +  
facet_wrap(~Label, ncol=1, scales = "free_y") + 
geom_vline(data = midpoints, aes(xintercept = h), col="red") + 
scale_x_continuous(labels = scales::percent_format(accuracy = 1), 
                   breaks = c(0, .25, .5, .75, 1)) +
scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
xlab("Percent reads spliced") + 
ylab("MT%") + 
theme_minimal() + theme(text=element_text(size=22))

#=========================================
# %SF
#=========================================

p_dens <- ggplot(all_sf, aes(x = SpliceFrctn)) + 
facet_wrap(~Label, ncol=1, scales = "free_y") +
geom_density() + 
geom_vline(data = midpoints, aes(xintercept = h), col="red") + 
scale_x_continuous(labels = scales::percent_format(accuracy = 1), 
                   breaks = c(0, .25, .5, .75, 1)) +
ylab("Density") + 
xlab("Percent reads spliced") + 
theme_minimal() + theme(text=element_text(size=22))

library(ggpubr)

# Arrange into plot
fl <- list(size = 20, face="bold")
dir_plot <- paste0("results/plots/"); dir.create(dir_plot, showWarnings=FALSE, recursive=TRUE)
pdfname <- paste0(dir_plot, "sf_relat.pdf")
jpgname <- paste0(dir_plot, "sf_relat.jpeg")
pdf(pdfname, width=18, height=14)
ggarrange(p_tc, p_dens, p_mt, p_malat, labels=c("a", "b", "c", "d"), font.label = fl,ncol=4)
dev.off()
system(paste("convert", "-density", "200", pdfname, jpgname))


