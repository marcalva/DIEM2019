
setwd("../")

library(diem)
suppressMessages(library(DropletUtils))
library(ggplot2)
library(ggpubr)
library(gridExtra)

labl <- "adpcyte"
dp <- paste0("data/processed/", labl, "/")

#=========================================
#=========================================

# Read in DIEM
sce <- readRDS("data/processed/adpcyte/diem/adpcyte.diem_sce.rds")

p1 <- barcode_rank_plot(sce, title="Barcode-Rank", ret=TRUE) + 
geom_vline(xintercept=1e4, linetype="dashed", color = "red") + 
theme(plot.title=element_text(face="bold")) + 
theme(text=element_blank())

pdfname <- "results/plots/fig3a_top.pdf"
pdf(pdfname, width=10, height=4)
p1
dev.off()

sigf <- function(x, m=1e4, a=5e-4) {
	ret <- 1/(1+exp(-a*(x-m)))
	return(ret)
	if (is.na(ret)) return(1)
	else return(ret)
}

base=1.025
umis <- seq(1, log(ncol(sce), base=base), 1)
pps <- sapply(base^umis, sigf, 1e3, 0.8e-3)
pps <- pps - min(pps)
pps <- pps * (1/max(pps))
# pps[base^umis > 1e4] <- 1
# pps[base^umis < 1e2] <- 0 + pps[base^umis < 1e2]/10


xd <- matrix(nrow=2, ncol=length(umis))
rownames(xd) <- c("Nucleus", "Background")

# colnames(xd) <- umis
xd[1,] <- 1-pps
xd[2,] <- pps

hmcol = colorRampPalette(c("white", "grey"))(100)

xdm <- reshape2::melt(xd)
p2 <- ggplot(xdm, aes(x=Var2, y=Var1, fill=value)) + 
geom_tile(width=2) + 
theme_classic() + 
scale_fill_gradient(low="white", high="black", name="Posterior\nProbability") + 
xlab("") + ylab("") + 
theme(axis.text.x=element_blank(), 
	  axis.ticks=element_blank(), 
	  axis.line=element_blank(), 
	  legend.position="none") + 
theme(text=element_blank())

pdfname <- "results/plots/fig3a_bottom.pdf"
pdf(pdfname, width=10, height=6)
ggarrange(p1, p2, nrow=2, heights=c(4,1), align="v")
dev.off()


# Create probability plot

n_cells <- 50
pd <- matrix(nrow=2, ncol=n_cells)
rownames(pd) <- c("Nucleus", "Background")
set.seed(1)
pd[1,] <- runif(n_cells)
pd[2,] <- runif(n_cells)

pdm <- reshape2::melt(pd)
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(100))
p3 <- ggplot(pdm, aes(x=Var2, y=Var1, fill=value)) + 
geom_tile(width=2) + 
theme_classic() + 
scale_fill_gradientn(colours=hmcol) + 
xlab("") + ylab("") + 
theme(axis.text.x=element_blank(), 
	  axis.ticks=element_blank(), 
	  axis.line=element_blank(), 
	  legend.position="none") + 
theme(text=element_blank())


pdfname <- "results/plots/fig3a_middle.pdf"
pdf(pdfname, width=7, height=2)
p3
dev.off()


