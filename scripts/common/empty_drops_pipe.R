
empty_drops_pipe <- function(counts, 
							 dir_label, 
							 project, 
							 method="emptydrops", 
							 fdr_thresh=0.05, 
							 seedn=100, 
                             ...){
	suppressMessages(require(DropletUtils))

	dp <- paste0("data/processed/", dir_label, "/", method, "/")
	dir.create(dp, recursive=TRUE, showWarnings=FALSE)
	dr <- paste0("results/", dir_label, "/", method, "/")
	dir.create(dr, recursive=TRUE, showWarnings=FALSE)
	dir_plot <- paste0(dr, "plots/")
	dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

	set.seed(seedn)
	e.out <- emptyDrops(counts, ...)
	write.table(e.out, paste0(dr, project, ".", method, ".out.txt"), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

	is_cell <- e.out$FDR < fdr_thresh
	is_cell[is.na(is_cell)] <- FALSE
	counts <- counts[,is_cell]

	pdf(paste0(dir_plot, project, ".", method, ".umi.lprob.pdf"))
	plot(log(e.out$Total), -e.out$LogProb, col=ifelse(is_cell, "red", "black"),
    	xlab="Total UMI count", ylab="-Log Probability")
	dev.off()

	ret <- list(counts=counts, e.out=e.out)
	saveRDS(ret, paste0(dp, project, ".emptydrops_counts.rds"))

	return(ret)
}

