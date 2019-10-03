
quantile_pipe <- function(counts, 
						  dir_label, 
						  project, 
						  top_n=3000, 
						  method="quantile"){

	dp <- paste0("data/processed/", dir_label, "/", method, "/")
	dir.create(dp, recursive=TRUE, showWarnings=FALSE)
	dr <- paste0("results/", dir_label, "/", method, "/")
	dir.create(dr, recursive=TRUE, showWarnings=FALSE)
	dir_plot <- paste0(dr, "plots/")
	dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

	counts <- counts[,Matrix::colSums(counts) > 0]
	tc <- Matrix::colSums(counts)
	tc_top <- tc[order(tc, decreasing=TRUE)[1:top_n]]
	m <- quantile(tc_top, probs=0.99)
	thresh <- m/10

	tco <- tc[order(tc, decreasing=TRUE)]
	tco[duplicated(tco)] <- NA
	dat <- data.frame(Rank=1:length(tco), Total=tco)
	# Plot
	pdf(paste0(dir_plot, project, ".quantile_thresh.pdf"))
	plot(dat$Rank, dat$Total, log="xy", xlab="Rank", ylab="Total", main="Quantile Threshold")
	abline(h=thresh, col="dodgerblue", lty=2)
	dev.off()

	keep <- tc >= thresh
	counts <- counts[,keep]
	ret <- list(counts=counts, thresh=thresh)
	saveRDS(ret, paste0(dp, project, ".quantile_counts.rds"))

	return(ret)
}
