
# Compare clusters and return data frame for plotting
#' Compare 2 cluster assignments
#'
#' Determine the overlap of observations in \code{clust1} assignments 
#' with \code{clust2} assignments. \code{clust1} and \code{clust2} are 
#' vectors that should be of the same length. The index elements of 
#' \code{clust1} and \code{clust2} must correspond with each other, e.g. 
#' the first element is the same observation in both. 
#' This returns a data frame with 
#' the percent of observations in \code{clust1} that are found in 
#' \code{clust2}. The rows are \code{clust1} and the columns are 
#' \code{clust2}. The rows sum to 1.
#' 
#' @param clust1 Vector of first cluster assignments.
#' @param clust2 Vector of second cluster assignments.
#' 
#' @return Data frame of percent overlap of clust2 in clust1.
comp_clust <- function(clust1, clust2){

    if (length(clust1) != length(clust2)){
        stop("Length of clust1 and clust2 should be the same.")
    }

	ct1 <- as.character(sort(as.numeric(unique(clust1))))
	ct2 <- as.character(sort(as.numeric(unique(clust2))))

	dfc <- as.data.frame(matrix(nrow=length(ct1), ncol=length(ct2)))
	rownames(dfc) <- ct1
	colnames(dfc) <- ct2
	for (c1 in ct1){
        a <- which(clust1 == c1)
        n_a <- sum(a)
		for (c2 in ct2){
            b <- which(clust2 == c2)
            ab <- a & b
			n_ab <- sum(ab)
			dfc[c1, c2] <- n_ab/n_a
		}
	}
	return(dfc)
}

#' Plot an overlap graph for clusters
#' 
#' Plot the overlap between cluster assignments. Each unique cluster 
#' is plotted as a circle, with the first cluster assignments on the 
#' left and the second on the right. Lines connect each circle, and the 
#' thickness of the cirle portrays the percent of observations in the 
#' the \code{clust1} assignment correspond to the connecting assingment 
#' in \code{clust2}.
# x is a frequency table data frame
# cols1 contains the intensity of colors to plot for the row values of x
# cols2 contains the intensity of colors to plot for the column values of x
overlap_graph <- function(x, 
                          cols1, 
                          cols2, 
                          labels1=NULL, 
                          labels2=NULL, 
                          col_title=NULL, 
                          main="", 
                          scale_size=1){
	# Order data frame
	r_max <- apply(x, 1, max)
	x <- x[order(r_max, decreasing=TRUE),]
	cols1 <- cols1[order(r_max, decreasing=TRUE)]
	
	c_ord <- apply(x, 2, function(x){ return(rev(1:length(x)) %*% x)})
	x <- x[,order(c_ord, decreasing=TRUE)]
	cols2 <- cols2[order(c_ord, decreasing=TRUE)]

	# Scale colors
	#minc <- min(cols1, cols2)
	#cols1 <- cols1 - minc
	#cols2 <- cols2 - minc
	#maxc <- max(cols1, cols2)
	#cols1 <- cols1/maxc
	#cols2 <- cols2/maxc

	ny <- max(nrow(x), ncol(x))
	xlim <- c(0.75,2.25)

	# sz <- 12*(6/ny)
	sz = 12

	p <- ggplot() + geom_point() + theme_void()

	# Plot lines
	for (i in 1:nrow(x)){
		for (j in 1:ncol(x)){
			# if (x[i,j] == 0) next
			dpoints <- data.frame(x1 = 1, x2 = 2,y1 = ny-i, y2 = ny-j, alpha=x[i,j], size=x[i,j])
			p <- p + geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2, alpha=alpha, size=size), 
				data=dpoints, show.legend=FALSE)
		}
	}
	p <- p + scale_alpha_continuous(range = c(0, 1))
	p <- p + scale_size_continuous(range = c(0, 2))

	# Annotate labels for c1
	for (i in 1:nrow(x)){
		xy <- c(1, ny-i)
		dpoints <- data.frame(x1 = 1, y1 = ny-i, col=cols1[i])
		p <- p + geom_point(aes(x=x1, y=y1, fill=col), shape = 21, size=sz*1.5, stroke=1, color="black", data=dpoints)
		p <- p + annotate("text", x=xy[1], y=xy[2], label=as.character(rownames(x)[i]), 
			size=sz)
	}
	# 

	# Annotate labels for c2
	for (i in 1:ncol(x)){
		xy <- c(2, ny-i)
		dpoints <- data.frame(x2 = 2, y2 = ny-i, col=cols2[i])
		p <- p + geom_point(aes(x=x2, y=y2, fill=col), shape = 21, size=sz*1.5, stroke=1, color="black", data=dpoints)
		p <- p + annotate("text", x=xy[1], y=xy[2], label=as.character(colnames(x)[i]), 
			size=sz)
	}

	p <- p + scale_x_continuous(position = "top", limits=xlim, 
		breaks=c(1,2),
		labels=c(labels1, labels2)) + 
	scale_fill_continuous(low="white", high="firebrick2", name=col_title) + 
	ylim(-0.5,ny-0.5) + 
	theme(legend.position="bottom", 
		legend.text=element_text(angle=270, hjust=0.5, vjust=0.5), 
		text=element_text(size=scale_size*16),
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.border = element_rect(fill=NA, color="black"), 
		axis.ticks=element_blank(), 
		axis.text.y=element_blank(), 
		axis.text.x=element_text(size=scale_size*20), 
		axis.title=element_blank(), 
		plot.title=element_text(size=scale_size*24, hjust=0.5, face="bold"),
		plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
	p <- p + ggtitle(main)

	return(p)
}
