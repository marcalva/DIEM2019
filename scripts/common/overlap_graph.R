
# Plot an ovelap graph for clusters
# x is a frequency table data frame
# cols1 contains the intensity of colors to plot for the row values of x
# cols2 contains the intensity of colors to plot for the column values of x
overlap_graph <- function(x, cols1, cols2, labels1=NULL, labels2=NULL, col_title=NULL, main="", 
                          scale_size=1){
	# Order data frame
	r_max <- apply(x, 1, max)
	x <- x[order(r_max, decreasing=TRUE),]
	cols1 <- cols1[order(r_max, decreasing=TRUE)]
	
	c_ord <- apply(x, 2, function(x){ return(rev(1:length(x)) %*% x)})
	x <- x[,order(c_ord, decreasing=TRUE)]
	cols2 <- cols2[order(c_ord, decreasing=TRUE)]

	# Scale colors
	minc <- min(cols1, cols2)
	cols1 <- cols1 - minc
	cols2 <- cols2 - minc
	maxc <- max(cols1, cols2)
	cols1 <- cols1/maxc
	cols2 <- cols2/maxc

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

	# Annotate labels for inflection
	for (i in 1:nrow(x)){
		xy <- c(1, ny-i)
		dpoints <- data.frame(x1 = 1, y1 = ny-i, col=cols1[i])
		p <- p + geom_point(aes(x=x1, y=y1, fill=col), shape = 21, size=sz*1.5, stroke=1, color="black", data=dpoints)
		p <- p + annotate("text", x=xy[1], y=xy[2], label=as.character(rownames(x)[i]), 
			size=sz)
	}
	# 

	# Annotate labels for DIEM
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
