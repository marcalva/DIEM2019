
# x is a seurat object
plot_umap_meta <- function(x, 
                           colname="percent.mt", 
                           legend_name="MT%", 
                           palette = "Spectral", 
                           size=1, 
                           alpha=1, 
						   color_limits=NULL, 
                           color_breaks=waiver(), 
                           order_col = FALSE){
	df <- data.frame(Mito=x@meta.data[,colname],
					 UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
					 UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
    if (order_col) df <- df[order(df$Mito, decreasing=FALSE),,drop=FALSE]
	p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Mito)) + 
	geom_point(size=size, alpha=alpha) + theme_classic() + 
	theme(axis.text=element_blank(), axis.ticks=element_blank()) +
	scale_colour_distiller(palette = palette, direction=1, 
                           name=legend_name, limits=color_limits, 
                           breaks=color_breaks) 
	return(p)
}

plot_umap_labels <- function(x, ct_col="RNA_snn_res.0.8", legend_title="Cell Type"){
	df <- data.frame(Cluster=x@meta.data[,ct_col], 
					 UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
					 UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
	# Get median UMAP points for each cluster
	clusters <- sort(unique(df$Cluster))
	meds <- lapply(clusters, function(ct) {
				   umap <- df[df$Cluster == ct, c("UMAP1", "UMAP2")]
				   meds <- apply(umap, 2, median)
				   return(meds)
					 })
	names(meds) <- clusters
	p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Cluster)) + 
	geom_point() + theme_classic() + scale_color_discrete(name = legend_title) + 
	theme(axis.text=element_blank(), axis.ticks=element_blank())
	for (ct in clusters){
		p <- p + annotate(geom="label", x=meds[[ct]][1], y=meds[[ct]][2], label=as.character(ct))
	}
	return(p)
}

# x is a seurat object
boxplot_meta <- function(x, colname="percent.mt", ct_col="RNA_snn_res.0.8", yname="MT%"){
	require(reshape2)
	df <- data.frame(Mito=x@meta.data[,colname], Cluster=x@meta.data[,ct_col])
	df <- reshape2::melt(df)
	p <- ggplot(df, aes(x=Cluster, y=value, fill=Cluster)) +
	geom_boxplot() + theme_classic() + 
	ylab(yname)
	return(p)
}

