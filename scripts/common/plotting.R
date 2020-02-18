
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

# x is a seurat object
boxplot_gene <- function(x, gene, ct_col="RNA_snn_res.0.8", assay = NULL){
	require(reshape2)
    expr <- GetAssayData(x, slot="data", assay=assay)
    feat <- expr[gene,]
    ct <- x@meta.data[,ct_col]
    datf <- data.frame("Gene" = feat, "CellType" = ct)

	df <- reshape2::melt(datf)
	p <- ggplot(datf, aes(x=CellType, y=Gene)) +
	geom_boxplot() + theme_classic() + 
    ggtitle(gene) + 
    ylab("Expression") + 
    theme(legend.position = "none")
	return(p)
}

# Plot gene expression values in a UMAP
# x is a seurat object
# returns a ggplot object
plot_umap_gene <- function(x, 
                           gene, 
                           legend_title="UMI", 
                           palette="GnBu", 
                           assay=NULL, 
                           alpha=1, 
                           shape = 16, 
                           size = 8, 
                           rand=TRUE){
    require(ggplot2)
    require(RColorBrewer)

    expr <- GetAssayData(x, slot="data", assay=assay)

    datf <- data.frame(feat=expr[gene,],
                       UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
                       UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
    if (rand) datf <- datf[sample(rownames(datf)),]
    p <- ggplot(datf, aes(x=UMAP1, y=UMAP2, color=feat)) + 
    geom_point(alpha=alpha, shape = shape, size = size) + 
    theme_classic() + ggtitle(gene) + 
    theme(axis.text=element_blank(), axis.ticks=element_blank()) +
    scale_colour_distiller(palette = palette, name=legend_title)
    return(p)
}

