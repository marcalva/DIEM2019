
library(diem)
library(Seurat)
library(ggplot2)

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
plot_umap_meta <- function(x, colname="percent.mt", legend_name="MT%", 
						   color_limits=NULL, color_breaks=waiver()){
	df <- data.frame(Mito=x@meta.data[,colname],
					 UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
					 UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
	p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Mito)) + 
	geom_point() + theme_classic() + 
	theme(axis.text=element_blank(), axis.ticks=element_blank()) +
	scale_colour_distiller(palette = "Spectral", name=legend_name, limits=color_limits, breaks=color_breaks) 
	return(p)
}

#==========================================================
# Normalization
#==========================================================

remove_multiplet <- function(x, lower=FALSE){
	nUMI <- x@meta.data[,"nCount_RNA"]
	upper_bound <- 10^(mean(log10(nUMI)) + 2*sd(log10(nUMI)))
	keep <- nUMI < upper_bound
	if (lower){
		lower_bound <- 10^(mean(log10(nUMI)) - 2*sd(log10(nUMI)))
		keep <- keep & (nUMI > lower_bound)
	}
	cells2keep <- rownames(x@meta.data)[keep]
	x <- subset(x, cells=cells2keep)
	return(x)
}

seurat_norm <- function(x, 
						vars.to.regress = NULL, 
						subset_mt=FALSE, 
						mt_thresh=5, 
                        method="LogNormalize", 
                        scale.factor=NULL, 
						lower=FALSE){
	x <- RenameCells(x, add.cell.id = sub("_", "-", x@project.name))
	mt_genes <- grep(pattern="^mt", rownames(x), value=TRUE, ignore.case=TRUE)
	x <- PercentageFeatureSet(x, features=mt_genes, col.name = "percent.mt")
    genes <- grep(pattern="^malat1$", rownames(x), value=TRUE, ignore.case=TRUE)
    x <- PercentageFeatureSet(x, features=genes, col.name = "MALAT1")

	# Subset by MT %
	if (subset_mt){
		keep <- (x@meta.data[,"percent.mt"] < mt_thresh)
		cells2keep <- rownames(x@meta.data)[keep]
		x <- subset(x, cells=cells2keep)
	}

    if (ncol(x) < 2) return(x)
	# x <- remove_multiplet(x, lower=lower)
	if (method == "sctransform"){
	    x <- SCTransform(x, vars.to.regress = vars.to.regress, verbose = FALSE)
    } else {
        if (is.null(scale.factor)){
            scale.factor <- median( Matrix::colSums(x@assays$RNA@counts) )
        }
	    x <- NormalizeData(x, normalization.method=method, scale.factor=scale.factor)
	    x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
	    x <- ScaleData(x)
    }
	return(x)
}

#==========================================================
# Integration
#==========================================================

seurat_merge <- function(x.list){
	merged <- merge(x.list[[1]], x.list[2:length(x.list)])
	merged <- FindVariableFeatures(merged)
	return(merged)
}

seurat_integrate <- function(x.list){
	anchors <- FindIntegrationAnchors(object.list = x.list, dims = 1:30)
	integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
	return(integrated)
}

#==========================================================
# Cluster
#==========================================================

seurat_cluster <- function(x){
	x <- ScaleData(x, verbose = TRUE)
	# vars = sapply(x@assays$RNA@var.features, function(i) var(x@assays$RNA@scale.data[i,]))
	# vf <- x@assays$RNA@var.features[vars != 0]
	# x@assays$RNA@var.features <- vf
	x <- RunPCA(x, npcs = 30, verbose = TRUE)
	
	require(reticulate)
	reticulate::use_python("/u/local/apps/python/3.6.1-shared/bin/python3.6")
	# umap <- import("umap")
	x <- FindNeighbors(x, dims = 1:30, verbose = TRUE)
	x <- FindClusters(x, verbose = TRUE)
	x <- RunUMAP(x, dims = 1:30, reduction = "pca", verbose = TRUE)
	# x <- merge_clust(x)
	return(x)
}

seurat_pipe_single <- function(counts, 
							   dir_label, 
							   project, 
							   method, 
							   meta.data=NULL, 
							   min.features=200, 
                               scale.factor=NULL){
	require(Seurat)

	# Create directories
	dp <- paste0("data/processed/", dir_label, "/", method, "/")
	dir.create(dp, recursive=TRUE, showWarnings=FALSE)
	dr <- paste0("results/", dir_label, "/", method, "/")
	dir.create(dr, recursive=TRUE, showWarnings=FALSE)
	dir_plot <- paste0(dr, "plots/")
	dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

	# Run Seurat
	seur <- CreateSeuratObject(counts, project = project, min.features=min.features, meta.data=meta.data)
	seur <- seurat_norm(seur, scale.factor=scale.factor)
	seur <- seurat_cluster(seur)

	# Save
	saveRDS(seur, paste0(dp, project, ".seur_obj.rds"))

	# Plot
	pdfname <- paste0(dir_plot, project, ".seur_clusters.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_clusters.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_labels(seur)
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	pdfname <- paste0(dir_plot, project, ".seur_mt_pct.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_mt_pct.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_meta(seur)
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	pdfname <- paste0(dir_plot, project, ".seur_malat1.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_malat1.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_meta(seur, colname="MALAT1", legend_name="MALAT1")
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	pdfname <- paste0(dir_plot, project, ".seur_ngene.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_ngene.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_meta(seur, colname="nFeature_RNA", legend_name="n_genes")
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	#pdfname <- paste0(dir_plot, project, ".seur_DebrisLogOdds.pdf")
	#jpgname <- paste0(dir_plot, project, ".seur_DebrisLogOdds.jpeg")
	#pdf(pdfname, width=9,height=9)
	#p <- plot_umap_meta(seur, colname="DebrisLogOdds", legend_name="Debris log odds")
	#print(p)
	#dev.off()
	#system(paste("convert", "-density", "200", pdfname, jpgname))

	markers <- FindAllMarkers(seur, only.pos=TRUE)
	write.table(markers, paste0(dr, project, ".seur_markers.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

seurat_pipe_list <- function(counts.l, 
							 dir_label, 
							 project.l, 
							 project, 
							 method, 
							 meta.data.l=NULL, 
							 min.features=200, 
                             scale.factor=NULL){
	require(Seurat)

	# Create directories
	dp <- paste0("data/processed/", dir_label, "/", method, "/")
	dir.create(dp, recursive=TRUE, showWarnings=FALSE)
	dr <- paste0("results/", dir_label, "/", method, "/")
	dir.create(dr, recursive=TRUE, showWarnings=FALSE)
	dir_plot <- paste0(dr, "plots/")
	dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

	# Run Seurat

	seur.list <- lapply(1:length(counts.l), function(i) {
		CreateSeuratObject(counts.l[[i]], project = project.l[[i]], min.features=min.features, meta.data=meta.data.l[[i]])
		})
	seur.list <- lapply(seur.list, seurat_norm, scale.factor=scale.factor)
	seur <- seurat_merge(seur.list)
	seur <- seurat_cluster(seur)

	# Save
	saveRDS(seur, paste0(dp, project, ".seur_obj.rds"))

	# Plot
	pdfname <- paste0(dir_plot, project, ".seur_clusters.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_clusters.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_labels(seur)
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	pdfname <- paste0(dir_plot, project, ".seur_mt_pct.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_mt_pct.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_meta(seur)
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	pdfname <- paste0(dir_plot, project, ".seur_malat1.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_malat1.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_meta(seur, colname="MALAT1", legend_name="MALAT1")
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	pdfname <- paste0(dir_plot, project, ".seur_ngene.pdf")
	jpgname <- paste0(dir_plot, project, ".seur_ngene.jpeg")
	pdf(pdfname, width=9,height=9)
	p <- plot_umap_meta(seur, colname="nFeature_RNA", legend_name="n_genes")
	print(p)
	dev.off()
	system(paste("convert", "-density", "200", pdfname, jpgname))

	#pdfname <- paste0(dir_plot, project, ".seur_DebrisLogOdds.pdf")
	#jpgname <- paste0(dir_plot, project, ".seur_DebrisLogOdds.jpeg")
	#pdf(pdfname, width=9,height=9)
	#p <- plot_umap_meta(seur, colname="DebrisLogOdds", legend_name="Debris log odds")
	#print(p)
	#dev.off()
	#system(paste("convert", "-density", "200", pdfname, jpgname))

	#markers <- FindAllMarkers(seur, only.pos=TRUE)
	#write.table(markers, paste0(dr, project, ".seur_markers.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}
