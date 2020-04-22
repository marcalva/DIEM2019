
diem_pipe1 <- function(counts, 
                       dir_label, 
                       project, 
                       method="diem", 
                       k_init = 20, 
                       threads = 1,
                       model = "mltn", 
                       ...){
    require(diem)

    dp <- paste0("data/processed/", dir_label, "/", method, "/")
    dir.create(dp, recursive=TRUE, showWarnings=FALSE)
    dr <- paste0("results/", dir_label, "/", method, "/")
    dir.create(dr, recursive=TRUE, showWarnings=FALSE)
    dir_plot <- paste0(dr, "plots/")
    dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

    sce <- create_SCE(counts, name=project)
    
    mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
    malat <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=malat, name="MALAT1")

    sce <- set_debris_test_set(sce)
    sce <- filter_genes(sce)
    sce <- get_pcs(sce)
    sce <- init(sce, k_init = k_init, model = model)
    sce <- run_em(sce, threads = threads, model = model)

    saveRDS(sce, paste0(dp, project, ".diem_sce.rds"))
    return(sce)
}

diem_pipe2 <- function(sce, 
                       dir_label, 
                       project, 
                       thresh = 0.5, 
                       p_method = "fdr", 
                       method="diem"){
    require(diem)

	dp <- paste0("data/processed/", dir_label, "/", method, "/")
	dir.create(dp, recursive=TRUE, showWarnings=FALSE)
	dr <- paste0("results/", dir_label, "/", method, "/")
	dir.create(dr, recursive=TRUE, showWarnings=FALSE)
	dir_plot <- paste0(dr, "plots/")
	dir.create(dir_plot, recursive=TRUE, showWarnings=FALSE)

    fn <- paste0(dp, project, ".diem_sce.rds")

    sce <- assign_clusters(sce)
    sce <- estimate_dbr_score(sce, p_method = p_method)
    sce <- call_targets(sce, thresh = thresh)

	# Save
	saveRDS(sce, paste0(dp, project, ".diem_sce.rds"))

    mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
    genes <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x = sce, genes=genes, name="MALAT1")


	# Plots
	prefix=paste0(dir_plot, project)
	pdf(paste0(prefix, ".umi_gene.call.pdf"), width=9)
	plot_umi_gene_call(sce)
	dev.off()

	# Plots
    # Plot clusters
    pdfname <- paste0(dir_plot, project, ".umi_gene.call.pdf")
    jpgname <- paste0(dir_plot, label, "umi_gene.call.jpeg")
    pdf(pdfname, width=9,height=9)
    plot_umi_gene_call(sce)
    dev.off()
    system(paste("convert", "-density", "200", pdfname, jpgname))

    pdfname <- paste0(dir_plot, project, ".umi_gene.mtpct.pdf")
    jpgname <- paste0(dir_plot, label, "umi_gene.mtpct.jpeg")
    pdf(pdfname, width=9,height=9)
    plot_umi_gene(sce, color="pct.mt", color_name="MT%")
    dev.off()
    system(paste("convert", "-density", "200", pdfname, jpgname))

    pdfname <- paste0(dir_plot, project, ".umi_gene.malat1.pdf")
    jpgname <- paste0(dir_plot, label, "umi_gene.malat1.jpeg")
    pdf(pdfname, width=9,height=9)
    plot_umi_gene(sce, color="MALAT1", color_name="MALAT1")
    dev.off()
    system(paste("convert", "-density", "200", pdfname, jpgname))

	return(sce)
}
