
setwd("../../")

source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "quantile"
dir10X <- c("~/sn_rnaseq2018/data/processed/TwinFatSN/110686_121444/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
            "~/sn_rnaseq2018/data/processed/TwinFatSN/54172_58356/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
            "~/sn_rnaseq2018/data/processed/TwinFatSN/57289_64605/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
            "~/sn_rnaseq2018/data/processed/TwinFatSN/62630_75306/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
            "~/sn_rnaseq2018/data/processed/TwinFatSN/63099_76245/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
            "~/sn_rnaseq2018/data/processed/TwinFatSN/63099_76246/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/")
snames <- c("110686_121444", "54172_58356", "57289_64605", "62630_75306", "63099_76245", "63099_76246")
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
names(dir10X) <- lab_ids

#=========================================
#=========================================

# Make directories

	# Create directories
dp <- paste0("data/processed/", label, "/", method, "/")

quantile.l <- lapply(lab_ids, function(x) {
				 readRDS(paste0(dp, x, ".quantile_counts.rds"))})
names(quantile.l) <- lab_ids

countsl <- lapply(lab_ids, function(x) {
                  print(x)
                  counts <- diem::read_10x(dir10X[x])
                  ng <- Matrix::colSums(counts > 0)
                  nc <- Matrix::colSums(counts)
                  keep <- colnames(counts)[ ng >= 200 & nc < quantile.l[[x]]$thresh ]

        return(counts[,keep])
    })

method <- "quantile_debris"
dp <- paste0("data/processed/", label, "/", method, "/")

seurat_pipe_list(countsl, dir_label=label, project.l=lab_ids, project_out=label, method=method, scale.factor=1e3)
