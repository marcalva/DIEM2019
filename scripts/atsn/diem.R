
setwd("../../")

library(diem)
source("scripts/common/diem_pipe.R")
source("scripts/common/standard_seurat.R")

#=========================================
# Set variables
#=========================================

label <- "atsn"
method <- "diem"
dir10X <- c("~/sn_rnaseq2018/data/processed/TwinFatSN/110686_121444/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
			"~/sn_rnaseq2018/data/processed/TwinFatSN/54172_58356/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
			"~/sn_rnaseq2018/data/processed/TwinFatSN/57289_64605/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
			"~/sn_rnaseq2018/data/processed/TwinFatSN/62630_75306/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
			"~/sn_rnaseq2018/data/processed/TwinFatSN/63099_76245/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/",
			"~/sn_rnaseq2018/data/processed/TwinFatSN/63099_76246/GRCh38_primary_assembly_pli_premrna_genome/outs/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/")
snames <- c("110686_121444", "54172_58356", "57289_64605", "62630_75306", "63099_76245", "63099_76246")
lab_ids <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

#=========================================
#=========================================

args = commandArgs(TRUE)
i = as.integer(args[1])

counts <- diem::read_10x(dir10X[i])
diem_pipe_out <- diem_pipe(counts, dir_label=label, project=lab_ids[i])
