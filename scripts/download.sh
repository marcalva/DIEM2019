#!/bin/bash

cd ../

mkdir -p data/raw
cd data/raw

# Download 10X adult brain snRNA-seq data
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/nuclei_2k/nuclei_2k_raw_gene_bc_matrices.tar.gz
tar xzf nuclei_2k_raw_gene_bc_matrices.tar.gz
rm nuclei_2k_raw_gene_bc_matrices.tar.gz
mkdir -p mouse_nuclei_2k
rm -rf mouse_nuclei_2k/*
mv -f raw_gene_bc_matrices mouse_nuclei_2k

# Download BAMs
cd mouse_nuclei_2k

mkdir -p bams
cd bams
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/nuclei_2k/nuclei_2k_possorted_genome_bam.bam
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/nuclei_2k/nuclei_2k_possorted_genome_bam_index.bam.bai
cd ../

# mouse ref
mkdir -p ref
cd ref
wget ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz
gunzip Mus_musculus.GRCm38.84.gtf.gz
cd ../

cd ../
# Download 68K
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_raw_gene_bc_matrices.tar.gz
tar xzf fresh_68k_pbmc_donor_a_raw_gene_bc_matrices.tar.gz
rm fresh_68k_pbmc_donor_a_raw_gene_bc_matrices.tar.gz
mkdir -p fresh_68k
rm -rf fresh_68k/*
mv -f matrices_mex fresh_68k

