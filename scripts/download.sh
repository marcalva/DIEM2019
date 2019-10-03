#!/bin/bash

cd ../data

mkdir -p raw
cd raw

# Download 10X adult brain snRNA-seq data
mkdir -p mouse_nuclei_2k
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/nuclei_2k/nuclei_2k_raw_gene_bc_matrices.tar.gz
tar xzvf nuclei_2k_raw_gene_bc_matrices.tar.gz
rm nuclei_2k_raw_gene_bc_matrices.tar.gz
mv raw_gene_bc_matrices mouse_nuclei_2k

wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/nuclei_2k/nuclei_2k_filtered_gene_bc_matrices.tar.gz
tar xzvf nuclei_2k_filtered_gene_bc_matrices.tar.gz
rm nuclei_2k_filtered_gene_bc_matrices.tar.gz
mv filtered_gene_bc_matrices mouse_nuclei_2k

# Download 68K
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_raw_gene_bc_matrices.tar.gz
tar xzvf fresh_68k_pbmc_donor_a_raw_gene_bc_matrices.tar.gz
rm fresh_68k_pbmc_donor_a_raw_gene_bc_matrices.tar.gz
mkdir -p fresh_68k
rm -rf fresh_68k/*
mv matrices_mex fresh_68k



