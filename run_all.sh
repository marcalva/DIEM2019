#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4000M,h_rt=12:00:00,highp
#$ -v QQAPP=openmp
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -V
#$ -o run_all.sh.log

. /u/local/Modules/default/init/modules.sh
module load R/3.5.1

#==================================================
# IMPORTANT
# These commands are system-dependent and required 
# to run the analysis on my system.
# 
# You may not need these and can comment them out
# 
#==================================================

export HOME=/u/project/pajukant/malvarez/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/local/compilers/gcc/4.7.2/lib:/u/local/compilers/gcc/4.7.2/lib64:/u/project/pajukant/malvarez/lib/
export R_LIBS_USER=/u/project/pajukant/malvarez/lib/R_%V:R_LIBS_USER

#==================================================
#==================================================


# Download
cd scripts
./download.sh

# spliced reads fraction
cd mouse_nuclei_2k
Rscript get_top_ids.R
qsub -sync y run10x.sh
python3.6 get_splice_frctn.py
cd ../

Rscript process_sf.R
Rscript get_sf_midpoint.R

# Test k_init values 1,5,10,...,50
cd atsn
qsub -sync y diem.k.R
cd ../

# Run filtering and clustering
# Adipocyte data
cd adpcyte
qsub -sync y diem.R
Rscript emptydrops.R
Rscript quantile.R
cd ../

# mouse brain
cd mouse_nuclei_2k
qsub -sync y diem.R
Rscript emptydrops.R
Rscript quantile.R
cd ../

# adipose tissue
cd atsn
qsub -sync y diem.R
Rscript diem.seurat.R
diem.integrated.R
qsub -sync y emptydrops.sh
Rscript emptydrops.seurat.R
./quantile.sh
Rscript quantile.seurat.R

Rscript diem.seurat.res2.R
Rscript emptydrops.seurat.res2.R
Rscript quantile.seurat.res2.R

# Debris for adipose
Rscript diem.debris.seurat.R
Rscript emptydrops.debris.seurat.R
Rscript quantile.seurat.R

# DE for adipose
Rscript edgeR_ct.R
Rscript edgeR_de.R
cd ../../

# Fresh 68K PBMCs
cd scripts/fresh_68k/
qsub -sync y diem.sh
qsub -sync y emptydrops.sh
qsub -sync y quantile.sh
cd ../../

cd scripts/

# Plots, quants, stats
Rscript melt_sf_data.R # Place data in one file

Rscript plot_barcode_rank.R # Fig 1a barcode rank plot
Rscript plot_quant_ndrop.bar.R # Fig 1b
Rscript plot_quant_umap_sf.R # Fig 1c UMAP %SF 
Rscript plot_quantile_sf_barplot.R # Fig 1d barplot cluster %SF
Rscript plot_sf_relat.R # Figure S1
Rscript plot_quant.ptct_mt.R # Figure S2
Rscript plot_de_ct.R # Figure 2
Rscript plot_cor_log2fc.R # Figure S3
Rscript plot_k.R # Figure S4 test range of k values
Rscript plot_cluster_sf_ngenes.R # Figure S5 n_genes by average %SR of clusters
Rscript plot_thresh.R # Figure S6 showing effect of varying t threshold
Rscript plot_diem_clusters.R # Figure S7 showing overlap between diem and seurat clusters
Rscript plot_drop.sf.score.R # Fig 3b
Rscript plot_umi_ngene.call.R # Fig 3c
Rscript plot_method.sample.n_bg_nuc.R # Fig 4a, 5a
Rscript plot_method.droplet_sf.R # Figure 4b, 5b
Rscript plot_pass_clust.R # Figure 4c-e
Rscript plot_cluster_sf_barplot.R # Figure S8a S11a cluster SF barplots
Rscript plot_overlap_heatmap.R # Figure S8b,d S11b cluster overlaps
Rscript diffpa_malat1 # Figure S9
Rscript plot_gene_ct.boxplot.R # Figure S10 boxplot of cluster expression of cell type markers
Rscript plot_debris_clust.R # Figure 5c-e
Rscript plot_removed.sf_clust.genes.R # Figure S12 boxplot of removed clusters of cell type markers
Rscript plot_score_umi.R # Figure S13
Rscript plot_fresh_68k.R # Figure 6
Rscript plot_umap_intg.R # Figure S14



