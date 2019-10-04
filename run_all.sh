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

# Run filtering and clustering
# Adipocyte data
cd adpcyte
Rscript diem.R
Rscript emptydrops.R
Rscript quantile.R
cd ../

# mouse brain
cd mouse_nuclei_2k
Rscript diem.R
Rscript emptydrops.R
Rscript quantile.R
cd ../

# adipose tissue
cd atsn
qsub -sync y diem.sh
Rscript diem.seurat.R
qsub -sync y emptydrops.sh
Rscript emptydrops.seurat.R
./quantile.sh
Rscript quantile.seurat.R

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

# Run analyses for paper
cd scripts/
Rscript plot_quant.fig1.R # Figure 1
Rscript plot_mt_umi.R # Figure S1
Rscript plot_de_ct.R # Figure 2
Rscript plot_cor_log2fc.R # Figure S2
Rscript plot_overview.R # Figure 3A
Rscript plot_umi_ngene.call.R # Fig 3B
Rscript plot_overlap_clusters.R # Figure 4 and Supp Figure S6
Rscript plot_removed.R # Supp Figure S3
Rscript plot_diem_clusters.R # Supp Figure S4
Rscript plot_umap_mt.R # Figure 5
Rscript diffpa_malat1.R # Figure S5
Rscript plot_fresh_68k.R # Figure 6
Rscript mt_boxplot.R # Supp Figure S7
cd ../

