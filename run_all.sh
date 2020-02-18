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

# Test k_init
cd atsn
qsub -sync y diem.k_init.t0.sh
qsub -sync y diem.k_init.sh
cd ../

# Test filter values for k_init = 30
cd adpcyte
qsub -sync y diem.filt.sh
cd ../
cd mouse_nuclei_2k
qsub -sync y diem.filt.sh
cd ../
cd atsn
qsub -sync y diem.filt.sh
cd ../


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
# Rscript diem.seurat_intg.R
qsub -sync y diem.integrated.fltr.1.R
qsub -sync y diem.integrated.sh
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
cd scripts/atsn
Rscript plot.k_init.t0 # Figure S3
cd ../../

cd scripts/
Rscript plot_sf_relat.R # Figure S1
Rscript plot_quant.fig1.R # Figure 1 and S2
Rscript plot_de_ct.R # Figure 2
Rscript plot_cor_log2fc.R # Figure S3
Rscript plot_umi_ngene.call.R # Fig 3B
cd atsn
    Rscript diem.seurat.k_init.t0.R # Figure S4c,d
cd ../
Rscript plot_params.R # Figure S5
Rscript plot_diem_clusters.R # Figure S6
Rscript plot_removed.R # Figure S7
Rscript plot_overlap_clusters.R # Figure S8 S10
Rscript diffpa_malat1.R # Figure S9
Rscript plot_subtype.R # Figure S11
Rscript mt_boxplot.R # Figure S12
Rscript plot_umap_intg.R # Figure S13
Rscript plot_init_filter_umi.R # Figure S4

Rscript plot_umap_sf.R # Figure 4a,b
Rscript overlap_droplets.R # Figure 4c
Rscript plot_umap_mt.R # Figure 5
Rscript plot_fresh_68k.R # Figure 6



