
# DIEM analysis

This repo contains the analysis scripts for evaluating 
the DIEM method.

## Downloading data sets

The public data sets can be downloaded with the following script:
```bashrc
cd scripts
./download.sh
cd ../
```

You can download the 10X CellRanger output for the 
differentiating preadipocytes (DffPA) and adipose tissue 
single-nucleus data sets here: XXX XXX

I placed the 10X CellRanger output in `data/raw/`, where the 
DiffPA data was placed in `adpcyte`, while the 6 adipose tissue 
samples were placed in `AT1` through `AT6`. Each of these 
directories should have the files `barcodes.tsv`, `genes.tsv`, 
and `matrix.mtx`.

## Run filtering methods

Run each of the three filtering for the three methods.

For the DiffPA
```bashrc
cd scripts/adpcyte/
Rscript diem.R
Rscript emptydrops.R
Rscript quantile.R
cd ../../
```

For the mouse brain
```bashrc
cd scripts/mouse_nuclei_2k
Rscript diem.R
Rscript emptydrops.R
Rscript quantile.R
cd ../../
```

For the adipose tissue
```bashrc
cd scripts/atsn
qsub -sync y diem.sh
Rscript diem.seurat.R
qsub -sync y emptydrops.sh
Rscript emptydrops.seurat.R
./quantile.sh
Rscript quantile.seurat.R
cd ../../
```

For running differential expression between nuclear, cell type, 
and debris RNA profiles
```bashrc
cd scripts/atsn
Rscript edgeR_ct.R
Rscript edgeR_de.R
cd ../../
```

For clustering the filtered out droplets in the adipose tissue
```bashrc
cd scripts/atsn
Rscript diem.debris.seurat.R
Rscript emptydrops.debris.seurat.R
Rscript quantile.debris.seurat.R
cd ../../
```

Analysis and figures
```bashrc
cd scripts
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
```

