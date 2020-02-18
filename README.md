
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

## Spliced reads fractions

I ran velocyto on the DiffPA and adipose tissue, and copied the results in. 
The script `process_sf.R` processes the raw data copied over into an 
analysis-ready table.

I ran velocyto on the mouse brain data here. The veloctyo results from this 
are processed with `process_sf.R`.

```bashrc
cd scripts/mouse_nuclei_2k
Rscript get_top_ids.R
qsub -sync y run10x.sh
python3.6 get_splice_frctn.py
cd ../

Rscript process_sf.R
cd ../
```

## Testing parameters

I tested different k_init values with these scripts

```bashrc
cd scripts/atsn
qsub -sync y diem.k_init.t0.sh
qsub -sync y diem.k_init.sh
cd ../../
```

I tested filtering values with these scripts:

```bashrc
cd scripts/adpcyte
qsub -sync y diem.filt.sh
cd ../
cd mouse_nuclei_2k
qsub -sync y diem.filt.sh
cd ../
cd atsn
qsub -sync y diem.filt.sh
cd ../../
```

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

# Integrated:
qsub -sync y diem.integrated.fltr.1.R
qsub -sync y diem.integrated.sh
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

For the PBMC single-cell data
```bashrc
cd scripts/fresh_68k/
qsub -sync y diem.sh
qsub -sync y emptydrops.sh
qsub -sync y quantile.sh
cd ../../
```

Analysis and figures
```bashrc
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
```

## Run everything at once

You can run all the analyses with the following script.
This requires qsub and you may have to change the headers. 

For the `run_all.sh` script, you will have to change or comment out 
the bash variable assignments since they are specific to my system 
and will cause your run to fail.

```bash
./run_all.sh
```

or better yet

```bash
qsub run_all.sh
```

