
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

You can download the adipocyte and adipose tissue single-nucleus 
data sets here: XXX XXX

I placed the 10X CellRanger output in 'data/raw/', where the 
DiffPA data was placed in 'adpcyte', while the 6 adipose tissue 
samples were placed in 'AT1' through 'AT6'. Each of these 
directories should have the files 'barcodes.tsv', 

