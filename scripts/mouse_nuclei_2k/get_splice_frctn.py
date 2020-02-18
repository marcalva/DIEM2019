# Get spliced ratios per droplet

import os
import loompy
import numpy as np
import pandas as pd

os.chdir("../../")

fn = "data/processed/mouse_nuclei_2k/velocyto/mouse_nuclei_2k.loom"

ds = loompy.connect(fn)
gene_ids = ds.ra['Gene']
mt_genes = ["MT-" in gene for gene in gene_ids]
keep = ["MT-" not in gene for gene in gene_ids]
spliced = ds.layers['spliced'].sparse()
sa = spliced.toarray()
unspliced = ds.layers['unspliced'].sparse()
ua = unspliced.toarray()

ua = ua[keep,:]
sa = sa[keep,:]

spliced_d = sa.sum(axis = 0)
unspliced_d = ua.sum(axis = 0)
all_d = spliced_d + unspliced_d
splice_f = spliced_d / all_d
splice_f = splice_f.transpose()
splice_f = np.array(splice_f)

# Update cell IDs from
# mouse_nuclei_2k:AAACCTGGTAAATGTGx
# to
# mouse-nuclei_2k_AAACCTGGTAAATGTG
cell_ids = ds.ca['CellID']
cell_ids = np.char.replace(cell_ids, "mouse_nuclei", "mouse-nuclei")
cell_ids = np.char.replace(cell_ids, ":", "_")
cell_ids = np.char.replace(cell_ids, "x", "")
df = pd.DataFrame(splice_f, index = cell_ids)
df.columns = pd.Series("SpliceFrctn")

ofn = "data/raw/splice_frctn/mouse_nuclei_2k.splice_fraction.txt"
df.to_csv(ofn)
ds.close()

