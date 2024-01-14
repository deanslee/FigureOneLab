
# Generate Seurat object and add metadata

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)

# Import Data -------------------------------------------------------------

# Full aggregated data from 8 chromium Wells
aggr_dat <- Read10X("~/Dropbox (Gladstone)/CRISPRa/Perturb-Seq/A046_CRISPRa-PerturbSeq-LZR06/data/cellranger_aggr/outs/count/filtered_feature_bc_matrix/")
aggr_dat

# Cell metadata (prefiltered guide-calls for singlets)
cell.metadata <- read_tsv("data/cell_metadata.txt", col_types = "cffff")
cell.metadata



# Prepare Cell Metadata ---------------------------------------------------

cell.metadata <- cell.metadata %>%
  as.data.frame() %>%
  column_to_rownames("cell_barcode")
cell.metadata %>% head()


# Prepare gene expression count data --------------------------------------

# Guide singlet cell filter (prefiltered in metadata table)
singlet_barcodes <- rownames(cell.metadata)

aggr_dat$`Gene Expression` <- aggr_dat$`Gene Expression`[,singlet_barcodes]
aggr_dat



# Generate Seurat Object --------------------------------------------------

tcells <- CreateSeuratObject(counts = aggr_dat$`Gene Expression`, 
                             min.cells = 3, 
                             min.features = 200, 
                             meta.data = cell.metadata)
tcells
tcells %>% head()


# Save Data ---------------------------------------------------------------

SaveH5Seurat(tcells, filename = "data/tcells_all_unfilt.h5Seurat")




