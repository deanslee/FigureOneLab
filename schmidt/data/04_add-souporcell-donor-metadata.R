# Donor demultiplexing metadata from souporcell

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)

# Import Seurat object ----------------------------------------------------

tcells <- LoadH5Seurat("data/tcells_all_filt.h5Seurat")


# Import SoupOrCell Data --------------------------------------------------

souporcell <- tibble(sample = list.files("data/souporcell/")) %>%
  mutate(directory = paste0("data/souporcell/", sample, "/clusters.tsv"),
         data = map(directory, read_tsv))
souporcell

# Fix well identifiers
souporcell <- souporcell %>%
  mutate(well_id = as.character(1:8),
         data = map2(data, 
                     well_id, 
                     function(x,y) 
                       mutate(x, 
                              barcode = str_replace(barcode, ".$", y)
                       )
         )
  )

# Combine to single dataframe

souporcell <- reduce(souporcell$data, bind_rows)
souporcell



# SoupOrCell Donor Values -------------------------------------------------

# Need to harmonize donor values across wells, as "0" or "1" assignment is 
# arbitrary each time. Donors were matched using output .vcf file

souporcell_donor_calls <- read_tsv(
  "data/souporcell_match_vcf_res/donor_calls.txt")
souporcell_donor_calls

# Add to souporcell dataframe
souporcell <- souporcell %>%
  mutate(Well_ID = str_sub(barcode, -1, -1) %>% as.double()) %>%
  left_join(souporcell_donor_calls,
            by = "Well_ID")
souporcell

# Add donor assignments to each droplet
souporcell <- souporcell %>%
  mutate(Donor = if_else(assignment == "0", 
                         Souporcell_call0_DonorA_or_B, 
                         Souporcell_call1_DonorA_or_B),
         Donor = if_else(status == "unassigned",
                         "unassigned",
                         Donor))
souporcell

# Check that all cell barcodes in seurat object are present in souporcell df
all(rownames(tcells@meta.data) %in% souporcell$barcode)

# Add to metadata
souporcell_filtered <- souporcell %>%
  dplyr::filter(barcode %in% rownames(tcells@meta.data))
souporcell_filtered

tcells <- AddMetaData(tcells,
                      souporcell_filtered %>%
                        as.data.frame() %>%
                        column_to_rownames("barcode") %>%
                        select(souporcell_donor = Donor))
tcells %>% head()


# Check male/female donors using RPS4Y1 (Y chromosome gene)
# From donor info sheets donor 1 was female, and donor 2 was male

VlnPlot(tcells, "RPS4Y1", group.by = "souporcell_donor")

# Donor "A" is likely the male, so reassign as "Donor2  

tcells@meta.data <- tcells@meta.data %>%
  mutate(donor = if_else(souporcell_donor == "A", "Donor2", "unassigned"),
         donor = if_else(souporcell_donor == "B", "Donor1", donor)) %>%
  select(-souporcell_donor)

tcells@meta.data %>% head()

# Save Seurat Object ------------------------------------------------------

SaveH5Seurat(tcells,
             "data/tcells_all_filt.h5Seurat",
             overwrite = T)

