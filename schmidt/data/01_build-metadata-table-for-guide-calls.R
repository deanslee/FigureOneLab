# Build cell metadata table using CellRanger guide calls

# Load packages -----------------------------------------------------------

library(tidyverse)

# Import data -------------------------------------------------------------

# Guide calls from CellRanger
guide_calls <- read_tsv("data/cellranger-guidecalls-aggregated-unfiltered.txt")
guide_calls


# Build metadata table ----------------------------------------------------

metadata <- guide_calls %>%
  # Filter for guide singlets
  dplyr::filter(num_features == 1) %>% 
  # Filter for number of UMIs > 5
  mutate(num_umis = as.double(num_umis)) %>% 
  dplyr::filter(num_umis >= 5) %>% 
  # Add condition, guide target ("gene"), perturbed/control (crispr) columns 
  mutate(condition = if_else(str_detect(cell_barcode, "-1|-2|-3|-4"), 
                             "Nostim", 
                             "Stim"),
         gene = str_sub(feature_call, end = -3),
         crispr = if_else(str_detect(gene, "NO-TARGET"), "NT", "perturbed")) %>%
  select(cell_barcode, condition, crispr, guide_id = feature_call, gene)
metadata


# Save metadata table -----------------------------------------------------

write_tsv(metadata, "data/cell_metadata.txt")



