# Filter out cells with high mitochondrial gene reads, and cells with sgRNAs
# with low representation

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)

# Import Seurat Object ----------------------------------------------------

tcells <- LoadH5Seurat("data/tcells_all_unfilt.h5Seurat")


# Plot RNA feature metrics and filter -------------------------------------

# store mitochondrial percentage in object meta data
tcells <- PercentageFeatureSet(tcells, 
                               pattern = "^MT-", 
                               col.name = "percent.mt")
head(tcells)

# store ribosomal percentage in object meta data
tcells <- PercentageFeatureSet(tcells, 
                               pattern = "^RP[SL]", 
                               col.name = "percent.ribo")
head(tcells)

# Plot statistics

VlnPlot(tcells, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 2,
        group.by = "condition",
        pt.size = 0) +
  patchwork::plot_annotation(title = "Pre-Filtering")


# Filter out dead cells (high MT) and multiplets (high N_features)
tcells <- subset(tcells, 
                 subset = nFeature_RNA > 400 & 
                   nFeature_RNA < 6000 & 
                   percent.mt < 25) 

# Plot statistics after filtering
VlnPlot(tcells, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 2,
        group.by = "condition",
        pt.size = 0) +
  patchwork::plot_annotation(title = "Post-Filtering")


# Plot CRISPR cell calls and filter ---------------------------------------

# Plot n cells per guide target 
guide_target_summary <- tcells@meta.data %>%
  group_by(condition, gene) %>%
  summarise(n_cells = n())
guide_target_summary

guide_target_summary %>%
  arrange(n_cells) %>%
  mutate(gene = as.character(gene) %>% as_factor()) %>%
  ggplot(aes(x = gene, y = n_cells, fill = condition)) +
  geom_col(position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Remove cells with guides targeting genes w < 100 cells in either group:
## IRX4, PRDM1, TCF7, and HELZ2

genes2filter <- guide_target_summary %>%
  dplyr::filter(n_cells < 100) %>%
  pull(gene) %>%
  unique() %>%
  as.character()
genes2filter

genes2keep <- guide_target_summary %>%
  dplyr::filter(!gene %in% genes2filter) %>%
  pull(gene) %>%
  unique() %>%
  as.character()
genes2keep

tcells <- subset(tcells, subset = gene %in% genes2keep) 
tcells

# Save Seurat Object ------------------------------------------------------

SaveH5Seurat(tcells, filename = "data/tcells_all_filt.h5Seurat")



