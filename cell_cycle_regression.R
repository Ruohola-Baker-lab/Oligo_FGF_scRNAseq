# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(ggridges)
  library(purrr)
  library(monocle3)
  library(RColorBrewer)
  library(stringr)
  library(Seurat)
  
  setwd("~/Documents/FGF_iEndo_1")
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 42)
})

cds_fgf <- readRDS("cds_fgf_raw.RDS")

cds_fgf <- cds_fgf %>%
  estimate_size_factors() %>%
  preprocess_cds(verbose = T) %>%
  reduce_dimension(verbose = T) %>%
  cluster_cells(resolution = 1e-5, verbose = T)

plot_cells(cds_fgf)

# Identify actively cycling cells
plot_cells(cds_fgf, genes = c("TOP2A", "PCNA", "MCM6", "MKI67"), cell_size = 1, label_cell_groups = F)
ggsave("plots_final/Cycling.gene.plot.png", dpi = 600)

# Score cell clusters based on cycling identity
cds_fgf <- cluster_cells(cds_fgf, k = 15, verbose = T)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

normalized_data <- normalized_counts(cds_fgf)

normalized_data@Dimnames[[1]] <- as.character(rowData(cds_fgf)$gene_short_name)
clusters <- monocle3::clusters(cds_fgf)
cell_id <- colnames(cds_fgf)
meta_clusters <- data.frame(row.names=cell_id,clusters)
seurat_object <- CreateSeuratObject(counts = normalized_data, meta.data = meta_clusters)

seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

data_phase  <- seurat_object[[]]
levels(data_phase$Phase) <- c("G1","S","G2M", "Undecided")
ggplot(data = data_phase, aes(x = clusters, fill = Phase)) +
  geom_bar(show.legend = T, position = "fill") + 
  scale_y_continuous(name="Percentage of cells") + 
  theme(text=element_text(size=10, face="bold"),aspect.ratio = 1)
ggsave("plots_final/Cycling.phase.split.png", dpi = 600)

# Transfer cell cycle scores to monocle CDS
colData(cds_fgf)$S_Score  <- data_phase$S.Score
colData(cds_fgf)$G2M_Score  <- data_phase$G2M.Score
colData(cds_fgf)$Phase  <- data_phase$Phase 

# Annotate cycling clusters
colData(cds_fgf)$cluster <- monocle3::clusters(cds_fgf)
coldata <- as.data.frame(colData(cds_fgf))
coldata <- coldata %>%
  mutate(cycling_cells = if_else(cluster %in% c(1,6,7,8,10,15,16,17,18,19,23), "non-cycling", "cycling"))

colData(cds_fgf) <- cbind(colData(cds_fgf), cycling = coldata$cycling_cells)
plot_cells(cds_fgf, color_cells_by = "cycling")
ggsave("plots_final/Cycling.UMAP.pre.regression.png", dpi = 600)

# Regress cell cycle specific differences
cds_fgf <- cds_fgf %>%
  preprocess_cds(verbose = T) %>%
  align_cds(residual_model_formula_str = "~S_Score+G2M_Score", verbose = T) %>%
  reduce_dimension(verbose = T) %>%
  cluster_cells(verbose = T, resolution = 1e-5)

plot_cells(cds_fgf, color_cells_by = "cycling")
ggsave("plots_final/Cycling.UMAP.post.regression.png", dpi = 600)

saveRDS(cds_fgf, file = "cds_fgf_raw.RDS")