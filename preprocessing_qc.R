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
  
  setwd("~/Documents/FGF_iEndo_1")
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 42)
})

cds_fgf <- readRDS("cds_fgf_raw.RDS")

# UMIs per cell
colData(cds_fgf) %>%
  as.data.frame() %>%
  ggplot() +
  geom_boxplot(aes(x = Tube,
                   y = n.umi),
               fill = "grey80") +
  scale_y_log10() +
  theme_classic() +
  ylab("Molecules per Cell") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(axis.title.x = element_blank())
ggsave("plots_final/UMIs.per.cell.png",
       dpi = 600)

# Cells per sample
colData(cds_fgf) %>%
  as.data.frame() %>%
  ggplot() +
  geom_bar(aes(x = Tube),
           fill = "grey90",
           color = "black") +
  scale_y_log10() +
  theme_classic() +
  ylab("Cells per Sample") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(axis.title.x = element_blank())
ggsave("plots_final/Cells.per.sample.png",
       dpi = 600)

# Determine UMI density and filter based on cutoff
qplot(colData(cds_fgf)$n.umi, geom = "density") + 
  scale_x_continuous(limits = c(0, 2500)) + 
  xlab("Number of cells") +
  ylab("UMI Density") +
  theme_classic()
ggsave("plots_final/UMI.density.png", dpi = 600)

cds_fgf <- cds_fgf[,Matrix::colSums(counts(cds_fgf)) < 1500] 

# Mitochondrial gene removal
List_of_MT_genes <- rowData(cds_fgf)$gene_short_name[grep("^MT-", rowData(cds_fgf)$gene_short_name)]
cds_fgf <- cds_fgf[!(rowData(cds_fgf)$gene_short_name %in% List_of_MT_genes),]

saveRDS(cds_fgf, file = "cds_fgf_raw.RDS")
