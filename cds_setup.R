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

# Read unprocessed CDS
cds <- readRDS(file = "cds.RDS")

# Read UMI per cell information
rna.umis.per.cell = 
  read.table(file = "UMIs.per.cell.barcode", sep = "\t",
             header = F,
             col.names = c("Sample","Cell","n.rna.umi"),
             colClasses = c("character","character","numeric"))

rna.umis.per.cell$Cell = 
  stringr::str_remove_all(rna.umis.per.cell$Cell," ")

# Create spatial hash table using oligo assignments
hashtable = 
  read.table(file = "hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("Sample","Cell","Oligo","Axis","Count"),
             colClasses = c("character","character","character","numeric","numeric"))

hashtable =
  hashtable %>% 
  left_join(rna.umis.per.cell %>% dplyr::select(-Sample))

hashtable$in_cds = hashtable$Cell %in% colnames(cds)

# Correct the spatial hash table by removing background estimated from debris
hashtable$Oligo = 
  factor(hashtable$Oligo, 
         levels = unique(hashtable$Oligo))

# Background Cells fall below RNA treshold
background_cells = 
  hashtable %>%
  filter(n.rna.umi < 10) %>%
  pull(Cell) %>%
  as.character() %>%
  unique()

# Test Cells that are in the CDS object
test_cells =
  hashtable %>%
  filter(in_cds) %>%
  pull(Cell) %>%
  as.character() %>%
  unique()

background_sample_hashes = 
  hashtable %>% 
  filter(Cell %in% background_cells) 

test_cell_hashes = 
  hashtable %>% 
  filter(Cell %in% test_cells) 

hash_df = assign_hash_labels(test_cell_hashes, 
                             background_sample_hashes, 
                             downsample_rate = 1)

corrected_hash_table = 
  left_join(hashtable,
            hash_df,
            by = c("Cell")) 


passing.cells = 
  corrected_hash_table %>%
  filter(in_cds,
         top_to_second_best_ratio > 4) %>%
  pull(Cell) %>%
  unique()

length(passing.cells)

new.coldata = 
  colData(cds[,passing.cells]) %>%
  as.data.frame()

new.coldata$Oligo = 
  left_join(new.coldata,
            corrected_hash_table %>%
              filter(in_cds,
                     top_to_second_best_ratio > 4) %>%
              group_by(Cell) %>%
              dplyr::top_n(1,Count) %>%
              dplyr::select(Cell,top_oligo),
            by = "Cell") %>%
  pull(top_oligo)

sample.key = 
  read.csv(file = "SamplesKey.csv") %>%
  dplyr::select(-X)

new.coldata = 
  left_join(new.coldata,sample.key)
rownames(new.coldata) = new.coldata$Cell

cells.to.keep = 
  new.coldata %>%
  filter(!is.na(Day)) %>%
  pull(Cell)

cds = 
  new_cell_data_set(expression_data = counts(cds)[,cells.to.keep],
                    cell_metadata = new.coldata[cells.to.keep,],
                    gene_metadata = rowData(cds))

rm(list = ls()[ls() != "cds" ])

saveRDS(cds, file = "cds_all_raw.RDS")

cds_myod <- cds[,colData(cds)$Tube %in% c("33A", "30", "31", "33B", "32B", "34A", "32A", "34B")]
cds_fgf <- cds[,!(colData(cds)$Tube %in% c("33A", "30", "31", "33B", "32B", "34A", "32A", "34B"))]
saveRDS(cds_myod, file = "cds_myod_raw.RDS")
saveRDS(cds_fgf, file = "cds_fgf_raw.RDS")
