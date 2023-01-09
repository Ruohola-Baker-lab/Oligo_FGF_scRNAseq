suppressPackageStartupMessages({
  library(monocle3)
  library(tidyr)
  library(ggplot2)
  setwd("~/Documents/FGF_iEndo_1")
  source("chiSq_test_functions.R")
  
})

cds = readRDS("cds.RDS")


rna.umis.per.cell = 
  read.table(file = "UMIs.per.cell.barcode", sep = "\t",
             header = F,
             col.names = c("Sample","Cell","n.rna.umi"),
             colClasses = c("character","character","numeric"))
  
rna.umis.per.cell$Cell = 
  stringr::str_remove_all(rna.umis.per.cell$Cell," ")

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
                             downsample_rate=1)

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


cds_fgf = 
  new_cell_data_set(expression_data = counts(cds)[,cells.to.keep],
                    cell_metadata = new.coldata[cells.to.keep,],
                    gene_metadata = rowData(cds))

rm(list = ls()[ls() != "cds_fgf" ])

cds_fgf = 
  cds_fgf %>%
  estimate_size_factors() %>%
  preprocess_cds() %>%
  reduce_dimension()

colData(cds_fgf)$umap1 = reducedDim(cds_fgf,type = "UMAP")[,1]
colData(cds_fgf)$umap2 = reducedDim(cds_fgf,type = "UMAP")[,2]

cd.for.plot = 
  colData(cds_fgf) %>% as.data.frame()

ggplot() +
  geom_point(data = cd.for.plot %>%
               dplyr::select(-Day),
             aes(x = umap1, 
                 y = umap2),
             color = "grey80",
             size = 0.50,
             stroke = 0) +
  geom_point(data = cd.for.plot,
             aes(x = umap1, 
                 y = umap2,
                 color = Oligo),
             size = 0.5,
             stroke = 0) +
  theme_classic() +
  facet_wrap(~Day)+
  theme(legend.position = "none")
ggsave(filename = "umap.by.day.png",
       dpi = 300,
       height = 6,
       width = 6)


saveRDS("cds_fgf",file = "fgf_cds2.RDS")

