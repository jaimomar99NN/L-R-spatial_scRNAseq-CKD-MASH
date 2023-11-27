library(Seurat)
library(tidyverse)

sn_rds_liver = "" 

sn <- readRDS(sn_rds_liver)

annotation <- read.csv(tsv_annotation, sep = "\t")
meta <- sn@meta.data %>% rownames_to_column("barcode")

sn$log10genesperUMI <- log10(metadata$nGenes) / log10(metadata$nCounts)

sn_filtered <- subset(x = sn, 
                           subset= nCounts > 500 &
                        nGenes > 250 &
                        log10genesperUMI > 0.80 &
                        pct_mito < 20)