library(rhdf5)
library(SeuratDisk)
library(Seurat)
library(dplyr)

pathh5seurat = "" 
path_metadata = "" 

sn  <-  Connect(pathh5seurat)
sn = as.Seurat(sn, slot = "counts", assay = "RNA")

metadata <- read.csv(path_metadata) %>%
  rename(specimen_id = Participant.ID) %>%
  dplyr:: filter(specimen_id %in% sn@meta.data$specimen_id)


metadata <- inner_join(sn@meta.data, metadata, by = "specimen_id")

metadata$log10genesperUMI <-  log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)

metadata$Tissue.Type[which(str_detect(metadata$Tissue.Type, "Healthy Reference"))] <- "Healthy"

sn@meta.data <- metadata
rownames(sn@meta.data) <- Cells(sn)

#filtering
sn_filtered <- subset(x = sn, 
                      subset= nCount_RNA > 500 &
                        nFeature_RNA > 250 &
                        log10genesperUMI > 0.80 &
                        percent.mt < 20)