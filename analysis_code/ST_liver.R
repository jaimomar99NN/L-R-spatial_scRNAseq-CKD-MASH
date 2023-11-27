library(tidyverse)
library(Seurat)


data.dir = "" 
metadata_ST = "" 
MASH_files = list.files(data.dir)

library(glue)
for (sample in MASH_files) {
  folder = paste(data.dir,sample, sep = "")
  file = glue("{sample} filtered_feature_bc_matrix.h5")
  assign(paste("MASH", sample, sep='_'),Load10X_Spatial(data.dir = folder, filename=file))
}

list_MASH <- list(MASH_1A, MASH_1B, MASH_1C, MASH_1D, MASH_2A, MASH_2B, MASH_2C, MASH_2D, MASH_3A, MASH_3B, MASH_3C, MASH_3D, MASH_4A, MASH_4B, MASH_4C,MASH_4D,MASH_5A,MASH_5B,MASH_5C,MASH_5D,MASH_6A, MASH_6B,MASH_6C,MASH_6D)

names(list_MASH) <- MASH_files

meta_MASH <- read_csv(metadata_ST)

#add metadata
for (i in 1:length(list_MASH)) {
  list_MASH[[i]]$disease_status <- meta_MASH %>% filter(linkID_spatial_ffpe == names(list_MASH[i])) %>% select("disease_status") %>% as.factor()
  list_MASH[[i]]$fibro_12_34  <- meta_MASH %>% filter(linkID_spatial_ffpe ==names(list_MASH[i])) %>% select("fibrosis_group_12_34") %>% as.factor()
  list_MASH[[i]]$orig.ident <- names(list_MASH[i]) %>% as.factor()
  list_MASH[[i]] <- SetIdent(list_MASH[[i]], value = "orig.ident")
}

#QC
for (i in 1:length(list_MASH)) { 
  list_MASH[[i]] <- list_MASH[[i]][, list_MASH[[i]]$nFeature_Spatial > 250 & list_MASH[[i]]$nCount_Spatial > 500]
}

#normalization
for (i in 1:length(list_MASH)) {
  list_MASH[[i]] <- SCTransform(list_MASH[[i]], assay = "Spatial", verbose = FALSE, method = "poisson")
  list_MASH[[i]] <- RunPCA(list_MASH[[i]], assay = "SCT", verbose= FALSE)
}

#integration
for (i in 1:length(list_MASH)) {
  list_MASH[[i]] <- SCTransform(list_MASH[[i]], assay = "Spatial", verbose = FALSE, method = "poisson")
}

integ_features <- SelectIntegrationFeatures(object.list = list_MASH, nfeatures = 3000) 

# Merge normalized samples
merge_MASH <- merge(x = list_MASH[[1]],
                    y = list_MASH[2:length(list_MASH)],
                    merge.data = TRUE)

DefaultAssay(merge_MASH) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(merge_MASH) <- integ_features

# Calculate PCs using variable features determined by SCTransform (3000 by default)
merge_MASH <- RunPCA(merge_MASH, assay = "SCT")

library(harmony)
harmonized_MASH <- RunHarmony(merge_MASH, 
                              group.by.vars = c("orig.ident", "fibro_12_34"), 
                              reduction = "pca", assay.use = "SCT", reduction.save = "harmony")