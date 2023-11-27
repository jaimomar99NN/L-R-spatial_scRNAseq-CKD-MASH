library(tidyverse)
library(Seurat)

CKD.dir <- ""
metadata_path <- "" 
files <- list.files(CKD.dir)

#read data
CKD_27_10066 = dior::read_h5(file=paste(CKD.dir,"90ac4420-ff53-47dc-acfc-829339a38688_12874867-d3ee-4de2-a474-8995c0f1029b_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_28_12265 = dior::read_h5(file=paste(CKD.dir,"456b07c1-cca6-40e0-8c9c-15029342191e_ef6fc180-8ecd-404f-ba44-d34b6cb74a57_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")
#two samples with ST_28_12265 ID, taking the biggest file

CKD_29_10012 = dior::read_h5(file=paste(CKD.dir,"fbdc27ed-a096-4d4b-8b22-187f6b6403d2_a241e119-4b97-4f65-b746-1904944e6bb6_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_29_10013 = dior::read_h5(file=paste(CKD.dir,"f45233d8-a03b-4f52-87d9-2139b3054ddd_e24e777e-ccc5-4e84-a2af-9b6a4c987934_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_29_10280 = dior::read_h5(file=paste(CKD.dir,"bf4a4ea9-ff59-4255-8d2b-ab833e9902d1_c60c8e3d-afb3-4de4-b9a1-9ca8b2025201_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_29_10282 = dior::read_h5(file=paste(CKD.dir,"1dc4a172-07ee-4a50-8e29-b8dd3f56465c_65cb3483-0b4a-4344-af7b-5d13c93191f7_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_29_10404 = dior::read_h5(file=paste(CKD.dir,"03d88385-cbf5-4e65-8964-777c72c7dad7_53893aab-179d-416d-9922-60b0ab69dfa4_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_30_10125 = dior::read_h5(file=paste(CKD.dir,"06e0f003-c04c-4ccd-a465-c4395839ab55_68c30f3a-ad68-4948-9270-32271f2a8185_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_30_10631 = dior::read_h5(file=paste(CKD.dir,"c6ac867d-531a-43ac-8ca2-1aff47c58dc2_6188d61d-d691-4d34-9add-4a7b9020ec85_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_30_10929 = dior::read_h5(file=paste(CKD.dir,"760da2e8-2087-42cf-9ea4-d6fdaa550628_816326d7-3315-48af-9387-4301c8acb690_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_31_10042 = dior::read_h5(file=paste(CKD.dir,"a81197e1-10ab-46be-a3f1-60d4e06c97dd_790af35b-462f-4daf-b756-7e3822253858_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_31_10221 = dior::read_h5(file=paste(CKD.dir,"27e4ddca-accc-4e3b-a835-e83f45ae99b3_591bd1a8-085a-4519-a3fe-0c243de5b72b_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_32_10003 = dior::read_h5(file=paste(CKD.dir,"48d2f1b2-3461-408b-92a2-0a879d74321e_6158725f-5ef8-445e-9a4b-4fb9f94e2f48_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_32_10074 = dior::read_h5(file=paste(CKD.dir,"742e3c6c-5f0f-4e6f-af36-578b9df1349b_c3d70bd4-a9c6-4a39-a9d5-ec0bb205edad_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

CKD_33_10331 = dior::read_h5(file=paste(CKD.dir,"b1cbf3fa-b8bc-48ac-af2b-8d5b4c0748eb_0643f336-2715-4cf5-98b3-a2603e560875_spatial-transcriptomics.h5ad.h5",sep=""), target.object = 'seurat', assay.name = "spatial")

meta_CKD <- read.csv(metadata_path, sep =",")

#create list objects
list_CKD <- list(CKD_27_10066, CKD_28_12265, CKD_29_10012, CKD_29_10013, CKD_29_10280, CKD_29_10282, CKD_29_10404, CKD_30_10125, CKD_30_10631, CKD_30_10929, CKD_31_10042, CKD_31_10221, CKD_32_10003, CKD_32_10074, CKD_33_10331)

CKD_ids <- list("27_10066", "28_12265", "29_10012", "29_10013", "29_10280", '29_10282', "29_10404", "30_10125", "30_10631", "30_10929", "31_10042", "31_10221", "32_10003", "32_10074", "33_10331")

names(list_CKD) <- CKD_ids
names(CKD_ids) <- CKD_ids

#QC
for (i in 1:length(list_CKD)) {
  list_CKD[[i]]$orig.ident <- names(list_CKD[i]) %>% as.factor()
  rna_counts <- GetAssay(list_CKD[[i]],"Spatial")
  list_CKD[[i]]$nCount_Spatial <- colSums(rna_counts) #ncounts: number of transcripts for each spot within the tissue sample.
  list_CKD[[i]]$nFeature_Spatial<- colSums(rna_counts@counts != 0) #number of genes
  list_CKD[[i]] <- PercentageFeatureSet(list_CKD[[i]], "^MT-", col.name = "percent_mito")
  list_CKD[[i]] <- PercentageFeatureSet(list_CKD[[i]], "(^HBA1)|(^HBA2)|(^HBB)|(^HBD)|(^HBE1)|(^HBG1)|(^HBG2)|(^HBM)|(^HBQ1)|(^HBZ)", col.name = "percent_hb")
  
  #assign state disease
  name_id <-  names(list_CKD[i]) %>% str_replace("_", "-" )
  list_CKD[[i]]$disease_status <- meta_CKD %>% filter(Participant.ID == name_id ) %>% select("Tissue.Type")  %>% as.factor()
}

for (i in 1:length(list_CKD)) {
  list_CKD[[i]] <- list_CKD[[i]][, list_CKD[[i]]$nFeature_Spatial > 250 & list_CKD[[i]]$nCount_Spatial > 500 & list_CKD[[i]]$percent_mito < 20]
}

#save spatial information
library(jsonlite)

for (i in 1:length(list_CKD)) {
  sample <- names(list_CKD[[i]])
  dir.create(paste0(path_sample, sample))
  tissue_position_list <- list_CKD[[i]]@images[[1]]@coordinates
  print(tissue_position_list)
  write.table(tissue_position_list, paste0(path_sample, sample, "/tissue_position_list.csv"), col.names=FALSE, sep=",")
  
  json <- list(
    "tissue_hires_scalef"      = unbox(list_CKD[[i]]@images[[1]]@scale.factors$hires),
    "tissue_lowres_scalef"     = unbox(list_CKD[[i]]@images[[1]]@scale.factors$lowres),
    "fiducial_diameter_fullres"= unbox(list_CKD[[i]]@images[[1]]@scale.factors$fiducial),
    "spot_diameter_fullres"    = unbox(list_CKD[[i]]@images[[1]]@scale.factors$spot)
    
  )
  
  write_json(json, paste0(path, sample, "/scalefactors_json.json"))
}


#normalization
for (i in 1:length(list_CKD)) {
  list_CKD[[i]] <- SCTransform(list_CKD[[i]], assay = "Spatial", verbose = FALSE, method = "poisson")
  list_CKD[[i]] <- RunPCA(list_CKD[[i]], assay = "SCT", verbose= FALSE)
}

#integration
integ_features <- SelectIntegrationFeatures(object.list = list_CKD, nfeatures = 3000) 

# Merge normalized samples
merge_CKD <- merge(x = list_CKD[[1]],
                   y = list_CKD[2:length(list_CKD)],
                   merge.data = TRUE)

DefaultAssay(merge_CKD) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(merge_CKD) <- integ_features

# Calculate PCs using variable features determined by SCTransform (3000 by default)
merge_CKD <- RunPCA(merge_CKD, assay = "SCT")

library(harmony)
harmonized_CKD <- RunHarmony(merge_CKD, 
                             group.by.vars = c("orig.ident", "disease_status"), 
                             reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

