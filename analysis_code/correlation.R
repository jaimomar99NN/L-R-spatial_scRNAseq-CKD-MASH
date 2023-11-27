library(Seurat)
library(corrplot)

ST_liver <- readRDS(path_ST_list_liver)

markers <-  c("LUM", "COL1A1", "PTGDS", "IGFBP7", "ACTA2")
ligand <- "PDGFC"
receptor <- "PDGFRA"
genes <- c(ligand,receptor,markers)

# Filter the list based on the disease
Healthy_list <- Filter(function(x) "Healthy" %in% levels(x$disease_status), ST_liver)
MASH_list <- Filter(function(x) "MASH" %in% levels(x$disease_status), ST_liver)

#CALCULATING CORRELATION TAKING EXPRESSION ALL SPOTS
correlation_spatial <- function(genes, ST_list){
  co_expression_df <- data.frame(matrix(nrow = length(genes), ncol = length(genes)), row.names=genes)
  colnames(co_expression_df) <- genes
  
  combinations <- combn(genes,2)
  
  for (i in 1:ncol(combinations)) {
    pair_genes <- combinations[,i]
    cor_v <- numeric(length(ST_list))
    
    for(j in 1:length(ST_list)){
      data_j <- ST_list[[j]]@assays$SCT@data[pair_genes, ] 
      cor_v[j] <- cor(data_j[1, ], data_j[2, ], method = c("pearson"))
    }
    co_expression_df[pair_genes[1], pair_genes[2]] <- mean(cor_v)
    co_expression_df[pair_genes[2], pair_genes[1]] <- mean(cor_v)
  }
  diag(co_expression_df) <- 1.0
  return(as.matrix(co_expression_df))
}
  
co_expression <- correlation_spatial(genes, ST_liver)
co_expression_Healthy <- correlation_spatial(genes, Healthy_list)
co_expression_MASH <- correlation_spatial(genes, MASH_list)
  