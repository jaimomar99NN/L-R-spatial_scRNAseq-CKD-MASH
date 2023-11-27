library(Seurat)
library(dplyr)
library(NMF)
library(ggalluvial)
library(CellChat)

sn_MASH_filtered <- readRDS(rds_sn_path_liver)
ST_MASH <- readRDS(rds_ST_path_liver)

Cellchat_pipeline <- function(seurat_object, disease_state_column, condition_status, annotation_column, organism="human", assay= "RNA") {
  Idents(seurat_object) <- disease_state_column # types of disease
  
  subset <- subset(seurat_object, idents = condition_status)
  meta <- subset@meta.data
  
  #create cellchat object
  cellchat <- createCellChat(object = subset, meta = meta, group.by = annotation_column, assay = assay)
  
  #obtain database
  CellChatDB <- paste("CellChatDB", organism, sep=".") %>% get()
  # CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat) #Necessary
  
  #Preprocessing the expression data for cell-cell communication analysis
  cellchat <- identifyOverExpressedGenes(cellchat) #p value = 0.05 and L2FC = 0
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  #Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat, seed.use = 1) 
  cellchat <- filterCommunication(cellchat, min.cells = 10) #filter fewer groups
  
  #Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05) #pvalue 0.05
  
  cellchat <- netAnalysis_computeCentrality(cellchat) #neccesary for the heatmaps and 2D 
  
  #Calculate the aggregated cell-cell communication network by counting the number of links 
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

Healthy <- Cellchat_pipeline(sn_MASH_filtered, "f_score_12", "Healthy", "cluster_annotation")
f12  <- Cellchat_pipeline(sn_MASH_filtered, "f_score_12", "f12", "cluster_annotation")
f34  <- Cellchat_pipeline(sn_MASH_filtered, "f_score_12", "f34", "cluster_annotation")

cellchat_whole_list <- list(Healthy = Healthy, f12 = f12, f34 = f34)
cellchat_sn_MASH <- mergeCellChat(cellchat_whole_list, add.names = names(cellchat_whole_list))

#obtain L-R matrix
L_R_matrix <- c()
for (state in c("Healthy", "f12", "f34")) {
  object <- cellchat_whole_list[[state]]
  pathways <- object@netP$pathways
  for (pathway in pathways) {
    matrix <- subsetCommunication(object, slot.name = "net", signaling = pathway) %>% arrange(desc(prob))
    matrix$state <- state
    L_R_matrix <- rbind(L_R_matrix, matrix)
  }
}

#compositional clusters
library(ISCHIA)
df_fibrosis <- read.csv(path_abundance_matrix, row.names = 1)
wss_fibrosis <- Composition.cluster.k(df_fibrosis, 10)

ST_CC <- Composition.cluster(ST_MASH, df, 4)

#calculate co-ocurrence
CC_co_occurrence <- list()
for (i in 1:length(ST_CC$CompositionCluster_CC %>% unique())) {
  CC_co_occurrence[[i]] <- spatial.celltype.cooccurence(spatial.object=ST_CC,
                                                        deconv.prob.mat=df_fibrosis,
                                                        COI=glue("CC{i}"), prob.th= 0.2, Condition=unique(ST_CC$orig.ident))
}

#obtain the matrix of the co-ocurrence per CC
matrix_co_ocurr <- function(x){
  dim <- x$species
  comat_pos <- comat_neg <- matrix(nrow=dim,ncol=dim)
  co_tab <- x$result
  for (i in 1:nrow(co_tab)){
    comat_pos[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_gt"]
    comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_gt"]
  
    row.names(comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]])
  }
  
  for (i in 1:nrow(co_tab)){
    comat_neg[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_lt"]
    comat_neg[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_lt"]
  }
  comat <- ifelse(comat_pos>=0.05,0,1) + ifelse(comat_neg>=0.05,0,-1)
  colnames(comat) <- 1:dim
  row.names(comat) <- 1:dim
  
  if ("spp_key" %in% names(x)){

  sp1_name <- merge(x=data.frame(order=1:length(colnames(comat)),sp1=colnames(comat)),y=x$spp_key,by.x="sp1",by.y="num",all.x=T)
  sp2_name <- merge(x=data.frame(order=1:length(row.names(comat)),sp2=row.names(comat)),y=x$spp_key,by.x="sp2",by.y="num",all.x=T)

  colnames(comat) <- sp1_name[with(sp1_name,order(order)),"spp"]
  row.names(comat) <- sp2_name[with(sp2_name,order(order)),"spp"]
  }
  comat[is.na(comat)] <- 0
  ind <- apply(comat, 1, function(x) all(x==0))
  comat <- comat[names(sort(ind)),]
  ind <- apply(comat, 2, function(x) all(x==0))
  comat <- comat[,names(sort(ind))]
  
  data.m = reshape2::melt(comat)
  colnames(data.m) <- c("X1","X2","value")
  data.m$X1 <- as.character(data.m$X1)
  data.m$X2 <- as.character(data.m$X2)
  
  dfids <- subset(data.m, X1 == X2)
  X1 <- data.m$X1
  X2 <- data.m$X2
  df.lower = subset(data.m[lower.tri(comat),],X1 != X2)
  return(df.lower)
}


co_oc_1 <- CC_co_occurrence[[1]] %>% 
  matrix_co_ocurr() %>% 
  filter(value ==1) %>% 
  mutate(CC = "CC1")
  
#filter co-ocurrence
df_filtered <- inner_join(L_R_matrix, co_oc_1, by = c("source" = "X1", "target" = "X2")) %>%
  bind_rows(inner_join(L_R_matrix, co_oc_1, by = c("source" = "X2", "target" = "X1"))) 


#DE

L_R_DE <- function(cellchat_object, cond1, cond2, L_up=0.1, R_up=0.1, L_down=-0.1, R_down=-0.1){
  # return L-R up and down in the second conditions
  
  # perform differential expression analysis
  cellchat_object <- identifyOverExpressedGenes(cellchat_object, group.dataset = "datasets", pos.dataset = cond2, features.name = cond2, only.pos = FALSE, thresh.fc = 0.1, thresh.p = 0.05)
  
  #> Use the joint cell labels from the merged CellChat object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(cellchat_object, features.name = cond2, thresh = 0.05)
  
  # extract the ligand-receptor pairs with upregulated ligands & upregulated receptors in second condition(disease)
  net.up <- subsetCommunication(cellchat_object, net = net, datasets = cond2, ligand.logFC = L_up, receptor.logFC = R_up, ligand.pct.1 = 0.25, ligand.pct.2 = 0.25, receptor.pct.1 = 0.25, receptor.pct.2 = 0.25)
  
  
  # extract the ligand-receptor pairs downregulated in second condition
  net.down <- subsetCommunication(cellchat_object, net = net, datasets = cond1 ,ligand.logFC = L_down, receptor.logFC = R_down, ligand.pct.1 = 0.25, ligand.pct.2 = 0.25, receptor.pct.1 = 0.25, receptor.pct.2 = 0.25)
  
  return(list("up" = net.up, "down"=  net.down))
}

H_f12 <- L_R_DE(cellchat_whole_list, "Healthy", "f12")
H_f34 <- L_R_DE(cellchat_whole_list, "Healthy", "f34")
f12_f34<- L_R_DE(cellchat_whole_list, "f12", "f34")












