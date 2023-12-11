library(ISCHIA)
library(extrafont)
library(wesanderson)


sn_NASH <- readRDS("")
sn_CKD <- readRDS("")


CKD <- sn_CKD@meta.data %>%  mutate(dataset = "Kidney") %>% as.tibble() %>% 
  dplyr::select(nCount_RNA, nFeature_RNA, percent.mt, specimen_id, dataset)

CKD$specimen_id <- CKD$specimen_id %>% as.factor()

NASH <- sn_NASH@meta.data %>% mutate(dataset= "Liver") %>% as.tibble() %>% 
  dplyr:: select(nCounts, nGenes, pct_mito, library_id, dataset)  %>% 
  set_names(c("nCount_RNA",
              "nFeature_RNA",
              "percent.mt",
              "specimen_id",
              "dataset"))

whole <- rbind(CKD, NASH)
my_colors <- c("#F9A825", "#C62828")


p1<- ggplot(whole, aes(x = specimen_id, color = dataset, y = nCount_RNA, fill=dataset)) + 
  geom_boxplot(alpha = 0.2) + ggtitle("Count") + 
  scale_x_discrete( labels = as.character(1:47)) +
  theme(
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_blank(),
  ) +
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors)



p2<- ggplot(whole, aes(x = specimen_id, color = dataset, y = nFeature_RNA, fill=dataset)) + 
  geom_boxplot(alpha = 0.2) + ggtitle("Number features") + 
  scale_x_discrete( labels = as.character(1:47)) +
  theme(
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_blank(),
  ) +
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors)


p3<- ggplot(whole, aes(x = specimen_id, color = dataset, y = percent.mt, fill=dataset)) + 
  geom_boxplot(alpha = 0.2) + ggtitle("Mitochondrial percent") + 
  scale_x_discrete( labels = as.character(1:47)) +
  theme(
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_blank(),
  ) +
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors)



grid.arrange(p1,p2,p3, ncol = 1, nrow = 3, top = textGrob("scRNAseq data",gp=gpar(fontsize=14)))



###QC ST
list_NASH <- readRDS("")
list_CKD <- readRDS("")

QC_NASH <- data.frame()

for (i in 1:length(list_NASH)){
  QC_NASH <- rbind(QC_NASH, data.frame(sample_id = list_NASH[[i]]$orig.ident, 
                                       nCount_Spatial = list_NASH[[i]]$nCount_Spatial,
                                       nFeature_Spatial =  list_NASH[[i]]$nFeature_Spatial,
                                       percent_mito = list_NASH[[i]]$percent_mito,
                                       percent_hb = list_NASH[[i]]$percent_hb,
                                       fibro_02_34 = list_NASH[[i]]$fibro_02_34,
                                       disease_status = list_NASH[[i]]$disease_status,
                                       spots = (Cells(list_NASH[[i]]) %>% length())
  ))
}

QC_CKD <- data.frame()

for (i in 1:length(list_CKD)){
  QC_CKD <- rbind(QC_CKD, data.frame(sample_id = list_CKD[[i]]$orig.ident, 
                                     nCount_Spatial = list_CKD[[i]]$nCount_Spatial,
                                     nFeature_Spatial =  list_CKD[[i]]$nFeature_Spatial,
                                     percent_mito = list_CKD[[i]]$percent_mito,
                                     percent_hb = list_CKD[[i]]$percent_hb,
                                     disease_status = list_CKD[[i]]$disease_status,
                                     spots = (Cells(list_CKD[[i]]) %>% length())
  ))
}

QC_CKD$Dataset <- "Kidney"
QC_NASH$Dataset <- "Liver"

QC_NASH$fibro_02_34<-NULL

Whole_QC <- rbind(QC_CKD, QC_NASH)


my_colors <- c("#F9A825", "#C62828")

p1<- ggplot(Whole_QC, aes(x = sample_id, color = Dataset, y = nCount_Spatial, fill=Dataset)) + 
  geom_boxplot(alpha = 0.2) + ggtitle("Spatial count") + 
  scale_x_discrete( labels = as.character(1:39)) +
  theme(
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_blank(),
  ) +
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors)

p2 <- ggplot(Whole_QC, aes(x = sample_id, color=Dataset, y = nFeature_Spatial, fill=Dataset)) + 
  geom_boxplot(alpha = 0.2) + ggtitle("Number features ")  +
  scale_x_discrete( labels = as.character(1:39)) +
  theme(
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors) 

p3 <- ggplot(Whole_QC, aes(x = sample_id, color=Dataset, y = spots, fill=Dataset)) +
  geom_col(width=0.5, alpha = 0.7) + ggtitle("Number spots ")  +
  scale_x_discrete( labels = as.character(1:39)) +
  theme(
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors) 



grid.arrange(p1,p2,p3, ncol = 1, nrow = 3, top = textGrob("Spatial Transcriptomics Data",gp=gpar(fontsize=14)))




#### sn kidney cell type
sn_kidney <- readRDS("")


n_cells <- FetchData(sn_kidney, 
                     vars = c("subclass.l1", "Tissue.Type")) %>%
  dplyr::count(subclass.l1, Tissue.Type) 


n_cells$Tissue.Type <- factor(n_cells$Tissue.Type, levels = c("Healthy", "AKI", "CKD"))

# Define the colors for each bar
bar_colors <- wes_palette("GrandBudapest1", n = 3)

ggplot(data = n_cells, aes(x = subclass.l1, y = n, fill = Tissue.Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), hjust= -0.1, size=3, angle=90) +
  
  # Customize the plot theme
  theme_classic() +
  xlab("Cell Type") +
  ylab("Number of Cells") +
  
  # Rotate x-axis labels for better readability
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  # Set the order of the bars and the bar colors
  scale_x_discrete(limits = rev(rownames(matrix))) +
  scale_fill_manual(values = bar_colors) +
  
  # Remove the legend title and change the legend position and appearance
  labs(fill = NULL) +
  theme(legend.position = "top", legend.background = element_rect(fill = "white", size = 0.5)) +
  guides(fill = guide_legend(title = "Tissue Type")) +
  ylim(0,22000)





#sn liver barplot by celltype


sn_liver <- readRDS("")



sn_liver$f_score_12 <- sn_liver$f_score_02 %>% as.character()

sn_liver$f_score_12 <- ifelse(sn_liver$f_score_12 == "f02", "F12",  sn_liver$f_score_12)
sn_liver$f_score_12 <- ifelse(sn_liver$f_score_12 == "f34", "F34",  sn_liver$f_score_12)

n_cells <- FetchData(sn_liver, 
                         vars = c("f_score_12", "cluster_annotation")) %>%
  dplyr::count(f_score_12, cluster_annotation)

n_cells$f_score_12 <- factor(n_cells$f_score_12, levels = c("Healthy", "F12", "F34"))

bar_colors <- wes_palette("GrandBudapest1", n = 3)

ggplot(data = n_cells, aes(x = cluster_annotation, y = n, fill = f_score_12)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), hjust= -0.1, size=3, angle=90) +
  
  # Customize the plot theme
  theme_classic() +
  xlab("Cell Type") +
  ylab("Number of Cells") +
  
  # Rotate x-axis labels for better readability
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  
  # Set the order of the bars and the bar colors
  scale_x_discrete(limits = rev(rownames(matrix))) +
  scale_fill_manual(values = bar_colors) +
  
  # Remove the legend title and change the legend position and appearance
  labs(fill = NULL) +
  theme(legend.position = "top", legend.background = element_rect(fill = "white", size = 0.5)) +
  guides(fill = guide_legend(title = "Tissue Type")) +
  ylim(0,80000)




#harmony
library(RColorBrewer)
ST_harmonized_kidney <- readRDS("")

my_palette <- c("#FFFF00", "#FF6600")

p1 <- DimPlot(object = ST_harmonized_kidney, reduction = "harmony", pt.size = .1, group.by = "disease_status") +
  ggtitle("Kidney") +
  scale_color_manual(values = my_palette) +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "top") +
  xlab("Component 1") +
  ylab("Component 2") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(arrow = arrow(angle = 15, length = unit(.2, "cm"), type = "closed")))



ST_harmonized_liver <- readRDS("")
my_palette <- "#C51B8A"
my_palette <- c(my_palette, (brewer.pal(3, "Reds")[2:3]))

ST_harmonized_liver$fibro_12_34 <- ST_harmonized_liver$fibro_02_34 %>% as.character()

ST_harmonized_liver$fibro_12_34 <- ifelse(ST_harmonized_liver$fibro_12_34 == "F02", "F12",  ST_harmonized_liver$fibro_12_34)
# ST_harmonized_liver$fibro_12_34 <- ifelse(ST_harmonized_liver$fibro_12_34 == "f34", "F34",  ST_harmonized_liver$fibro_12_34)
ST_harmonized_liver$fibro_12_34 <- factor(ST_harmonized_liver$fibro_12_34, levels = c("Healthy", "F12", "F34"))

p2 <- DimPlot(object = ST_harmonized_liver , reduction = "harmony", pt.size = .1, group.by = "fibro_12_34") +
  ggtitle("Liver") +
  scale_color_manual(values = my_palette) +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "top") +
  xlab("Component 1") +
  ylab("Component 2") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(arrow = arrow(angle = 15, length = unit(.2, "cm"), type = "closed"))) 

library(ggpubr)
grid.arrange(p1, p2, nrow=1, top= text_grob("Harmony Integration Spatial Data", size = 15, face = "bold"))




#MASH CC
ggplot(df_times_cluster, aes(x = CompositionCluster_CC, y = Freq, fill = disease_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Composition Cluster", y = "Spots", fill = "Disease state") +
  # ggtitle("Number of spots per disease state per CC ") +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(face = "bold", size = 16),
        legend.position = "top") 


##Elbow
library(ISCHIA)
df_fibrosis <- read.csv("cell_abundance_q5.csv", row.names = 1)
df_kidney <- read.csv("cell_abundance_q5.csv", row.names = 1)

fix(Composition.cluster.k) #return wss

wss_fibrosis <- Composition.cluster.k(df_fibrosis, 15)
wss_kidney <- Composition.cluster.k(df_kidney, 15) 

#normalize values between 0-1 to compare between datasets

normalize_minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
} 


ggplot() +
  geom_line( aes(1:15, normalize_minmax(wss_fibrosis), color = "1"), size = 1) +
  geom_point(aes(1:15, normalize_minmax(wss_fibrosis), color = "1"), size = 2) +
  geom_point(aes(4, normalize_minmax(wss_fibrosis)[4]), color = "black", size = 2) + 
  geom_line( aes(1:15, normalize_minmax(wss_kidney),  color = "2"), size = 1) +
  geom_point(aes(1:15, normalize_minmax(wss_kidney), color = "2"), size = 2) +
  geom_point(aes(6, normalize_minmax(wss_kidney)[6]),   color = "black", size = 2) +
  labs(title = "Elbow Plots",
       x = "Number of Clusters",
       y = "WSS") +
  scale_color_manual(values = c("1" = "#C62828", "2" = "#F9A825"), labels = c("Liver", "Kidney"), name = "Dataset") +
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 



#correlation
markers <-  c("LUM", "IGFBP7", "COL1A1")
ligand <- "PDGFC"
receptor <- "PDGFRA"
genes <- c(ligand,receptor,markers)


list_ST_NASH <- readRDS()
Healthy_list <- Filter(function(x) "Healthy" %in% levels(x$disease_status), list_ST_NASH)
NASH_list <- Filter(function(x) "NASH" %in% levels(x$disease_status), list_ST_NASH)

correlation_spatial <- function(genes, ST_list, assay = "SCT"){
  co_expression_df <- data.frame(matrix(nrow = length(genes), ncol = length(genes)), row.names=genes)
  colnames(co_expression_df) <- genes
  
  combinations <- combn(genes,2)
  
  for (i in 1:ncol(combinations)) {
    pair_genes <- combinations[,i]
    cor_v <- numeric(length(ST_list))
    
    for(j in 1:length(ST_list)){
      data_j <- ST_list[[j]]@assays[[assay]]@data[pair_genes, ] 
      cor_v[j] <- cor(data_j[1, ], data_j[2, ], method = c("pearson"))
    }
    co_expression_df[pair_genes[1], pair_genes[2]] <- mean(cor_v)
    co_expression_df[pair_genes[2], pair_genes[1]] <- mean(cor_v)
  }
  diag(co_expression_df) <- 1.0
  
  return(as.matrix(co_expression_df))
}


co_expression <- correlation_spatial(genes, list_ST_NASH)
co_expression_Healthy <- correlation_spatial(genes, Healthy_list)
co_expression_NASH <- correlation_spatial(genes, NASH_list)


corrplot(co_expression, method="number", type = "lower", col = rev(COL2('RdBu', 100)), tl.col= "black", title = "All samples", mar=c(0,0,1,0))
corrplot(co_expression_Healthy, method="number", type = "lower", col = rev(COL2('RdBu', 100)), tl.col= "black", title = "Healthy samples", mar=c(0,0,1,0))
corrplot(co_expression_NASH, method="number", type = "lower", col = rev(COL2('RdBu', 100)), tl.col= "black", title = "MASH samples", mar=c(0,0,1,0))