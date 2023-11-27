library(ISCHIA)
library(extrafont)
library(wesanderson)

sn_liver <- readRDS("")


#sn liver barplot by celltype
n_cells <- FetchData(sn_liver, 
                         vars = c("f_score_12", "cluster_annotation")) %>%
  dplyr::count(f_score_02, cluster_annotation)

n_cells$f_score_12 <- factor(n_cells$f_score_12, levels = c("Healthy", "f12", "f34"))

ggplot(data = n_cells_old, aes(x = cluster_annotation, y = n, fill = f_score_12)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), hjust= -0.1, size=3, angle=90) +
  
  # Customize the plot theme
  theme_classic() +
  xlab("Cell Type") +
  ylab("Number of Cells") +
  ggtitle("Liver single nuclei dataset") +
  
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
ST_harmonized_kidney <- readRDS()

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



ST_harmonized_liver <- readRDS()
my_palette <- "#C51B8A"
my_palette <- c(my_palette, (brewer.pal(3, "Reds")[2:3]))

ST_harmonized_liver$F_score_12 <- factor(ST_harmonized_liver$F_score_12, levels = c("Healthy", "F12", "F34"))

p2 <- DimPlot(object = ST_harmonized_liver , reduction = "harmony", pt.size = .1, group.by = "F_score_12") +
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




#correlation
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

corrplot(co_expression, method="number", type = "lower", col = rev(COL2('RdBu', 100)), tl.col= "black")
corrplot(co_expression_Healthy, method="number", type = "lower", col = rev(COL2('RdBu', 100)), tl.col= "black")
corrplot(co_expression_NASH, method="number", type = "lower", col = rev(COL2('RdBu', 100)), tl.col= "black")