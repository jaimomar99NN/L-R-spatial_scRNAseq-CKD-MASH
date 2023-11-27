library(ISCHIA)
library(extrafont)
library(wesanderson)


ST_liver_list <- readRDS("")
ST_liver_integrated <- readRDS("")
abundance_matrix <- read.csv("", row.names = 1)



#KIDNEY
#composition clusters
# ST_liver_integrated <- Composition.cluster(ST_liver_integrated, abundance_matrix, 4)  #be careful the order of the cluster changes every time is run# ST_liver_integrated <- Composition.cluster(ST_liver_integrated, df, 6)  #be careful the order of the cluster changes every time is run
CC_boxplots <- list()
for (i in 1:length(ST_liver_integrated$CompositionCluster_CC %>% unique())) {
  CC_boxplots[[i]] <- Composition_cluster_enrichedCelltypes(ST_liver_integrated, glue("CC{i}"), as.matrix(abundance_matrix)) +
    ggtitle(glue("Composition Cluster {i}")) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}
grid.arrange(grobs = CC_boxplots)


#contribution per condition
df_times_cluster <- data.frame(table(ST_liver_integrated@meta.data %>% dplyr:: select(c("CompositionCluster_CC", "fibro_02_34"))))
df_times_cluster$fibro_02_34 <- factor(df_times_cluster$fibro_02_34, levels = c("Healthy", "F02", "F34"))

colors <- wes_palette("GrandBudapest1", n = 3)

ggplot(df_times_cluster, aes(x = CompositionCluster_CC, y = Freq, fill = fibro_02_34)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Composition Cluster", y = "Spots", fill = "Disease state") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(face = "bold", size = 16),
        legend.position = "top") 

#co-ocurrence plot
CC_co_occurrence <- list()
CC_co_occurrence_plot <- list()
for (i in 1:length(unique(ST_liver_integrated$CompositionCluster_CC))) {
  CC_co_occurrence[[i]] <- spatial.celltype.cooccurence(spatial.object=ST_liver_integrated,
                                                        deconv.prob.mat=abundance_matrix,
                                                        COI=glue("CC{i}"), prob.th= 0.2, Condition=unique(ST_liver_integrated$orig.ident))
  
  CC_co_occurrence_plot[[i]] <- plot.celltype.cooccurence(CC_co_occurrence[[i]]) + 
    ggtitle(glue(" Co-occurrence CC{i}")) + 
    theme(legend.position="bottom",
          panel.border = element_blank())
}

grid.arrange(CC_co_occurrence_plot[[2]])




#PROS1-AXL 
for (i in 1:length(ST_liver_list)) {
  ST_correlated_neighboring_quantity(ST_liver_list[[i]], c("PROS1","AXL"))
}

# PDGFC-PDGFRA quantitative & qualitative & fibrotic-inflammation markers
for (i in 1:length(ST_liver_list)) {
  ST_correlated_coexpression_quantity(ST_liver_list[[i]], c("PDGFC","PDGFRA"))
  ST_correlated_coexpression_quality(ST_liver_list[[i]], c("PDGFC","PDGFRA"))
  ST_correlated_coexpression_quality(ST_liver_list[[i]], c("LUM", "IGFBP7", "AKAP9","USP33"))
}

