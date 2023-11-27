ST_correlated_neighboring_quantity <- function(ST_object, genes, variable = "disease_status"){
  library(glue)
  condition <- ifelse(!is.null(variable), pull(unique(ST_object[[variable]])), "")
  
  for (gene in genes) {
    if (!gene %in% rownames(ST_object)) {
      print(glue("{gene} is not sequenced in {unique(ST_object$orig.ident)} – {condition}"))
      return()
    }
  }
  
  exp_mat = GetAssayData(ST_object)[genes,] %>% data.frame() %>% t()
  
  for (i in 1:length(genes)) {
    if (all(exp_mat[,i] %>% unique() == 0)) {
      print(glue("{genes[i]} is not expressed in {unique(ST_object$orig.ident)} – {condition}"))
      return()
    }
  }
  
  
  # Normalize gene expression between 0 and 1
  normalizing_data = function(values){return((values - min(values)) / (max(values) - min(values)))}
  exp_mat_norm = apply(exp_mat,2,normalizing_data) %>% data.frame()
  
  #taking min  
  exp_mat_norm$min_both = as.vector(apply(exp_mat_norm,1,min))
  
  #str column, which one we pick
  exp_mat_norm$selected <- ifelse(exp_mat_norm$min_both > 0, 'Both',
                                  ifelse(exp_mat_norm[[genes[1]]] > 0 , genes[1], ifelse(exp_mat_norm[[genes[2]]] > 0, genes[2],
                                                                                         'None')))
  
  max_both <- max(exp_mat_norm$min_both) * 100
  
  if (is.nan(max_both)){max_both = 0}
  
  # Create color gradients based on expression values
  gene1_colors <- colorRampPalette(c("#C9FFC9", "#00A000"))(100 + 1) #green
  gene2_colors <- colorRampPalette(c("#C6E0FF", "#015CC8"))(100 + 1) #blue
  both_colors <- colorRampPalette(c("#F97979", "#DE0418"))(max_both + 1) #red
  
  index_vector <- function(val, color_pal) {
    return(color_pal[val * 100 + 1])
  }
  
  exp_mat_norm$color_L <- sapply(exp_mat_norm[[genes[1]]], index_vector, color_pal = gene1_colors)
  exp_mat_norm$color_R <- sapply(exp_mat_norm[[ genes[2]]], index_vector, color_pal = gene2_colors)
  exp_mat_norm$color_both <- sapply(exp_mat_norm[["min_both"]], index_vector, color_pal = both_colors)
  
  #final color, which one we pick
  exp_mat_norm$final_color <- ifelse(exp_mat_norm$selected == "Both", exp_mat_norm$color_both,
                                     ifelse(exp_mat_norm$selected == genes[1] , exp_mat_norm$color_L, ifelse(exp_mat_norm$selected ==  genes[2], exp_mat_norm$color_R,
                                                                                                             "#DAD9D9")))
  # ST_object map to metadata object
  ST_object$colors <- exp_mat_norm$final_color 
  
  #important to be detected by seurat
  colors <- exp_mat_norm$final_color
  names(colors) <- exp_mat_norm$final_color 
  
  
  
  #DUMB data to take legend
  #No coexpression at all
  if (max_both == 0) {
    p <- ggplot() +
      #gene2
      geom_tile(aes(x = 1:101, y = 101:201, fill = rnorm(101))) +
      scale_fill_gradient2( genes[2], low = gene2_colors[1], high= gene2_colors[length(gene2_colors)], limits = c(0, 1), breaks = c(0, 1)) +
      new_scale("fill") +
      
      #gene1
      geom_tile(aes(x = 1:101, y = 1:101, fill = rnorm(101))) +
      scale_fill_gradient2( genes[1], low = gene1_colors[1], high= gene1_colors[length(gene1_colors)], limits = c(0, 1), breaks = c(0, 1)) +
      new_scale("fill") +
      theme(legend.direction = "horizontal")+
      labs(fill = "Value") + 
      guides(fill = guide_colorbar(title.position = "left"))
  } else{
    p <- ggplot() +
      #gene2
      geom_tile(aes(x = 1:101, y = 101:201, fill = rnorm(101))) +
      scale_fill_gradient2( genes[2], low = gene2_colors[1], high= gene2_colors[length(gene2_colors)], limits = c(0, 1), breaks = c(0, 1)) +
      new_scale("fill") +
      
      #gene1
      geom_tile(aes(x = 1:101, y = 1:101, fill = rnorm(101))) +
      scale_fill_gradient2( genes[1], low = gene1_colors[1], high= gene1_colors[length(gene1_colors)], limits = c(0, 1), breaks = c(0, 1)) +
      new_scale("fill") +
      
      #both
      geom_tile(aes(x = 1:max_both, y = 1:max_both, fill = rnorm(max_both))) +
      scale_fill_gradient2("BOTH", low = both_colors[1], high= both_colors[length(both_colors)], limits = c(0, round(max_both/100, 2)), breaks = c(0, round(max_both/100, 2))) +
      theme(legend.direction = "horizontal")+
      labs(fill = "Value") + 
      guides(fill = guide_colorbar(title.position = "left"))
  }
  
  p1 <- SpatialDimPlot(ST_object, group.by = "colors", crop = FALSE, pt.size.factor = 1.2) + 
    scale_fill_manual(values = colors) + 
    NoLegend() +
    theme(plot.margin = margin(t = 0, r=-1.5, b=0, l=0, "cm"))
  
  grid.arrange(p1, get_legend(p), nrow=1, widths = c(2,1))
  
}

ST_correlated_neighboring_quality <- function(ST_object, genes, cut_off = 0.1, variable = "disease_status"){
  
  condition <- ifelse(!is.null(variable), ST_object[[variable]] %>% unique() %>% pull(), "")
  
  for (gene in genes) {
    if (!gene %in% rownames(ST_object)) {
      print(glue("{gene} is not sequenced in {ST_object$orig.ident %>% unique()} – {condition}"))
      return()
    }
  }
  
  exp_mat = GetAssayData(ST_object)[genes,] %>% data.frame() %>% t()
  
  for (i in c(1,2)) {
    if (all(exp_mat[,i] %>% unique() == 0)) {
      print(glue("one of the genes is not expressed in {ST_object$orig.ident %>% unique()} – {condition}"))
      return()
    }
  }
  
  # Normalize gene expression between 0 and 1
  normalizing_data = function(values){return((values - min(values)) / (max(values) - min(values)))}
  exp_mat_norm = apply(exp_mat,2,normalizing_data) %>% data.frame()
  
  #taking min  
  exp_mat_norm$min_both = as.vector(apply(exp_mat_norm,1,min))
  
  #str column, which one we pick
  exp_mat_norm$selected <- ifelse(exp_mat_norm$min_both > cut_off, 'Both',
                                  ifelse(exp_mat_norm[[genes[1]]] > cut_off , genes[1], ifelse(exp_mat_norm[[genes[2]]] > cut_off, genes[2],
                                                                                         'None')))
  
  #assign color
  exp_mat_norm$final_color <- ifelse(exp_mat_norm$selected == "Both", "#DE0418",
                                     ifelse(exp_mat_norm$selected == genes[1] , "#00A000", ifelse(exp_mat_norm$selected == genes[2], "#015CC8",
                                                                                               "#DAD9D9")))
  
  # ST_object map to metadata object
  ST_object$colors <- exp_mat_norm$final_color 
  
  #important to be detected by seurat
  colors <- exp_mat_norm$final_color
  names(colors) <- exp_mat_norm$final_color 
  
  
  p1 <- SpatialDimPlot(ST_object, group.by = "colors", crop = FALSE, pt.size.factor = 1.2) + 
    scale_fill_manual(values = colors) + 
    NoLegend() +
    theme(plot.margin = margin(t = 0, r=-1.5, b=0, l=0, "cm"))
  
  
  #DUMB data for the legend
  p <- ggplot() +
    geom_point(aes(c(1:3), c(1:3), color = c("#DE0418", "#00A000", "#015CC8"))) +
    scale_color_manual(values = c("#DE0418", "#00A000", "#015CC8"),
                       labels = c("Both", genes[1], genes[2]),
                       name = glue("Threshold: {cut_off}")) 
  
  grid.arrange(p1, get_legend(p), nrow=1, widths = c(4,1), top = textGrob(glue("{ST_object$orig.ident} – {condition}")))
}

ST_correlated_coexpression_quantity <- function(ST_object, genes, cut_off = 0.1, variable = "disease_status"){
  
  condition <- ifelse(!is.null(variable), ST_object[[variable]] %>% unique() %>% pull(), "")
  
  for (gene in genes) {
    if (!gene %in% rownames(ST_object)) {
      print(glue("{gene} is not sequenced in {ST_object$orig.ident %>% unique()} – {condition}"))
      return()
    }
  }
  
  exp_mat = GetAssayData(ST_object)[genes,] %>% data.frame() %>% t()
  
  for (i in c(1,2)) {
    if (all(exp_mat[,i] %>% unique() == 0)) {
      print(glue("one of the genes is not expressed in {ST_object$orig.ident %>% unique()} – {condition}"))
      return()
    }
  }
  
  # Normalize gene expression between 0 and 1
  normalizing_data = function(values){return((values - min(values)) / (max(values) - min(values)))}
  exp_mat_norm = apply(exp_mat,2,normalizing_data) %>% data.frame()
  
  #taking min  
  exp_mat_norm$min_both = as.vector(apply(exp_mat_norm,1,min))
  
  #str column, which one we pick
  exp_mat_norm$selected <- ifelse(exp_mat_norm$min_both > cut_off, 'Both',
                                  ifelse(exp_mat_norm[[gene1]] > cut_off , gene1, ifelse(exp_mat_norm[[gene2]] > cut_off, gene2,
                                                                                         'None')))
  
  #assign color
  exp_mat_norm$final_color <- ifelse(exp_mat_norm$selected == "Both", "#DE0418",
                                     ifelse(exp_mat_norm$selected == gene1 , "#00A000", ifelse(exp_mat_norm$selected == gene2, "#015CC8",
                                                                                               "#DAD9D9")))
  
  # ST_object map to metadata object
  ST_object$colors <- exp_mat_norm$final_color 
  
  #important to be detected by seurat
  colors <- exp_mat_norm$final_color
  names(colors) <- exp_mat_norm$final_color 
  
  
  p1 <- SpatialDimPlot(ST_object, group.by = "colors", crop = FALSE, pt.size.factor = 1.2) + 
    scale_fill_manual(values = colors) + 
    NoLegend() +
    theme(plot.margin = margin(t = 0, r=-1.5, b=0, l=0, "cm"))
  
  
  #DUMB data for the legend
  p <- ggplot() +
    geom_point(aes(c(1:3), c(1:3), color = c("#DE0418", "#00A000", "#015CC8"))) +
    scale_color_manual(values = c("#DE0418", "#00A000", "#015CC8"),
                       labels = c("Both", gene1, gene2),
                       name = glue("Threshold: {cut_off}")) 
  
  grid.arrange(p1, get_legend(p), nrow=1, widths = c(4,1), top = textGrob(glue("{ST_object$orig.ident} – {condition}")))
}

ST_correlated_coexpression_quality <- function(ST_object, genes, slot="data", cut_off = 0.1, variable = "disease_status"){
  condition <- ifelse(!is.null(variable), ST_object[[variable]] %>% 
                        unique() %>%
                        pull() %>%
                        as.character(),
                      "")
  
  for (gene in genes) {
    if (!gene %in% rownames(ST_object)) {
      print(glue("{gene} is not sequenced in {ST_object$orig.ident %>% unique()} – {condition}"))
      return()
    }
  } 
  exp_mat = GetAssayData(ST_object, slot = slot)[genes,] %>% data.frame() %>% t()
  
  for (i in 1:length(genes)) {
    if (all(exp_mat[,i] %>% unique() == 0)) {
      print(glue("{genes[i]} is not expressed in {ST_object$orig.ident %>% unique()} – {condition}"))
      return()
    }
  }
  
  if (slot == "data") {
    # Normalize gene expression between 0 and 1
    normalizing_data = function(values){return((values - min(values)) / (max(values) - min(values)))}
    exp_mat_norm = apply(exp_mat,2,normalizing_data) %>% data.frame()
  } else if (slot == "counts") { 
    exp_mat_norm = exp_mat
  } else { print("review function")}
  
  
  #taking min  
  exp_mat_norm$min_all = as.vector(apply(exp_mat_norm,1,min))
  
  #str column, which one we pick
  exp_mat_norm$selected <- ifelse(exp_mat_norm$min_all > cut_off, 'ALL', 'None')
  
  max_all <- max(exp_mat_norm$min_all) * 100
  
  if (is.nan(max_all)){max_all = 0}
  
  #final color, which one we pick
  exp_mat_norm$final_color <- ifelse(exp_mat_norm$selected == "ALL", "#DE0418", "#DAD9D9")
  
  # ST_object map to metadata object
  ST_object$colors <- exp_mat_norm$final_color 
  
  #important to be detected by seurat
  colors <- exp_mat_norm$final_color
  names(colors) <- exp_mat_norm$final_color 
  
  #DUMB data to take legend
  #No coexpression at all
  if (max_all == 0) {
    print(glue("No coexpression of {paste(genes, collapse = ',')} in {ST_object$orig.ident %>% unique()} – {condition}"))
    return()
  } 
  p <- SpatialDimPlot(ST_object, group.by = "colors", crop = FALSE, pt.size.factor = 1.2) + 
    scale_fill_manual(values = colors) + 
    NoLegend() +
    ggtitle(glue("{ST_object$orig.ident} – {condition}")) +
    labs(subtitle = glue(" {paste(genes, collapse=',')}. Th:{cut_off}")) +
    theme(plot.subtitle=element_text(size=8,color="black"),
          plot.title= element_text(size = 12, face="bold"))
  
  print(p)
  
}