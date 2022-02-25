# Plot average gene expression per cluster (color) and visualize the fraction of cells passing 
# a given expression threshold (dot size) using. Works on sc objects from the RaceID package.


fraction_dotplot <- function(sc_object, 
                          genes, 
                          logsc = T, 
                          scaled = F, 
                          cex = 11, 
                          colorscale = c('black','darkblue', 'royalblue','turquoise', 'lightgreen', 'gold','orange', 'red2', 'red4'),
                          title_size = 15,
                          min_transcript = 1){
  
  # UMAP coordinates and cluster partition
  # data.frame with columns cellid (character), cluster (factor) , X (numeric) , Y (numeric) , is_medoid (bool)
  plotmap_df <- sc_object@cpart %>% 
    as.data.frame() %>%
    rownames_to_column() %>%
    # Bind columns since cells in sc_object@cpart and @umap slot should be ordered the same way
    bind_cols(sc_object@umap) %>%
    magrittr::set_colnames(c("cellid", "cluster", "X", "Y")) %>%
    mutate("is_medoid"= cellid %in% sc_object@medoids, 
           cluster = factor(cluster, levels = 1:length(unique(cluster)) ) ) 
  
  # Gene expression values for selected genes
  # columns: cellid, cluster X, Y, is_medoid, variable(genesymbol; character), value (numeric), expr_thr (bool)
  genes_df <- cbind(plotmap_df, 
                    (sc_object@ndata[genes, rownames(plotmap_df), drop = F] * sc_object@counts) %>% 
                      as.matrix() %>% t() + 0.1) %>% 
    data.table::data.table() %>%
    data.table::melt(., 
                     id.vars = c('cellid', 'cluster', 'X', 'Y', 'is_medoid'), 
                     measure.vars = genes) %>% 
    mutate(expr_thr = value >= min_transcript)
  
  # Summary statistic for selected genes per cluster
  ggplot_df <- genes_df %>% 
    group_by(cluster, variable) %>% 
    summarise(count_expr_thr = sum(expr_thr),
              count          = n(),
              fraction       = count_expr_thr/count,
              mean_expr      = mean(value),
              median_expr    = median(value),
              sd_expr        = sd(value),
              norm_sd        = sd_expr/mean_expr
              )
  
  
  ggplot(ggplot_df, aes(x= variable, y = cluster, color = mean_expr, size = fraction)) + 
    geom_point() + 
    scale_color_gradientn(colors = colorscale,
                          trans = ifelse(logsc, 'log2', 'identity')) +
    xlab("Gene") +
    ylab("Cluster") +
    theme_classic() +
    theme(axis.text.y  = element_text( size = cex ),
          axis.text.x  = element_text(angle = 40, vjust = 1.1, hjust=1.1, size = cex),
          axis.title   = element_text( size = 18, face = "bold" ))
  
  
  }











