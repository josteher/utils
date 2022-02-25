# Plot multiple UMAP expression plots next to each other
# Works on sc objects from the RaceID package.


library(purrr)
library(cowplot)
library(magrittr)
library(ggplot2)
multi_plotexp <- function(sc_object, 
                          genes, 
                          logsc = T, 
                          scaled = F, 
                          cex = 0.2, 
                          # colorscale = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                          colorscale = c('black','darkblue', 'royalblue','turquoise', 'lightgreen', 'gold','orange', 'red2', 'red4'), # Moonlight rainbow
                          # colorscale = c('darkgrey', 'white', 'pink', 'purple'), # Purple
                          # colorscale = c('darkgrey', 'white', 'pink', 'red2'), # Red
                          # colorscale = c('darkgrey', 'white', 'lightblue', 'darkblue'), # Deep Arctic Sea
                          title_size = 15){
  
                plotmap_df <- sc_object@cpart %>% 
                  as.data.frame() %>%
                  rownames_to_column() %>%
                  # Bind columns since cells in sc_object@cpart and tsne should be ordered the same way
                  bind_cols(sc_object@umap) %>%
                  magrittr::set_colnames(c("cellid", "cluster", "X", "Y")) %>%
                  mutate("is_medoid"= cellid %in% sc_object@medoids, 
                         cluster = factor(cluster, levels = 1:length(unique(cluster)) ) ) 
                
                
                genes_df <- cbind(plotmap_df, 
                                  (sc_object@ndata[genes, rownames(plotmap_df), drop = F] * sc_object@counts) %>% 
                                    as.matrix() %>% t() + 0.1) %>% 
                                 data.table::data.table() %>%
                                 data.table::melt(., 
                                                  id.vars = c('cellid', 'cluster', 'X', 'Y', 'is_medoid'), 
                                                  measure.vars = genes)
                
                if (scaled){
                  gplot <- genes_df %>% 
                    ggplot(aes(X, Y, color = value )) +
                    geom_point(size = cex ) + # Adjust point size if more plots are plotted
                    scale_color_gradientn(colors = colorscale,
                                          trans = ifelse(logsc, 'log2', 'identity')  ) + # only change
                    theme_void() +
                    facet_wrap(facets = ~variable) +
                    coord_fixed(ratio = 1) +
                    theme(strip.text.x = element_text(size = title_size, 
                                                      face = "bold.italic", 
                                                      colour = "black", 
                                                      angle = 0,
                                                      vjust = 1),
                          legend.key.size = unit(1,"line"),
                          legend.key.width = unit(0.25,"line"),
                          legend.key.height = unit(0.7,"line")
                          ) 
                }

                else{
                 gplot <-  genes_df %>% 
                            group_split(variable) %>% 
                            purrr::map(~arrange(., value ) ) %>% # plot lowest expression first
                            purrr::map(
                              ~ggplot(. , aes(X, Y, color = value)) + 
                                geom_point(size = cex ) + # Adjust point size if more plots are plotted
                                scale_color_gradientn(colors = colorscale,
                                                      trans = ifelse(logsc, 'log2', 'identity')) +
                                theme_void()  +
                                facet_grid(~ variable, labeller = function(x) label_value(x, multi_line = FALSE)) +
                                coord_fixed(ratio = 1) +
                                theme(strip.text.x = element_text(size = title_size, 
                                                                  face = "bold.italic", 
                                                                  colour = "black", 
                                                                  angle = 0,
                                                                  vjust = 1),
                                      legend.key.size = unit(1,"line"),
                                      legend.key.width = unit(0.25,"line"),
                                      legend.key.height = unit(0.7,"line")
                                      ) 
                            ) %>% 
                            cowplot::plot_grid(plotlist = ., 
                                               align = 'hv', 
                                               ncol = length(genes) %>% sqrt() %>% round() )
                  
                }
                
                return(gplot)
                
}


# # Example:
#
# multi_plotexp(sc, c('Car1', 'Cxcl12', 'Cxcl9', 'Rag1', 'Il7r', 'Flt3', 'Mpo', 'Siglech', 'Kit'), 
#               scaled =F,       # show expression on common scale
#.              logsc =T,        # show as log2 expression
#               cex = 0.1,       # change pointsize
#               colorscale = c('black','darkblue', 'royalblue','turquoise', 'lightgreen', 'gold','orange', 'red2', 'red4') ) # colorscale
# 


