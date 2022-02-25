library(RaceID)
library(ggplot2)
library(ggrepel)
library(tidyverse)

# Caluculate differentially expressed genes and save result
# res data.frame structure
# rownames: genesymbols,
# colnames: baseMean, baseMeanA, baseMeanB, foldChange, log2FoldChange, pval, padj
cluster1 <- 8
cluster2 <- 2

A <- names(sc@cpart[sc@cpart %in% cluster1])
B <- names(sc@cpart[sc@cpart %in% cluster2])

res <- diffexpnb(getfdata(sc,n=c(A, B)), A=A, B=B )$res


#
MA_plot <- function(res){

          # Head and tail function
          head_tail <- function(x, top, bottom) {
            rbind(head(x, top), tail(x, bottom) ) %>% .[unique(rownames(.)), ]
          }
          
          
          res %>%
            rownames_to_column(var = "gene") %>%
            mutate(up_thr   = foldChange < 1, 
                   down_thr = foldChange > 1) %>%
            mutate(significance = ifelse( (padj < my_padj & up_thr), "sig_upin_A", 
                                          ifelse(
                                            (padj < my_padj & down_thr), "sig_downin_A", "ns") ) 
            ) %>%
            ggplot( aes(y = log2(baseMeanB) - log2(baseMeanA), x = log2(baseMeanA + baseMeanB), color = significance) ) +
            geom_point() +
            geom_hline( yintercept = 0, color= "black", linetype = "dashed") +
            geom_label_repel( aes( label = gene), 
                              data     = . %>% filter(padj < 0.05) %>% arrange(desc(abs(log2FoldChange))) %>% head_tail(50, 50), 
                              size     = 2.5, 
                              alpha    = 0.8, 
                              col      = "black", 
                              fontface = "bold.italic",
                              max.overlaps = 20) +
            scale_color_manual( values = c("grey","#B31B21", "#1465AC") ) +
            scale_y_continuous(breaks=seq(-200,200, 0.5)) +
            scale_x_continuous(breaks=seq(-200,200, 1)) +
            theme_classic()
}
