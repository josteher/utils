library(RaceID)
library(doParallel)
library(purrr)




# Provide an sc object with cluster partitions

#####
#####  cdiffgenes on multiple cores
#####
# Register number of cores for parallel computing using doParallel
# Specify cores according to your system
registerDoParallel(cores =  detectCores() %/% 6 )
s <- Sys.time()

# Specify unique clusters
uq_cluster        <- sort(unique(sc@cpart))
names(uq_cluster) <- paste("cl", uq_cluster, sep=".")

cl_diff <- foreach(i = uq_cluster, .final = function(x) stats::setNames(x, names(uq_cluster))  ) %dopar% {
  clustdiffgenes(sc, i, pvalue=0.01)
}
# identical(cdiff, cl_diff)
Sys.time() - s





#####
##### Pairwise comparison of clusters to one another using diffexpnb
#####
##### Choose clusters to compare
my_clusters <- unique(sc@cpart)
# my_clusters <- c(1,2,3,4,5,6)
comb_cl <- combn(my_clusters,2)
colnames(comb_cl) <- unlist(lapply(1:dim(comb_cl)[2], function(i) {paste(comb_cl[[1,i]], comb_cl[[2,i]] ,sep = ".v.") })) 

# Parallel computation of differnetail expressed genes in pairwise Comparisons using doParallel
# Specify cores according to your system
registerDoParallel(cores =  detectCores() %/% 6 )

s <- Sys.time()

d <- foreach(i = 1:dim(comb_cl)[2], .final = function(x) setNames(x, colnames(comb_cl)) ) %dopar% {
  cluster1 <- comb_cl[[1,i]]
  cluster2 <- comb_cl[[2,i]]
  
  cl1 <- names(sc@cpart[sc@cpart %in% cluster1])
  cl2 <- names(sc@cpart[sc@cpart %in% cluster2])
  
  list(diffexpnb(getfdata(sc,n=c(cl1,cl2)), A=cl1, B=cl2 ) )
}
# names(d) <- unlist(lapply(1:dim(comb_cl)[2], function(i) {paste(comb_cl[[1,i]], comb_cl[[2,i]] ,sep = ".vs.") }))
fw_dexl <- lapply(d, function(x) {x[[1]]  })

Sys.time() - s


# Function to obtain output from diffexpnb when reversed comparison is desired e.g. comparison of 1 vs 2  =>  2 vs 1
rev_diffexp_comp <- function(x){
  res                <- x$res
  
  res$baseMeanA      <- x$res$baseMeanB
  res$baseMeanB      <- x$res$baseMeanA
  res$foldChange     <- 1/res$foldChange
  res$log2FoldChange <- log2(res$foldChange)
  
  return(list(vf1 = x$vf2 , vf2 = x$vf1 , res = res ))
}

# Create list of reversed comparisons and name them
rev_dexl        <- lapply(fw_dexl, rev_diffexp_comp)
names(rev_dexl) <- stringi::stri_reverse(names(fw_dexl))

# Concatenate into a common list
dexl <- c(fw_dexl, rev_dexl) 


# up_genes   <- lapply(dexl, function(x) {x$res[x$res$padj < 0.05 & x$res$foldChange > 1.1,]}  )
# down_genes <- lapply(dexl, function(x) {x$res[x$res$padj < 0.05 & x$res$foldChange < 1/1.1,]}  )
# 

# Create one data frame for all pairwise comparisons
dex_df <- map(dexl, "res") %>%
            map(rownames_to_column, var = "gene") %>%
            bind_rows(.id = "comparison") %>%
            separate(comparison, c("Cluster_A", "Cluster_B"), sep=".v.", remove = F) %>%
            mutate(upin_A = foldChange < 1, upin_B = foldChange > 1) %>%
            # Count how many times a gene is upregulated in cluster A in all comparisons
            group_by(Cluster_A, gene, upin_A) %>%
            mutate(count_gene_upin_A = n() )



# SAve differential gene expression result to res
cl1 <- 1
cl2 <- 2
res <- dexl[[paste(cl1, ".v.", cl2, sep="")]][["res"]] 








