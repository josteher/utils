####
####
#### Hypergeometric test

# Extract knockout and wt cells
wt <- colnames(sc@ndata)[grep("WT|wt|Wt", colnames(sc@ndata))]
ko <- colnames(sc@ndata)[grep("KO|ko|Ko", colnames(sc@ndata))]

# Take a sample from wt, number of wt == number of ko
# wt_samp <- sample(wt,1201,replace = F)


# Population size
pop_sz    <- length(sc@cpart)
cenr <- list()
pv_thr <- 0.05


d <- data.frame(row.names = c( "#WT_cells", "#KO_cells", "fraction_of_wt", "pv"))
for (i in unique(sc@cpart)){
  print(i)
  # Sample size
  sample_sz <- length(sc@cpart[sc@cpart == i]) 
  
  # Number of successes in the population
  n_suc_pop <- sum( names(sc@cpart) %in% ko )
  
  # Number of successes in the sample
  n_suc_samp <- sum( names(sc@cpart[sc@cpart == i]) %in% ko )
  
  cenr[paste("C",i,sep="")] <- phyper(n_suc_samp,n_suc_pop, pop_sz - n_suc_pop, sample_sz )
  
  # Using Fisher test function
  nw  <- length(wt)
  nk  <- length(ko)
  niw <- sum( names(sc@cpart[sc@cpart == i]) %in% wt )
  nik <- sum( names(sc@cpart[sc@cpart == i]) %in% ko )
  
  ct <- matrix(c(nk - nik, nw - niw, nik, niw), ncol=2)
  pv <- fisher.test(ct)$p.value
  x <- c( niw, nik, niw/(niw + nik), pv)
  
  d[[as.character(i)]] <- x
  
}

# Check significance
sig_d <- d[,d["pv",] < pv_thr]

# Write Fishertest
# write.csv(d, "FisherTest_LSK_HSC_WT-KO.csv")
