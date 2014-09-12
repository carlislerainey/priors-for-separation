# set working directory
setwd("~/Dropbox/projects/priors-for-separation")

# clear working directory
rm(list = ls())

# load packages
library(parallel)

# function to iterate
do.it.once <- function(unused) {
  x <- unused
  source("bm-replication/R/mcmc.R")
}


# run chains
n.chains <- 4
mclapply(1:n.chains, do.it.once, mc.cores = 4)


