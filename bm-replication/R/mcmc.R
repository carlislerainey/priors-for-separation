# clear working directory
rm(list = ls())

# set working directory
setwd("~/Dropbox/projects/priors-for-separation")

# load firth function
source("br-replication/R/mcmc-firth.R")
source("br-replication/R/mcmc-cauchy.R")


# load packages
library(compactr)
library(arm)
library(coda)

# set seed
set.seed(8742570)

d <- read.csv("bm-replication/data/bm.csv")
d <- d[, c("warl2", "onenukedyad", "twonukedyad", "logCapabilityRatio", "Ally",
           "SmlDemocracy", "SmlDependence", "logDistance", "Contiguity",
           "MajorPower", "NIGOs")]
d <- na.omit(d)

d$c.onenukedyad <- rescale(d$onenukedyad)
d$c.twonukedyad <- rescale(d$twonukedyad)
d$c.logCapabilityRatio <- rescale(d$logCapabilityRatio)
d$c.Ally <- rescale(d$Ally)
d$c.SmlDemocracy <- rescale(d$SmlDemocracy)
d$c.SmlDependence <- rescale(d$SmlDependence)
d$c.logDistance <- rescale(d$logDistance)
d$c.Contiguity<- rescale(d$Contiguity)
d$c.MajorPower <- rescale(d$MajorPower)
d$c.NIGOs <- rescale(d$NIGOs)

# set formula
f <- warl2 ~ c.onenukedyad + c.twonukedyad + c.logCapabilityRatio + 
  c.Ally + c.SmlDemocracy + c.SmlDependence + c.logDistance + 
  c.Contiguity + c.MajorPower + c.NIGOs

tune <- rep(.8, 11)
tune[3] <- 1.2

# cauchy prior + mcmc
m.cauchy.mcmc <- cauchy(f, d, mcmc = 50000, burnin = 10000, thin = 1, tune = tune, verbose = 100)
save(m.cauchy.mcmc, file = paste("output/cauchy", Sys.time()))

# firth prior + mcmc
m.firth.mcmc <- firth(f, d, mcmc = 50000, burnin = 10000, thin = 1, tune = tune, verbose = 100)
save(m.firth.mcmc, file = paste("output/firth", Sys.time()))


