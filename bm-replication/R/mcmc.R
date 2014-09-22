# set working directory
setwd("~/Dropbox/projects/priors-for-separation")

# clear working directory
rm(list = ls())

# load packages
library(compactr)
library(arm)
library(coda)

# set seed
# set.seed(8742570)

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



# load firth function
source("R/fn-mcmc-firth.R")
source("R/fn-mcmc-cauchy.R")
source("R/fn-mcmc-scaled-t.R")
source("R/fn-mcmc-normal.R")


# # my prior
# m.me <- normal(f, d, sep.var = "c.twonukedyad", sd = 4.5, 
#                n.sims = 40000, n.burnin = 10000, n.chains = 4)
# save(m.me, file = "safe/m-me.RData")
# 
# # # skeptical prior
# # m.skep <- normal(f, d, sep.var = "c.twonukedyad", sd = 2, 
# #                  n.sims = 1000, n.burnin = 100, n.chains = 4)
# 
# # enthusiastic prior
# m.enth <- normal(f, d, sep.var = "c.twonukedyad", sd = 8, 
#                  n.sims = 40000, n.burnin = 10000, n.chains = 4); plot(m.enth$mcmc.chains)
# save(m.enth, file = "safe/m-enth.RData")

# Zorn's default
#m.zorn <- firth(f, d, mcmc = 500, burnin = 100, thin = 1, tune = tune, verbose = 10)

# Gelman's default
#m.gelman <- cauchy(f, d, mcmc = 5000, burnin = 1000, thin = 1, tune = tune, verbose = 10)
