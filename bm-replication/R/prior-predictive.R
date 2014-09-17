# set working directory
setwd("~/Dropbox/projects/priors-for-separation")

# clear working directory
rm(list = ls())

# load packages
library(compactr)
library(arm)
library(metRology)

# set seed
# set.seed(8742570)

d <- read.csv("bm-replication/data/bm-small.csv")
d <- d[, c("warl2", "onenukedyad", "twonukedyad", "logCapabilityRatio", "Ally",
           "SmlDemocracy", "SmlDependence", "logDistance", "Contiguity",
           "MajorPower", "NIGOs")]
d <- na.omit(d)

# set formula (w/o twonekedyad)
f <- warl2 ~ onenukedyad + twonukedyad + logCapabilityRatio + 
  Ally + SmlDemocracy + SmlDependence + logDistance + 
  Contiguity + MajorPower + NIGOs

# read in pppd()
source("R/fn-pppd.R")

ps <- rnorm(10000, 0, 4.5)
pppd.out <- pppd(f = f, d = d, prior.sims = ps, s = "twonukedyad", s.at = 0, s.at.lo = FALSE)
