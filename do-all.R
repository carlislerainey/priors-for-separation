# set working directory to the priors-for-separation subdirectory
# setwd("~/Dropbox/projects/priors-for-separation")

# clear workspace
rm(list = ls())

# install all needed packages from CRAN and load
install.packages(c("arm", "coda", "logistf", "lme4", "xtable"))
library(arm)
library(coda)
library(logistf)

# install MCMCpack dependencies
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
biocLite("Rgraphviz")

# install packages from github and load 
devtools::install_github("carlislerainey/separation", ref = "b6a35bafe33b7850ca787e6e14a96ba961629da0")
devtools::install_github("carlislerainey/compactr", ref = "67a0df74d497a18f36affdc54533ef141c0e2bde")
library(separation)
library(compactr)

# create subdirectories (if missing) to store output
if (!file.exists("doc/figs")){
  dir.create(file.path("doc", "figs"))
}
if (!file.exists("doc/tabs")){
  dir.create(file.path("doc", "tabs"))
}
if (!file.exists("output")){
  dir.create("output")
}

# set seed
set.seed(8742570)

# theorem 1 illustrated
source("R/thm-1-illustrated.R")

# barrilleaux and rainey replication
source("br-replication/R/br.R")

# bell and miller replication
source("bm-replication/R/partial-prior.R")
source("bm-replication/R/mcmc.R")
source("bm-replication/R/plots.R")

