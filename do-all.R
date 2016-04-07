# set working directory to the priors-for-separation subdirectory
# setwd("~/Dropbox/projects/priors-for-separation")

# clear workspace
rm(list = ls())

# paramaters
install_packages <- FALSE
n_sims <- 200000  # 200,000 in paper---use 200 as the initial run
n_burnin <- 50000 # 50,0000 in paper---use 50 as the initial run
n_chains <- 4  # 4 in paper

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

# install packages
if (install_packages) {
  # install all needed packages from CRAN
  install.packages(c("devtools", "arm", "coda", "logistf", "lme4", "xtable"))
  # install MCMCpack dependencies not on CRAN
  source("https://bioconductor.org/biocLite.R")
  biocLite("graph")
  biocLite("Rgraphviz")
  # install packages from github
  devtools::install_github("carlislerainey/separation", ref = "b6a35bafe33b7850ca787e6e14a96ba961629da0")
  devtools::install_github("carlislerainey/compactr", ref = "67a0df74d497a18f36affdc54533ef141c0e2bde")
}

# load packages
library(arm)
library(coda)
library(logistf)
library(separation)
library(compactr)
library(xtable)

# set seed
set.seed(8742570)

# theorem 1 illustrated 
source("R/thm-1-illustrated.R")  # note: glm.fit warning expected

# barrilleaux and rainey replication
source("br-replication/R/br.R")

# bell and miller replication
source("bm-replication/R/partial-prior.R")  # (note: glm.fit warning expected)
source("bm-replication/R/mcmc.R")  # (note: glm.fit warning expected)
source("bm-replication/R/plots.R")

