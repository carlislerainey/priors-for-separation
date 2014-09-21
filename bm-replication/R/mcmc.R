# # set working directory
# setwd("~/Dropbox/projects/priors-for-separation")
# 
# # clear working directory
# rm(list = ls())
# 
# # load packages
# library(compactr)
# library(arm)
# library(coda)
# 
# # set seed
# # set.seed(8742570)
# 
# d <- read.csv("bm-replication/data/bm.csv")
# d <- d[, c("warl2", "onenukedyad", "twonukedyad", "logCapabilityRatio", "Ally",
#            "SmlDemocracy", "SmlDependence", "logDistance", "Contiguity",
#            "MajorPower", "NIGOs")]
# d <- na.omit(d)
# 
# d$c.onenukedyad <- rescale(d$onenukedyad)
# d$c.twonukedyad <- rescale(d$twonukedyad)
# d$c.logCapabilityRatio <- rescale(d$logCapabilityRatio)
# d$c.Ally <- rescale(d$Ally)
# d$c.SmlDemocracy <- rescale(d$SmlDemocracy)
# d$c.SmlDependence <- rescale(d$SmlDependence)
# d$c.logDistance <- rescale(d$logDistance)
# d$c.Contiguity<- rescale(d$Contiguity)
# d$c.MajorPower <- rescale(d$MajorPower)
# d$c.NIGOs <- rescale(d$NIGOs)
# 
# # set formula
# f <- warl2 ~ c.onenukedyad + c.twonukedyad + c.logCapabilityRatio + 
#   c.Ally + c.SmlDemocracy + c.SmlDependence + c.logDistance + 
#   c.Contiguity + c.MajorPower + c.NIGOs
# 


# load firth function
source("R/fn-mcmc-firth.R")
source("R/fn-mcmc-cauchy.R")
source("R/fn-mcmc-scaled-t.R")
source("R/fn-mcmc-normal.R")


# my prior
#m.me <- normal(f, d, sep.var = "c.twonukedyad", sd = 4.5, n.sims = 20000, n.burnin = 10000, verbose = 100, n.chains = 4)

# # skeptical prior
# m.skep <- normal(f, d, sep.var = "c.twonukedyad", sd = 2, 
#                  n.sims = 20000, n.burnin = 10000, 
#                  verbose = 100, n.chains = 4)

# enthusiastic prior
m.enth <- normal(f, d, sep.var = "c.twonukedyad", sd = 8, 
                 n.sims = 10000, n.burnin = 1000, 
                 verbose = 100, n.chains = 4)

# Zorn's default
#m.zorn <- firth(f, d, mcmc = 500, burnin = 100, thin = 1, tune = tune, verbose = 10)

# Gelman's default
#m.gelman <- cauchy(f, d, mcmc = 5000, burnin = 1000, thin = 1, tune = tune, verbose = 10)

# load mcmc simulations
load("safe/cauchy1.RData"); m.gelman <- m.cauchy.mcmc
load("safe/firth1.RData"); m.zorn <- m.firth.mcmc
load("safe/m-skep.RData")
rm(m.cauchy.mcmc, m.firth.mcmc)

###########
## plots ##
###########

my.sims <- m.me$mcmc[, "c.twonukedyad"]
skep.sims <- m.skep$mcmc[, "c.twonukedyad"]
enth.sims <- m.enth$mcmc[, "c.twonukedyad"]

zorn.sims <- m.zorn[, 3]
gelman.sims <- m.gelman[, 3]

par(mfrow = c(1,1), mar = c(3,4,1,1), oma = c(0,0,0,0))

d.me <- density(my.sims)
d.skep <- density(skep.sims)
d.enth <- density(enth.sims)
d.zorn <- density(zorn.sims)
d.gelman <- density(gelman.sims)

eplot(xlim = c(-20, max(c(d.me$x, d.zorn$x, d.gelman$x))), mm(c(d.me$y, d.zorn$y, d.gelman$y)),
      xlab = "Coefficient for Two-Nuke Dyad",
      ylab = "Posterior Density",
      ylabpos = 2.5)
abline(v = 0, col = "grey50")
lines(d.me, lwd = 2, col = 1)
lines(d.skep, lwd = 2, col = 2)
lines(d.enth, lwd = 2, col = 3)
lines(d.zorn, lwd = 2, col = 4)
lines(d.gelman, lwd = 2, col = 5)
legend(x = par("usr")[1], y = par("usr")[4], xjust = 0, yjust = 1,
       legend = c("Informative Prior",
                  "Skeptical Prior",
                  "Enthusiastic Prior",
                  "Zorn's Default Prior",
                  "Gelman's Default Prior"),
       lty = 1,
       lwd = 2,
       col = 1:5, 
       cex = .8,
       bty = "n")

sumry <- function(sims, ht, lab = NA) {
  q <- quantile(sims, c(.05, .95))
  lines(q, c(ht, ht), lwd = 2)
  d <- density(sims, adjust = 1.5)
  i.max <- which.max(d$y)
  post.mode <- d$x[i.max]
  post.mean <- mean(sims)
  post.median <- median(sims)
  points(post.mode, ht, pch = 19)
  points(post.median, ht, pch = 21, bg = "white")
  points(post.mean, ht, pch = 4)
  text(post.mean, ht, lab, pos = 3, cex = .7)
  print(paste("CI = ", round(q, 1)))
  print(paste("mode = ", round(post.mode, 1)))
  print(paste("median = ", round(post.median, 1)))
  print(paste("mean = ", round(post.mean, 1)))
  
}

par(mfrow = c(1,1), mar = c(3,1,1,1), oma = c(0,0,0,0))
eplot(xlim = c(-15, 1), ylim = c(.75, 6.5),
      anny = FALSE,
      xlab = "Coefficient for Two-Nuke Dyad")
abline(v = 0, lty = 2)
sumry(my.sims, 6, "my prior")
sumry(skep.sims, 5, "skeptical prior")
sumry(enth.sims, 4, "enthusiastic prior")
sumry(zorn.sims, 3, "Zorn's default")
sumry(gelman.sims, 2, "Gelman's default")


legend(x = par("usr")[2], y = par("usr")[4], xjust = 1, yjust = 1,
       legend = c("mode",
                  "median",
                  "mean"),
       pch = c(19, 21, 4),
       bg = "white", 
       cex = .8,
       bty = "n")
