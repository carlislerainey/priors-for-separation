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

d <- read.csv("bm-replication/data/bm-small.csv")
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

# firth prior + mcmc
m.firth.mcmc <- firth(f, d, mcmc = 100, burnin = 10, thin = 1, tune = tune, verbose = 100)

#pdf("doc/figs/bm-density.pdf", height = 3, width = 4)
par(mfrow = c(1,1), mar = c(3,4,1,1), oma = c(0,0,0,0))
d.firth <- density(m.firth.mcmc[, 3], adjust = 1.5)
d.cauchy <- density(m.cauchy.mcmc[, 3], adjust = 1.5)
eplot(xlim = mm(c(d.firth$x, d.cauchy$x)), mm(c(d.firth$y, d.cauchy$y)),
      xlab = "Coefficient for Two-Nuke Dyad",
      ylab = "Posterior Density",
      ylabpos = 2.5)
lines(d.firth, lwd = 2, col = 1)
lines(d.cauchy, lwd = 2, col = 2, lty = 1)
legend(x = par("usr")[1], y = par("usr")[4], xjust = 0, yjust = 1,
       legend = c("Jeffreys'",
                  "Cauchy(2.5)"),
       lty = 1,
       lwd = 2,
       col = 1:2, 
       cex = .8,
       bty = "n")
#dev.off()

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

#pdf("doc/figs/bm-ci.pdf", height = 3, width = 4)
par(mfrow = c(1,1), mar = c(3,1,1,1), oma = c(0,0,0,0))
eplot(xlim = c(-15, 1), ylim = c(.75, 4.5),
      anny = FALSE,
      xlab = "Coefficient for Two-Nuke Dyad")
abline(v = 0, lty = 2)
# Jeffreys
sims <- m.firth.mcmc[, 3]
sumry(sims, 4, "Jeffreys'")
# Cauchy(2.5)
sims <- m.cauchy.mcmc[, 3]
sumry(sims, 2, "Cauchy(2.5)")

legend(x = par("usr")[2], y = par("usr")[4], xjust = 1, yjust = 1,
       legend = c("mode",
                  "median",
                  "mean"),
       pch = c(19, 21, 4),
       bg = "white", 
       cex = .8,
       bty = "n")
#dev.off()