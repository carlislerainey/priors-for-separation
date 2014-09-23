# set working directory
setwd("~/Dropbox/projects/priors-for-separation")

# clear working directory
# rm(list = ls())

# load packages
library(coda)
library(arm)
library(devtools)
#install_github("carlislerainey/compactr")
library(compactr)


# create functions
## find median and 90% hpd and add to plot
sumry <- function(sims, ht, lab = NA) {
  q <- quantile(sims, c(.05, .95))
  hpd <- HPDinterval(mcmc(sims), prob = .9)[1:2]
  lines(hpd, c(ht, ht), lwd = 2)
  #points(hpd, c(ht, ht), pch = 4)
  d <- density(sims, adjust = 1.5)
  i.max <- which.max(d$y)
  post.mode <- d$x[i.max]
  post.mean <- mean(sims)
  post.median <- median(sims)
  #points(post.mode, ht, pch = 22, bg = "white")
  #points(post.median, ht, pch = 24, bg = "grey")
  points(post.median, ht, pch = 21, bg = "black", cex = .7)
  text(post.median, ht, lab, pos = 3, cex = .7)
  print(paste("HPD = ", round(hpd, 1)))
  print(paste("CI = ", round(q, 1)))
  print(paste("mode = ", round(post.mode, 1)))
  print(paste("median = ", round(post.median, 1)))
  print(paste("mean = ", round(post.mean, 1)))
  s <- round(c(hpd[1], post.median, hpd[2]), 1)
  return(s)
}
# find and return median and 90% hpd
num.sumry <- function(sims, ht, lab = NA) {
  hpd <- HPDinterval(mcmc(sims), prob = .9)[1:2]
  post.median <- median(sims)
  s <- round(c(hpd[1], post.median, hpd[2]), 1)
  return(s)
}
# plot posteior density and HPD for coefficients
plot.posterior.density <- function(sims) {
  dens <- density(sims)
  hpd <- HPDinterval(mcmc(sims), prob = .9)
  lo <- hpd[1]
  hi <- hpd[2]
  x <- dens$x[dens$x < hi & dens$x > lo]
  y <- dens$y[dens$x < hi & dens$x > lo]
  zeros <- rep(0, length(y))
  polygon(c(x, rev(x)), c(y, zeros), col = "grey80", lty = 0)
  abline(v = 0, col = "grey40")
  lines(dens, lwd = 2)
}
## simulate quantities of interest
simulate.qi <- function(x.hi, x.lo, mcmc) {
  p.hi <- plogis(c(x.hi%*%t(mcmc)))
  p.lo <- plogis(c(x.lo%*%t(mcmc)))
  fd <- p.hi - p.lo
  rr <- p.hi/p.lo
  ret <- list(p.hi = p.hi,
              p.lo = p.lo,
              fd = fd,
              rr = rr)
}

# load mcmc simulations
load("safe/cauchy1.RData"); m.gelman <- m.cauchy.mcmc
load("safe/firth1.RData"); m.zorn <- m.firth.mcmc
load("safe/m-me.RData"); m.inf <- m.me
load("safe/m-skep.RData")
load("safe/m-enth.RData")
rm(m.cauchy.mcmc, m.firth.mcmc, m.me)

# pull out twonukedyad coefficients
inf.sims <- m.inf$mcmc[, "c.twonukedyad"]
skep.sims <- m.skep$mcmc[, "c.twonukedyad"]
enth.sims <- m.enth$mcmc[, "c.twonukedyad"]
zorn.sims <- m.zorn[, 3]
gelman.sims <- m.gelman[, 3]

# calculate densities
d.inf <- density(inf.sims)
d.skep <- density(skep.sims)
d.enth <- density(enth.sims)
d.zorn <- density(zorn.sims)
d.gelman <- density(gelman.sims)

# simulate quantities of interest
# ## load and clean data
# d <- read.csv("bm-replication/data/bm.csv")
# d <- d[, c("warl2", "onenukedyad", "twonukedyad", "logCapabilityRatio", "Ally",
#            "SmlDemocracy", "SmlDependence", "logDistance", "Contiguity",
#            "MajorPower", "NIGOs")]
# d <- na.omit(d)
## rescale variables
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
## set formula
f <- warl2 ~ c.onenukedyad + c.twonukedyad + c.logCapabilityRatio + 
  c.Ally + c.SmlDemocracy + c.SmlDependence + c.logDistance + 
  c.Contiguity + c.MajorPower + c.NIGOs
## build model matrix
mf <- model.frame(f, d)
y <- model.response(mf)
X <- model.matrix(f, d)
x.nonukes <- x.twonukes <- apply(X, 2, median)
x.twonukes["c.twonukedyad"] <- max(d$c.twonukedyad)
## compute risk-ratios
inf.rr.sims <- simulate.qi(x.nonukes, x.twonukes, m.inf$mcmc)$rr
skep.rr.sims <- simulate.qi(x.nonukes, x.twonukes, m.skep$mcmc)$rr
enth.rr.sims <- simulate.qi(x.nonukes, x.twonukes, m.enth$mcmc)$rr
gelman.rr.sims <- simulate.qi(x.nonukes, x.twonukes, m.gelman)$rr
zorn.rr.sims <- simulate.qi(x.nonukes, x.twonukes, m.zorn)$rr

# coefficient plot 
pdf("doc/figs/bm-coef.pdf", height = 3.5, width = 6)
par(mfrow = c(1,1), mar = c(4,1,1,1), oma = c(0,0,0,0))
eplot(xlim = c(-15, 1), ylim = c(1.75, 6.5),
      anny = FALSE,
      xlabpos = 2.5,
      xlab = "Posterior Median and 90% HPD for\nCoefficient of Symmetric Nuclear Dyads")
abline(v = 0, col = "grey80")
sumry(inf.sims, 6, "Informative Normal(0, 4.5) Prior")
sumry(skep.sims, 5, "Skeptical Normal(0, 2) Prior")
sumry(enth.sims, 4, "Enthusiastic Normal(0, 8) Prior")
sumry(zorn.sims, 3, "Zorn's Default Jefferys' Invariant Prior")
sumry(gelman.sims, 2, "Gelman's Default Cauchy(0, 2.5) Prior")
dev.off()

# density plot
pdf("doc/figs/bm-posterior-density.pdf", height = 3.75, width = 8)
par(mfrow = c(2, 3), mar = c(1,1,1,1), oma = c(2,3,1,1))
eplot(xlim = c(-20, max(c(d.inf$x, d.skep$x, d.enth$x, d.zorn$x, d.gelman$x))), 
      mm(c(d.inf$y, d.skep$y, d.enth$y, d.zorn$y, d.gelman$y)),
      xlab = "Coefficient of Symmetric Nuclear Dyads",
      ylab = "Posterior Density",
      ylabpos = 2.5,
      main = "Informative Normal(0, 4.5) Prior")
plot.posterior.density(inf.sims)
aplot("Skeptical Normal(0, 2) Prior")
plot.posterior.density(skep.sims)
aplot("Enthusiastic Normal(0, 8) Prior")
plot.posterior.density(enth.sims)
addxaxis()
aplot("Zorn's Default Jeffreys' Prior")
plot.posterior.density(zorn.sims)
aplot("Gelman's Default Cauchy(0, 2.5) Prior")
plot.posterior.density(gelman.sims)
dev.off()

# risk-ratio plot 
pdf("doc/figs/bm-rr.pdf", height = 3.5, width = 6)
par(mfrow = c(1,1), mar = c(4,1,1,1), oma = c(0,0,0,0))
eplot(xlim = c(0.1, 1000000), ylim = c(1.75, 6.5),
      anny = FALSE,
      xlab = "Posterior Distribution of Risk-Ratio of War in Nonnuclear Dyads\nCompared to Symmetric Nuclear Dyads",
      xlabpos = 2.5,
      log = "x",
      xat = 10^(-1:6),
      xticklab = c("0.1", "1", "10", "100", "1,000", "10,000", "100,000", "1,000,000"))
abline(v = 1, col = "grey80")
s <- sumry(inf.rr.sims, 6, "Informative Normal(0, 4.5) Prior")
text(s[1], 5.85, s[1], cex = .5)
text(s[2], 5.85, s[2], cex = .5)
text(s[3], 5.85, s[3], cex = .5)
s <- sumry(skep.rr.sims, 5, "Skeptical Normal(0, 2) Prior")
text(s[1], 4.85, s[1], cex = .5)
text(s[2], 4.85, s[2], cex = .5)
text(s[3], 4.85, s[3], cex = .5)
s <- sumry(enth.rr.sims, 4, "Enthusiastic Normal(0, 8) Prior")
text(s[1], 3.85, s[1], cex = .5)
text(s[2], 3.85, s[2], cex = .5)
text(s[3], 3.85, s[3], cex = .5)
s <- sumry(zorn.rr.sims, 3, "Zorn's Default Jefferys' Prior")
text(s[1], 2.85, s[1], cex = .5)
text(s[2], 2.85, s[2], cex = .5)
text(s[3], 2.85, s[3], cex = .5)
s <- sumry(gelman.rr.sims, 2, "Gelman's Default Cauchy(0, 2.5) Prior")
text(s[1], 1.85, s[1], cex = .5)
text(s[2], 1.85, s[2], cex = .5)
text(s[3], 1.85, s[3], cex = .5)
dev.off()

# table summarizing hpd and median
S <- matrix(NA, nrow = 5, ncol = 3)
S[1, ] <- num.sumry(inf.rr.sims)
S[2, ] <- num.sumry(skep.rr.sims)
S[3, ] <- num.sumry(enth.rr.sims)
S[4, ] <- num.sumry(zorn.rr.sims)
S[5, ] <- num.sumry(gelman.rr.sims)
colnames(S) <- c("Lower-Bound of 90% HPD",
                 "Posterior Median",
                 "Upper-Bound of 90% HPD")
rownames(S) <- c("Informative Normal(0, 4.5) Prior",
                 "Skeptical Normal(0, 2) Prior",
                 "Enthusiastic Normal(0, 8) Prior",
                 "Zorn's Default Jeffreys' Prior",
                 "Gelman's Default Cauchy(0, 2.5) Prior")
pretty.S <- matrix(prettyNum(S, big.mark = ",", format = "fg", flag = " "), nrow = nrow(S)); pretty.S
rownames(pretty.S) <- rownames(S)
colnames(pretty.S) <- colnames(S)
library(xtable)
tab <- xtable(pretty.S, align = c("|", rep("c", ncol(S) + 1), "|"),
              caption = "This table provides lower- and upper-bounds for the 90% HPD and posterior medians 
              for the five prior distributions I consider.",
              label = "tab:bm-pppd-deciles")
print(tab, table.placement = "H", size = "scriptsize",
      file = "doc/tabs/bm-posterior-rr-summary.tex")

# posterior probability plot
labs <- rev(c("Informative Normal(0, 4.5) Prior",
  "Skeptical Normal(0, 2) Prior",
  "Enthusiastic Normal(0, 8) Prior",
  "Zorn's Default Jeffreys' Prior",
  "Gelman's Default Cauchy(0, 2.5) Prior"))
pdf("doc/figs/bm-pr-hypothesis.pdf", height = 3.5, width = 6)
par(mfrow = c(1,1), mar = c(4,1,1,1), oma = c(0,0,0,0), xaxs = "i")
eplot(xlim = c(0, 1), ylim = c(0.5, 5.5),
      anny = FALSE,
      xlab = "Pr(RR > 1)")
pr.h <- function(sims, ht) {
  p <- mean(sims > 1)
  lines(c(0, p), c(ht, ht), lwd = 3)
  points(p, ht, pch = 19, cex = .8)
  text(0, ht + .2, labs[ht], pos = 4, cex = .7)
  text(p, ht, round(p, 2), pos = 3, cex = .6)
}
pr.h(inf.rr.sims, 5)
pr.h(skep.rr.sims, 4)
pr.h(enth.rr.sims, 3)
pr.h(zorn.rr.sims, 2)
pr.h(gelman.rr.sims, 1)
dev.off()


