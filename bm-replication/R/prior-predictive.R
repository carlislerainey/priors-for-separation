# set working directory
setwd("~/Dropbox/projects/priors-for-separation")

# clear working directory
rm(list = ls())

# load packages
library(compactr)
library(arm)

# set seed
# set.seed(8742570)

d <- read.csv("bm-replication/data/bm.csv")
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

# my prior
ps <- rnorm(10000, 0, 4.5)
my.p <- pppd(f = f, d = d, prior.sims = ps, s = "twonukedyad", s.at = 0, s.at.lo = FALSE)

# skeptical prior
ps <- rnorm(10000, 0, 2)
skep.p <- pppd(f = f, d = d, prior.sims = ps, s = "twonukedyad", s.at = 0, s.at.lo = FALSE)

# enthusiastic prior
ps <- rnorm(10000, 0, 8)
enth.p <- pppd(f = f, d = d, prior.sims = ps, s = "twonukedyad", s.at = 0, s.at.lo = FALSE)

# gelman's prior
ps <- rcauchy(10000, 0, 2.5)
gel.p <- pppd(f = f, d = d, prior.sims = ps, s = "twonukedyad", s.at = 0, s.at.lo = FALSE)

x.max <- 100000
trunc.my.rr <- my.p$rr[my.p$rr < x.max]
trunc.gel.rr <- gel.p$rr[gel.p$rr < x.max]
trunc.skep.rr <- skep.p$rr[skep.p$rr < x.max]
trunc.enth.rr <- enth.p$rr[enth.p$rr < x.max]
my.p.trunc <- sum(my.p$rr > x.max)/length(my.p$rr)
gel.p.trunc <- sum(gel.p$rr > x.max)/length(gel.p$rr)
skep.p.trunc <- sum(skep.p$rr > x.max)/length(skep.p$rr)
enth.p.trunc <- sum(enth.p$rr > x.max)/length(enth.p$rr)


breaks <- seq(log(1), log(x.max), length.out = 20)
h1 <- hist(log(trunc.my.rr), breaks = breaks)
h2 <- hist(log(trunc.gel.rr), breaks = breaks)
h3 <- hist(log(trunc.skep.rr), breaks = breaks)
h4 <- hist(log(trunc.enth.rr), breaks = breaks)

add.arrow <- function(p) {
  x0 <- .88*log(x.max)
  x1 <- .98*log(x.max)
  ht <- 200
  arrows(x0 = x0, y0 = ht, x1 = x1, length = .1)
  text <- paste(round(100*p), "% of\nsimulations", sep = "")
  if (round(100*p) == 0) {
    text <- "< 1% of\nsimulations"
  }
  text((x0 + x1)/2, ht, text, pos = 3, cex = .7)
}

library(compactr)
xlim0 <- log(c(1, 100000))
ylim0 <- c(0, max(c(h1$counts, h2$counts, h3$counts, h4$counts)))
pdf("doc/figs/bm-pppd-hist.pdf", height = 4, width = 7)
par(mfrow = c(2,2), mar = c(.75,.75,.75,.75), oma = c(2,3,1,1))
eplot(xlim = xlim0, ylim = ylim0,
      xlab = "Risk-Ratio (Log Scale)",
      ylab = "Counts",
      ylabpos = 2.3,
      xat = log(c(1, 10, 100, 1000, 10000, 100000)),
      xticklab = c(1, 10, 100, 1000, 10000, 100000),
      main = "Informative Normal(0, 4.5) Prior")
plot(h1, add = TRUE, col = "grey50", border = NA)
add.arrow(my.p.trunc)

aplot("Gelman's Default Cauchy(0, 2.5) Prior")
plot(h2, add = TRUE, col = "grey50", border = NA)
add.arrow(gel.p.trunc)

aplot("Skeptical Normal(0, 2) Prior")
plot(h3, add = TRUE, col = "grey50", border = NA)
add.arrow(skep.p.trunc)

aplot("Enthusiastic Normal(0, 8) Prior")
plot(h4, add = TRUE, col = "grey50", border = NA)
add.arrow(enth.p.trunc)
dev.off()


Q <- round(rbind(my.p$Q[, "risk-ratio"],
      gel.p$Q[, "risk-ratio"],
      skep.p$Q[, "risk-ratio"],
      enth.p$Q[, "risk-ratio"]), 1)
rownames(Q) <- c("Informative Normal(0, 4.5) Prior",
                 "Gelman's Default Cauchy(0, 2.5) Prior",
                 "Skeptical Normal(0, 2) Prior",
                 "Enthusiastic Normal(0, 8) Prior")
pretty.Q <- matrix(prettyNum(Q, big.mark = ",", digits = 8, format = "fg", flag = " "), nrow = nrow(Q)); pretty.Q
rownames(pretty.Q) <- rownames(Q)
colnames(pretty.Q) <- colnames(Q)
library(xtable)
tab <- xtable(pretty.Q, align = c("|", rep("c", ncol(Q) + 1), "|"),
              caption = "This table provides the deciles prior predictive distribution for the 
                  risk-ratio of war in nonnuclear and nuclear dyads. The risk-ratio is 
                  tells us how many times more likely war is in non-nuclear dyads compared 
                  to nuclear dyads. Notice that that the, informative prior suggests a median 
                  risk-ratio of about 20, which is a large, but plausible effect. Gelman's 
                  proposed default prior suggests a slightly smaller median ratio of about 13, 
                  but also allows very large effects. The skeptical prior suggests a median 
                  ratio of about 4 and the enthusiastic prior suggests a median ratio of over 
                  200.",
              label = "tab:bm-pppd-deciles")
print(tab, table.placement = "H", size = "scriptsize",
      file = "doc/tabs/bm-pppd-deciles.tex")

