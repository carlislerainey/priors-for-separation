# clear workspace
rm(list = ls())
set.seed(761835)

# set working directory
setwd("~/Dropbox/projects/priors-for-separation/")

# read code 
source("br-replication/R/mcmc-firth.R")

# load packages
library(arm)
library(compactr)
library(brglm)
library(logistf)
library(texreg)
library(R2jags)

d <- read.csv("br-replication/data/need-rescale.csv")

################################################################################
## Estimate the Models
################################################################################


f <- oppose_expansion ~ gop_governor + obama_win + gop_leg + percent_uninsured + 
  income + percent_nonwhite + percent_metro

m2 <- firth.prior(f, d, mcmc = 1000000, verbose = 1000)

plot(density(m2[, 2]))
curve(dcauchy(x, 0, 2.5), col = "red", add = TRUE)

# maximum likelihood
m0 <- glm(f, family = binomial, data = d)

# jeffreys prior
#m.jeffreys <- brglm(f, family = binomial, data = d)
m1 <- logistf(f, alpha = 0.1, data = d)
#plot(profile(m.jeffreys, variable = "gop_governor"))

m2 <- firth(f, data = d, mcmc = 500000, burnin = 100000, thin = 25, verbose = 10000)

# cauchy prior + asymptotic
m3 <- bayesglm(f, family = binomial, data = d)

# cauchy prior + mcmc
#y <- rep(NA, 50)
y <- d$oppose_expansion
X <- model.matrix(f, d)
n <- length(y)
K <- ncol(X)
jags.data <- list("y", "X", "n", "K", "scale")
jags.params <- "beta"
jags.inits <- function() {
  list("beta" = rnorm(ncol(X), 0, 3))
}
scale <- 1
m4a <- jags(data = jags.data,
                      param = jags.params,
                      inits = jags.inits,
                      DIC = FALSE,
                      model = "br-replication/R/cauchy.bugs",
                      n.chains = 4,
                      n.iter = 50000)
plot(m4a)

scale <- 2.5
m4b <- jags(data = jags.data,
           param = jags.params,
           inits = jags.inits,
           DIC = FALSE,
           model = "br-replication/R/cauchy.bugs",
           n.chains = 4,
           n.iter = 50000)
plot(m4b)

scale <- 5
m4c <- jags(data = jags.data,
            param = jags.params,
            inits = jags.inits,
            DIC = FALSE,
            model = "br-replication/R/cauchy.bugs",
            n.chains = 4,
            n.iter = 50000)
plot(m4c)

rm(y, X, n, K, jags.data, jags.params, jags.inits)

## plot posteriors

pdf("doc/figs/matters-post.pdf", height = 3, width = 4)
par(mfrow = c(1,1), mar = c(3,4,1,1), oma = c(0,0,0,0))
eplot(xlim = c(-1, 30), ylim = c(0, .20),
      xlab = "Coefficient for GOP Governor",
      ylab = "Posterior Density",
      ylabpos = 2.5)
lines(density(m2[, 2], adjust = 1.5), lwd = 2, col = 1)
lines(density(m4a$BUGSoutput$sims.matrix[, 2]), lwd = 2, col = 2, lty = 1)
lines(density(m4b$BUGSoutput$sims.matrix[, 2]), lwd = 2, col = 3, lty = 1)
lines(density(m4c$BUGSoutput$sims.matrix[, 2]), lwd = 2, col = 4, lty = 1)
legend(x = par("usr")[2], y = par("usr")[4], xjust = 1, yjust = 1,
       legend = c("Jeffreys'",
                  "Cauchy(1)",
                  "Cauchy(2.5)",
                  "Cauchy(5)"),
       lty = 1,
       lwd = 2,
       col = 1:4, 
       cex = .8,
       bty = "n")
dev.off()

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

pdf("doc/figs/matters-ci.pdf", height = 3, width = 4)
par(mfrow = c(1,1), mar = c(3,1,1,1), oma = c(0,0,0,0))
eplot(xlim = c(-1, 36), ylim = c(.75, 4.5),
      anny = FALSE,
      xlab = "Coefficient for GOP Governor")
abline(v = 0, lty = 2)
# Jeffreys
sims <- m2[, 2]
sumry(sims, 4, "Jeffreys'")
# Cauchy(1)
sims <- m4a$BUGSoutput$sims.matrix[, 2]
sumry(sims, 3, "Cauchy(1)")
# Cauchy(2.5)
sims <- m4b$BUGSoutput$sims.matrix[, 2]
sumry(sims, 2, "Cauchy(2.5)")
# Cauchy(5)
sims <- m4c$BUGSoutput$sims.matrix[, 2]
sumry(sims, 1, "Cauchy(5)")

legend(x = par("usr")[2], y = par("usr")[4], xjust = 1, yjust = 1,
       legend = c("mode",
                  "median",
                  "mean"),
       pch = c(19, 21, 4),
       bg = "white", 
       cex = .8,
       bty = "n")
dev.off()
      