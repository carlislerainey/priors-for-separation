setwd("~/Dropbox/projects/priors-for-separation/present")

# load package
library(separation)
library(texreg)
library(logistf)
library(separation)
library(arm)

# load and recode data
data(politics_and_need)  # load data set
d <- politics_and_need  # rename data set
d$dem_governor <- 1 - d$gop_governor  # create dem. gov. indicator
d$st_percent_uninsured <- rescale(d$percent_uninsured)  # standardize 

f <- oppose_expansion ~ dem_governor + percent_favorable_aca + gop_leg +
  st_percent_uninsured + bal2012 + multiplier + percent_nonwhite + percent_metro

m <- glm(f, data = d, family = binomial)
display(m)
screenreg(list(m), single.row = TRUE, ci.force = TRUE)


m.tol <- glm(f, data = d, family = binomial, control = list(epsilon = 10e-100))
screenreg(list(m.tol), single.row = TRUE, ci.force = TRUE)

xtabs(~ dem_governor + oppose_expansion, data = d)

m.firth <- logistf(f, d)
summary(m.firth)

mf <- model.frame(f, data = d)
X <- model.matrix(mf, data = d)
y <- d$oppose_expansion
b <- coef(m)
b.dif <- b[2] - b[1]

ll_fn <- function(b, X, y) {
  p <- plogis(X%*%b)
  ll <- sum(y*log(p) + (1 - y)*log(1 - p))
  return(ll)
}

n_pts <- 200
b1 <- seq(-10, 1, length.out = n_pts)

ll <- numeric(n_pts)
for (i in 1:n_pts) {
  b.star <- b
  b.star[2] <- b1[i]
  ll[i] <- ll_fn(b.star, X, y)
}

library(compactr)
png("figs/log-likelihood.png", height = 3, width = 5, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0))
eplot(xlim = mm(b1), ylim = mm(ll),
      xlab = "Coefficient for Democratic Governor",
      ylab = "Log-Likelihood")
lines(b1, ll, lwd = 3)
dev.off()

###############################
## theorem 1 figures
###############################

col1a <- rgb(170, 86, 57, 150, maxColorValue = 255)
col2a <- rgb(39, 118, 80, 150, maxColorValue = 255)
col1 <- rgb(170, 86, 57, 255, maxColorValue = 255)
col2 <- rgb(39, 118, 80, 255, maxColorValue = 255)

lik_fn <- function(b, X, y) {
  p <- plogis(X%*%b)
  ll <- sum(y*log(p) + (1 - y)*log(1 - p))
  lik <- exp(ll)
  return(lik)
}

post_fn1 <- function(b, X, y) {
  p <- plogis(X%*%b)
  lp <- sum(y*log(p) + (1 - y)*log(1 - p)) + log(dcauchy(b[2], 0, 2.5))
  post <- exp(lp)
  return(post)
}

post_fn2 <- function(b, X, y) {
  p <- plogis(X%*%b)
  lp <- sum(y*log(p) + (1 - y)*log(1 - p)) + log(dnorm(b[2], 0, 2.5))
  post <- exp(lp)
  return(post)
}

n_pts <- 200
b1 <- seq(-20, 10, length.out = n_pts)
y <- d$oppose_expansion
b <- coef(m)
lik <- post1 <- post2 <- numeric(n_pts)
for (i in 1:n_pts) {
  b.star <- b
  b.star[2] <- b1[i]
  lik[i] <- lik_fn(b.star, X, y)
  post1[i] <- post_fn1(b.star, X, y)
  post2[i] <- post_fn2(b.star, X, y)  
}

prior1 <- dcauchy(b1, 0, 2.5)
prior2 <- dnorm(b1, 0, 2.5)
zeros <- rep(0, length(b1))

library(compactr)
# empty
png("figs/thm1-empty.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
dev.off()

# lik only
png("figs/thm1-lik-only.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
dev.off()

# cauchy only
png("figs/thm1-cauchy-only.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
polygon(c(b1, rev(b1)), c(prior1/max(prior1), zeros), col = col1a, lty = 0)
text(3, .75, "Cauchy\nPrior", cex = 1, pos = 3, col = col1)
dev.off()

# normal only
png("figs/thm1-normal-only.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
polygon(c(b1, rev(b1)), c(prior2/max(prior2), zeros), col = col1a, lty = 0)
text(3, .75, "Normal\nPrior", cex = 1, pos = 3, col = col1)
dev.off()

# normal + lik
png("figs/thm1-normal-lik.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
polygon(c(b1, rev(b1)), c(prior2/max(prior2), zeros), col = col1a, lty = 0)
text(3, .75, "Normal\nPrior", cex = 1, pos = 3, col = col1)
dev.off()

# normal + lik + post
png("figs/thm1-normal-lik-post.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
polygon(c(b1, rev(b1)), c(prior2/max(prior2), zeros), col = col1a, lty = 0)
text(3, .75, "Normal\nPrior", cex = 1, pos = 3, col = col1)
polygon(c(b1, rev(b1)), c(post2/max(post2), zeros), col = col2a, lty = 0)
text(b1[which(post2 == max(post2))], 1, "Posterior", cex = 1, pos = 3, col = col2)
dev.off()

# cauchy + lik
png("figs/thm1-cauchy-lik.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
polygon(c(b1, rev(b1)), c(prior1/max(prior1), zeros), col = col1a, lty = 0)
text(3, .75, "Cauchy\nPrior", cex = 1, pos = 3, col = col1)
dev.off()

# cauchy + lik + post
png("figs/thm1-cauchy-lik-post.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
polygon(c(b1, rev(b1)), c(prior1/max(prior1), zeros), col = col1a, lty = 0)
text(3, .75, "Cauchy\nPrior", cex = 1, pos = 3, col = col1)
polygon(c(b1, rev(b1)), c(post1/max(post1), zeros), col = col2a, lty = 0)
text(b1[which(post1 == max(post1))], 1, "Posterior", cex = 1, pos = 3, col = col2)
dev.off()

# cauchy-post-only
png("figs/thm1-cauchy-post-only.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
polygon(c(b1, rev(b1)), c(post1/max(post1), zeros), col = col2a, lty = 0)
text(b1[which(post1 == max(post1))], 1, "Posterior", cex = 1, pos = 3, col = col2)
dev.off()

# 2post
png("figs/thm1-2post.png", height = 2.5, width = 6, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
polygon(c(b1, rev(b1)), c(post1/max(post1), zeros), col = col2a, lty = 0)
text(b1[which(post1 == max(post1))], 1, "Posterior", cex = 1, pos = 3, col = col2)
polygon(c(b1, rev(b1)), c(post2/max(post2), zeros), col = col2a, lty = 0)
#text(b1[which(post2 == max(post2))], 1, "Posterior", cex = 1, pos = 3, col = col2)
dev.off()

library(compactr)
#png("figs/log-likelihood.png", height = 3, width = 5, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0))
# eplot(xlim = mm(b1), ylim = c(0, 1),
#       xlab = "Coefficient for Democratic Governor",
#       ylab = "Log-Likelihood")
# lines(b1, lik/max(lik), lwd = 3, col = col1)
prior <- dnorm(b1, sd = 2.5)
lines(b1, prior/max(prior), lwd = 3, lty = 3, col = col2)
lines(b1, post/max(post), lwd = 3, lty = 3, col = "purple")

#dev.off()

library(logistf)
m.firth <- logistf(f, d)

lpl_fn <- function(beta, X, y) {
  beta <- Matrix::Matrix(beta)
  p <- plogis(matrix(X%*%beta))
  W <- Matrix::Matrix(0, nrow = length(p), ncol = length(p))
  diag(W) <- p*(1 - p)
  I <- Matrix::crossprod(X, W)%*%X
  lp <- sum(y*log(p)) + sum((1-y)*log(1 - p)) + 0.5*log(Matrix::det(I))
  return(lp)
}

lpl <- numeric(n_pts)
for (i in 1:n_pts) {
  b.star <- b
  b.star[2] <- b1[i]
  lpl[i] <- lpl_fn(b.star, X, y)
}

library(compactr)
png("figs/log-penalized-likelihood.png", height = 3, width = 5, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1, 1), oma = c(0, 0, 0, 0))
eplot(xlim = mm(b1), ylim = mm(lpl),
      xlab = "Coefficient for Democratic Governor",
      ylab = "Log-Likelihood")
lines(b1, lpl, lwd = 3)
dev.off()

# cauchy prior + mcmc
#y <- rep(NA, 50)
data(politics_and_need_rescaled)
d_rsd <- politics_and_need_rescaled
d_rsd$dem_governor <- -d_rsd$gop_governor
y <- d_rsd$oppose_expansion
X <- model.matrix(f, d_rsd)
n <- length(y)
K <- ncol(X)
scale <- 2.5
jags.data <- list("y", "X", "n", "K", "scale")
jags.params <- "beta"
m.cauchy.mle <- glm(f, data = d_rsd, family = "binomial")
b.hat <- coef(m.cauchy.mle)
V.hat <- vcov(m.cauchy.mle)
jags.inits <- function() {
  list("beta" = mvrnorm(1, b.hat, 10*V.hat))
}

set.seed(1234)
library(R2jags)
# m.cauchy <- jags(data = jags.data,
#            param = jags.params,
#            inits = jags.inits,
#            DIC = FALSE,
#            model = "R/cauchy_logit.bugs",
#            n.chains = 3,
#            n.iter = 100)

mcmc <- m.cauchy$BUGSoutput$sims.matrix[, 2]
library(coda)
ci <- HPDinterval(as.mcmc(mcmc))
find_mode <- function(x) {
  dens <- density(x, n = 5000)
  mode <- dens$x[which(dens$y == max(dens$y))]
  mode <- median(mode)
  return(mode)
}
est <- find_mode(mcmc)

abs(ci[1] - ci[2])/abs(m.firth$ci.lower[2] - m.firth$ci.upper[2])


png("figs/cauchy-only.png", height = 3, width = 5, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0))
eplot(xlim = c(-28, 1), ylim = c(.2, 2.8),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
# Cauchy prior
points(est, 1, pch = 19)
lines(ci, c(1, 1), lwd = 3)
text(est, 1, "Cauchy", pos = 3, cex = .8)
dev.off()

png("figs/cauchy-jeffreys.png", height = 3, width = 5, units = "in", res = 300)
par(mfrow = c(1, 1), mar = c(3.5, 1, 1, 1), oma = c(0, 0, 0, 0))
eplot(xlim = c(-28, 1), ylim = c(.2, 2.8),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
# Cauchy prior
points(est, 1, pch = 19)
lines(ci, c(1, 1), lwd = 3)
text(est, 1, "Cauchy", pos = 3, cex = .8)
# Jeffreys' prior
points(m.firth$coefficients[2], 2, pch = 19)
lines(c(m.firth$ci.lower[2], m.firth$ci.upper[2]), c(2, 2), lwd = 3)
text(m.firth$coefficients[2], 2, "Jeffreys'", pos = 3, cex = .8)
dev.off()

# estimated RR for Cauchy
X_hi_list <- X_lo_list <- set_at_median(f, d_rsd)
X_hi_list$dem_governor <- max(d_rsd$dem_governor)
X_lo_list$dem_governor <- min(d_rsd$dem_governor)
X_hi <- list_to_matrix(X_hi_list, f)
X_lo <- list_to_matrix(X_lo_list, f)

b_hat <- coef(m.cauchy.mle)
p_hi <- plogis(X_hi%*%b_hat)
p_lo <- plogis(X_lo%*%b_hat)
rr_cauchy <- p_lo/p_hi

# estimated RR for Jeffreys'
X_hi_list <- X_lo_list <- set_at_median(f, d)
X_hi_list$dem_governor <- max(d$dem_governor)
X_lo_list$dem_governor <- min(d$dem_governor)
X_hi <- list_to_matrix(X_hi_list, f)
X_lo <- list_to_matrix(X_lo_list, f)

b_hat <- coef(m.firth)
p_hi <- plogis(X_hi%*%b_hat)
p_lo <- plogis(X_lo%*%b_hat)
rr_firth <- p_lo/p_hi


100*(rr_cauchy/rr_firth - 1)



