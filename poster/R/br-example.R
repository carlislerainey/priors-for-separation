setwd("~/Dropbox/projects/priors-for-separation/poster")

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
pdf("figs/log-likelihood.pdf", height = 3, width = 5)
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

# normal + lik + post
pdf("figs/bm-post.pdf", height = 4, width = 6)
par(mfrow = c(2, 1), mar = c(.5, .5, .5, .5), oma = c(2, 0, 0, 0), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Democratic Governor",
      anny = FALSE)
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
polygon(c(b1, rev(b1)), c(prior2/max(prior2), zeros), col = col1a, lty = 0)
text(3, .75, "Normal\nPrior", cex = 1, pos = 3, col = col1)
polygon(c(b1, rev(b1)), c(post2/max(post2), zeros), col = col2a, lty = 0)
text(b1[which(post2 == max(post2))], 1, "Posterior", cex = 1, pos = 3, col = col2)

# cauchy + lik + post
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
