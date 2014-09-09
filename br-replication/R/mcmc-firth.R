rm(list = ls())

library(logistf)

n <- 10
x1 <- runif(n, -1, 1)
X <- cbind(1, x1)
beta <- c(0, 4)
y <- rbinom(n, 1, plogis(X%*%beta))


ll <- function(beta, X, y) {
  p <- plogis(X%*%beta)
  W <- matrix(0, nrow = length(p), ncol = length(p))
  diag(W) <- p*(1 - p)
  I <- t(X)%*%W%*%X
  loglik <- sum(y*log(p)) + sum((1-y)*log(1 - p)) + 0.5*log(det(I))
  return(loglik)
}

lprior <- function(beta, X, y) {
  p <- plogis(X%*%beta)
  W <- matrix(0, nrow = length(p), ncol = length(p))
  diag(W) <- p*(1 - p)
  I <- t(X)%*%W%*%X
  loglik <- 0.5*log(det(I))
  return(loglik)
}

b <- rep(1, 1, ncol(X))
optim(par = b, fn = ll, X = X, y = y, 
      control = list(fnscale = -1,
                     reltol = 1e-8))$par

logistf(y ~ x1)$coef

firth <- function(f, data, burnin = 500, mcmc = 1000, thin = 1,
                  tune = 1, verbose = 1000) {
  require(MCMCpack)
  # pull out the y and X
  mf <- model.frame(f, data) 
  y <- model.response(mf)
  X <- model.matrix(f, data)
  init <- rep(0, ncol(X))
  mcmc <- MCMCmetrop1R(fun = ll, theta.init = init, 
                       X = X, y = y, burnin = burnin,
                       mcmc = mcmc, thin = thin,
                       tune = tune, verbose = verbose)
  return(mcmc)
}

firth.prior <- function(f, data, burnin = 500, mcmc = 1000, thin = 1,
                  tune = 1, verbose = 1000) {
  require(MCMCpack)
  # pull out the y and X
  mf <- model.frame(f, data) 
  y <- model.response(mf)
  X <- model.matrix(f, data)
  init <- rep(0, ncol(X))
  mcmc <- MCMCmetrop1R(fun = lprior, theta.init = init, 
                       X = X, y = y, burnin = burnin,
                       mcmc = mcmc, thin = thin,
                       tune = tune, verbose = verbose)
  return(mcmc)
}

d <- data.frame(x1, y)
#m1 <- firth(y ~ x1, d, mcmc = 1000)
m2 <- firth.prior(y ~ x1, d, mcmc = 100, verbose = 10)

#plot(density(m2[, 2]))
#curve(dcauchy(x, 0, 2.5), col = "red", add = TRUE)
