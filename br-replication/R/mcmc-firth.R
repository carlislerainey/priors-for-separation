lp <- function(beta, X, y) {
  p <- plogis(matrix(X%*%beta))
  W <- Matrix(0, nrow = length(p), ncol = length(p))
  diag(W) <- p*(1 - p)
  I <- crossprod(X, W)%*%X
  loglik <- sum(y*log(p)) + sum((1-y)*log(1 - p)) + 0.5*log(det(I))
  return(loglik)
}

firth <- function(f, data, burnin = 500, mcmc = 1000, thin = 1,
                  tune = 1, verbose = 1000) {
  require(MCMCpack)
  require(Matrix)
  mf <- model.frame(f, data) 
  y <- model.response(mf)
  X <- model.matrix(f, data)
  X <- Matrix(X)
  init <- rep(0, ncol(X))
  mcmc <- MCMCmetrop1R(fun = lp, theta.init = init, 
                       X = X, y = y, burnin = burnin,
                       mcmc = mcmc, thin = thin,
                       tune = tune, verbose = verbose,
                       optim.control = list(reltol = 1e-2))
  return(mcmc)
}

# # test
# set.seed(1234)
# n <- 1000
# x1 <- runif(n, -1, 1)
# X <- cbind(1, x1)
# beta <- c(-1, 1)
# y <- rbinom(n, 1, plogis(X%*%beta))
# d <- data.frame(x1, y)
# m1 <- firth(y ~ x1, d, mcmc = 1000, verbose = 10)
# plot(m1)

