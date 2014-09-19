lp.scaled.t <- function(beta, scale.coef, scale.int, X, y) {
  p <- plogis(matrix(X%*%beta))
  loglik <- sum(y*log(p)) + sum((1-y)*log(1 - p)) + 
    log(dcauchy(beta[1], scale = scale.int)) + 
    sum(log(dcauchy(beta[2:length(beta)], scale = scale.coef)))
  return(loglik)
}

scaled.t <- function(f, data, df.coef = 7, scale.coef = 2.5,
                     df.int = 1, scale.int = 10, 
                   burnin = 500, mcmc = 1000, thin = 1,
                   tune = 1, verbose = 100) {
  require(MCMCpack)
  require(Matrix)
  require(arm)
  mf <- model.frame(f, data) 
  y <- model.response(mf)
  X <- model.matrix(f, data)
  X <- Matrix(X)
  init <- rnorm(ncol(X), 0, 1)
  cat("Computing proposal distribution...\n")
  mle <- bayesglm(f, data = d, family = binomial,
                  prior.scale.for.intercept = scale.int,
                  prior.scale = scale.coef,
                  prior.df.for.intercept = df.int,
                  prior.df = df.coef)
  V <- vcov(mle)
  print(summary(mle))
  seed <- round(runif(1, 0, 100000))
  mcmc <- MCMCmetrop1R(fun = lp.scaled.t, theta.init = init, 
                       scale.int = scale.int,
                       scale.coef = scale.coef,
                       X = X, y = y, burnin = burnin,
                       mcmc = mcmc, thin = thin, V = V,
                       tune = tune, seed = seed, 
                       verbose = verbose)
  return(mcmc)
}

# # test
# set.seed(1234)
# n <- 100
# x1 <- runif(n, -1, 1)
# X <- cbind(1, x1)
# beta <- c(-1, 1)
# y <- rbinom(n, 1, plogis(X%*%beta))
# d <- data.frame(x1, y)
# m1 <- scaled.t(y ~ x1, d, mcmc = 100, verbose = 10)
# plot(m1)

