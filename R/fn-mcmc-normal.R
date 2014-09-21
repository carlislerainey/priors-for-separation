# rm(list = ls())
# 
# # test
# set.seed(1234)
# n <- 100
# x1 <- runif(n, -1, 1)
# X <- cbind(1, x1)
# beta <- c(-1, 1)
# y <- rbinom(n, 1, plogis(X%*%beta))
# d <- data.frame(x1, y)
# f <- y ~ x1
# data <- d
# n.chains <- 3
# sep.var <- "x1"
# sd <- 4


lp.normal <- function(beta, which.var.sep, sd, X, y) {
  p <- plogis(matrix(X%*%beta))
  loglik <- sum(y*log(p)) + sum((1-y)*log(1 - p))
  logprior <- log(dnorm(beta[which.var.sep], mean = 0, sd = sd))
  logpost <- loglik + logprior
  return(logpost)
}

normal <- function(f, data, sep.var, sd = 4.5,
                   n.sims = 1000, n.burnin = 500,
                   n.thin = 1, n.chains = 3,
                   tune = 1) {
  require(MCMCpack)
  require(Matrix)
  require(arm)
  require(parallel)
  mf <- model.frame(f, data) 
  y <- model.response(mf)
  X <- model.matrix(f, data)
  X <- Matrix(X)
  if (sum(colnames(X) == sep.var) < 1) {
   stop(paste(sep.var,  " is not a variable in the model", sep = ""))
  }
  which.var.sep <- which(colnames(X) == sep.var)
  init <- rnorm(ncol(X), 0, 2)
  cat("\nComputing proposal distribution...\n")
  mle <- bayesglm(f, data = d, family = binomial)
  V <- vcov(mle)
  #print(summary(mle))
  mcmc.chains <- mcmc.list()
  mcmc <- NULL
  l.seed <- runif(6, 100000, 999999)
  init.seed <- runif(n.chains, 100000, 999999)
  run.mcmc <- function(x) {
    set.seed(init.seed[x])
    init <- rnorm(ncol(X), coef(mle), 2)
    mcmc <- MCMCmetrop1R(fun = lp.normal, 
                         theta.init = init,  V = V,
                         sd = sd, which.var.sep = which.var.sep, X = X, y = y, 
                         thin = n.thin, burnin = n.burnin, mcmc = n.sims,
                         tune = tune, verbose = 0,
                         seed = list(l.seed, x))
    return(mcmc)
  }
  cat(paste("\nRunning ", n.chains, " chains in parallel of ", n.sims + n.burnin, " iterations each--this may take a while...", sep = ""))
  mcmc.chains <- mclapply(1:n.chains, run.mcmc, mc.cores = 3)
  cat(paste("\nFinished running chains!\n", sep = ""))
  mcmc.chains <- as.mcmc.list(mcmc.chains)
  for (i in 1:n.chains) {  
    mcmc <- rbind(mcmc, mcmc.chains[[i]])
  }
  colnames(mcmc) <- colnames(X)
  R.hat <- gelman.diag(mcmc.chains)
  cat(paste("\nChecking convergence...\n", sep = ""))
  if (R.hat[[2]] <= 1.02) {
    cat(paste("\nThe multivariate R-hat statistic of ", round(R.hat[[2]], 2), 
              " suggests that the chains have converged.\n\n", sep = ""))
  }
  if (R.hat[[2]] > 1.02) {
    cat(paste("\n######## WARNING: #########\n\nThe multivariate R-hat statistic of ", round(R.hat[[2]], 2), 
              " suggests that the chains have NOT converged.\n\n", sep = ""))
  }
  res <- list(mcmc.chains = mcmc.chains,
              mcmc = mcmc,
              R.hat = R.hat)
  return(res)
}

# # test
# set.seed(1234)
# n <- 100
# x1 <- rbinom(n, 1, .5)
# X <- cbind(1, x1)
# beta <- c(-1, 1)
# y <- rbinom(n, 1, plogis(X%*%beta))
# y[x1 == 1] <- 1
# d <- data.frame(x1, y)
# m1 <- normal(y ~ x1, d, sep.var = "x1", n.sims = 1000, n.burnin = 10, verbose = 0, n.chains = 4)
# plot(m1$mcmc.chains)

