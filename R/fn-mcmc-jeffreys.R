lp.jeffreys <- function(beta, X, y) {
  p <- plogis(matrix(X%*%beta))
  W <- Matrix(0, nrow = length(p), ncol = length(p))
  diag(W) <- p*(1 - p)
  I <- crossprod(X, W)%*%X
  loglik <- sum(y*log(p)) + sum((1-y)*log(1 - p)) + 0.5*log(det(I))
  return(loglik)
}

jeffreys <- function(f, data, n.sims = 1000, n.burnin = 100, n.thin = 1,
                  tune = 1, n.chains = 3, n.cores = n.chains) {
  require(MCMCpack)
  require(Matrix)
  require(logistf)
  require(parallel)
  mf <- model.frame(f, data) 
  y <- model.response(mf)
  X <- model.matrix(f, data)
  X <- Matrix(X)
  cat("\nComputing proposal distribution...\n")
  mle <- logistf(f, d)
  V <- vcov(mle)
  #print(summary(mle))
  mcmc.chains <- mcmc.list()
  mcmc <- NULL
  l.seed <- runif(6, 100000, 999999)
  init.seed <- runif(n.chains, 100000, 999999)
  run.mcmc <- function(x) {
    set.seed(init.seed[x])
    init <- rnorm(ncol(X), coef(mle), 1)
    mcmc <- MCMCmetrop1R(fun = lp.jeffreys, 
                         theta.init = init,  V = V,
                         X = X, y = y, 
                         thin = n.thin, burnin = n.burnin, mcmc = n.sims,
                         tune = tune, verbose = 0,
                         seed = list(l.seed, x))
    return(mcmc)
  }
  cat(paste("\nRunning ", n.chains, " chains in parallel of ", n.sims + n.burnin, " iterations each--this may take a while...", sep = ""))
  mcmc.chains <- mclapply(1:n.chains, run.mcmc, mc.cores = n.cores)
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

# test
set.seed(1234)
n <- 100
x1 <- rbinom(n, 1, .5)
X <- cbind(1, x1)
beta <- c(-1, 1)
y <- rbinom(n, 1, plogis(X%*%beta))
y[x1 == 1] <- 1
d <- data.frame(x1, y)
m1 <- jeffreys(y ~ x1, d, n.sims = 100, n.burnin = 10, n.chains = 4)
plot(m1$mcmc.chains)

