## Basic parameters
n <- 1000
B <- 1000
## Generate multivariate data
x <- cbind(runif(n,2,8),
           runif(n,2,3),
           runif(n,-3,3),
           runif(n,1,2),
           runif(n,2,6))
y <- 12+3.6*x[,1]+4*x[,2]-0.7*x[,3]+1.1*x[,4]-x[,5]+rnorm(n)

## Run linear model
lm <- lm(y~x)

## Collect residuals and fitted values
resid <- residuals(lm)
fit <- fitted(lm)

## Bootstrap function
boot <- function(B){
  
  ## Sample residuals
  resid.star <- resid[sample(seq(1:n),replace=TRUE)]
  
  ## Construct bootstrap y
  y.star <- fit + resid.star
  
  ## Run bootstrap regression
  lm.boot <- lm(y.star~x)
  
  ## Return bootstrap coefficients
  return(coefficients(lm.boot))
}

## Load library
library(snowfall)

## Determine available processors
detectCores()

## Set L'Ecuyer seed (runs from `parallel' package)
set.seed(123,"L'Ecuyer")

## Run parallel computation and record computation time
time.par <- system.time(mclapply(1:B,boot,mc.cores=cores)) 
  
## Shut down parallel environment
multicore:::closeAll()

## Run sequentially and record computation time
set.seed(123)
time.seq <- system.time(lapply(1:B,boot))