# clear workspace
rm(list = ls())

# set working directory
setwd("~/Dropbox/projects/priors-for-separation/")

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

# maximum likelihood
m0 <- glm(f, family = binomial, data = d)

# jeffreys prior
#m.jeffreys <- brglm(f, family = binomial, data = d)
m1 <- logistf(f, alpha = 0.1, data = d)
#plot(profile(m.jeffreys, variable = "gop_governor"))

# cauchy prior + asymptotic
m3 <- bayesglm(f, family = binomial, data = d)

# cauchy prior + mcmc
#y <- rep(NA, 50)
y <- d$oppose_expansion
X <- model.matrix(f, d)
n <- length(y)
K <- ncol(X)
scale <- 2.5
jags.data <- list("y", "X", "n", "K", "scale")
jags.params <- "beta"
jags.inits <- function() {
  list("beta" = rnorm(ncol(X), 0, 3))
}
m4 <- jags(data = jags.data,
                      param = jags.params,
                      inits = jags.inits,
                      DIC = FALSE,
                      model = "BR_Replication/R_Code/cauchy.bugs",
                      n.chains = 3,
                      n.iter = 1000)
plot(m4)

scale <- 5
m4b <- jags(data = jags.data,
           param = jags.params,
           inits = jags.inits,
           DIC = FALSE,
           model = "BR_Replication/R_Code/cauchy.bugs",
           n.chains = 3,
           n.iter = 1000)
plot(m4b)

scale <- 10
m4c <- jags(data = jags.data,
            param = jags.params,
            inits = jags.inits,
            DIC = FALSE,
            model = "BR_Replication/R_Code/cauchy.bugs",
            n.chains = 3,
            n.iter = 1000)
plot(m4c)

rm(y, X, n, K, jags.data, jags.params, jags.inits)


# # firth's penalty + bootstrap
# n.bs <- 1000
# m5 <- matrix(NA, nrow = n.bs, ncol = length(all.vars(f))) 
# pb <- txtProgressBar(min = 0, max = n.bs, style = 3)
# for (bs.iter in 1:n.bs) {
#   bs.samp <- sample(1:nrow(d), nrow(d), replace = TRUE)
#   bs.data <- d[bs.samp, ]
#   m5[bs.iter, ] <- coef(logistf(f, alpha = 0.1, data = d))
#   setTxtProgressBar(pb, bs.iter)
# } 
# rm(n.bs, bs.samp, bs.data)
# 
# # cauchy prior + bootstrap
# n.bs <- 1000
# m6 <- matrix(NA, nrow = n.bs, ncol = length(all.vars(f))) 
# pb <- txtProgressBar(min = 0, max = n.bs, style = 3)
# for (bs.iter in 1:n.bs) {
#   bs.samp <- sample(1:nrow(d), nrow(d), replace = TRUE)
#   bs.data <- d[bs.samp, ]
#   m6[bs.iter, ] <- coef(bayesglm(f, family = binomial, data = bs.data))
#   setTxtProgressBar(pb, bs.iter)
# } 
# rm(n.bs, bs.samp, bs.data)

################################################################################
## Plot the Estimates
################################################################################

model.names <- c("Maximum Likelihood", 
                 "Firth's Penalty\n+ Asymptotics", 
                 "Firth's Penalty\n+ Likelihood Profiling", 
                 "Firth's Penalty\n+ Bootstrap", 
                 "Cauchy Prior\n+ Asymptotics",
                 "Cauchy Prior\n+ Bootstrap",
                 "Cauchy Prior\n+ MCMC")
n.models <- length(model.names)

var.names <- c('gop_governor', 'obama_win', 'gop_leg', 'percent_uninsured',
               'income', 'percent_nonwhite', 'percent_metro', '(Intercept)')
bugs.coef.names <- c("beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]", 
                     "beta[7]", "beta[8]", "beta[1]")
label.var.names <- c('GOP Governor', 'Obama Win in 2012', 'GOP Controlled Legislature', 'Percent Uninsured',
                     'Income', 'Percent Nonwhite', 'Percent Metropolitan', 'Constant')
state.names <- sort(unique(d$state_abbr))

n.vars <- length(all.vars(f))
n.states <- length(state.names)

#svg("Figures/coefs.svg", height = 5, width = 10, family = "Georgia")
par(mfrow = c(2,4), oma = c(3,1,2,1), mar = c(.75, .75, 1, .5))
for (var.index in 1:n.vars) {
  eplot(xlim = c(-10, 50), ylim = c(-n.models - .2, -.2), main = label.var.names[var.index], anny = FALSE, 
        text.size = 1, ylab = "State Removed from the Analysis", ylabpos = 2.7,
        xlab = "Logistic Regression Coefficient\nand 90% Confidence Interval", xlabpos = 2.5)
  abline(v = 0, lty = 3)
 
  ### MLE
  est <- coef(m0)[var.names[var.index]]
  se <- sqrt(diag(vcov(m0)))[var.names[var.index]]
  lines(c(est + 1.64*se, est - 1.64*se), c(-1, -1), lty = 1.5)
  points(est, -1, pch = 19)
  xpos <- est
  if (abs(est) > 4) xpos <- 0
  text(xpos, -1, model.names[1], pos = 3)
  
  ### Firth's + Asymptotics
  est <- coef(m1)[var.names[var.index]]
  se <- sqrt(diag(vcov(m1)))
  names(se) <- names(coef(m1))
  se <- se[var.names[var.index]]
  lines(c(est + 1.64*se, est - 1.64*se), c(-2, -2), lty = 1.5)
  points(est, -2, pch = 19)
  text(est, -2, model.names[2], pos = 3)
  
  ### Firth's + Likelihood Profiling
  est <- coef(m1)[var.names[var.index]]
  lwr <- m1$ci.lower[var.names[var.index]]
  upr <- m1$ci.upper[var.names[var.index]]
  lines(c(lwr, upr), c(-3, -3), lty = 1.5)
  points(est, -3, pch = 19)
  text(est, -3, model.names[3], pos = 3)
  
  # a function to find the mode
  dmode <- function(x) {
    den <- density(x, kernel=c("gaussian"))
    ( den$x[den$y==max(den$y)] )   
  }
  
  ### Cauchy + MCMC
  est <- mean(m4c$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  est2 <- median(m4c$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  est3 <- dmode(m4c$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  lwr <- quantile(m4c$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]], .05)
  upr <- quantile(m4c$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]], .95)
  lines(c(lwr, upr), c(-7, -7), lty = 1.5)
  points(est, -7, pch = 19)
  points(est2, -7, pch = 21, bg = "white")
  points(est3, -7, pch = 4, bg = "white")
  text(est, -7, model.names[7], pos = 3)
  
  ### Second Cauchy + MCMC
  est <- mean(m4b$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  est2 <- median(m4b$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  est3 <- dmode(m4b$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  lwr <- quantile(m4b$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]], .05)
  upr <- quantile(m4b$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]], .95)
  lines(c(lwr, upr), c(-6, -6), lty = 1.5)
  points(est, -6, pch = 19)
  points(est2, -6, pch = 21, bg = "white")
  points(est3, -6, pch = 4, bg = "white")
  text(est, -6, model.names[7], pos = 3)
  
  ### Second Cauchy + MCMC
  est <- mean(m4$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  est2 <- median(m4$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  est3 <- dmode(m4$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]])
  lwr <- quantile(m4$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]], .05)
  upr <- quantile(m4$BUGSoutput$sims.matrix[, bugs.coef.names[var.index]], .95)
  lines(c(lwr, upr), c(-5, -5), lty = 1.5)
  points(est, -5, pch = 19)
  points(est2, -5, pch = 21, bg = "white")
  points(est3, -5, pch = 4, bg = "white")
  text(est, -5, model.names[7], pos = 3)
  
#   ### Cauchy + Bootstrap
#   q <- quantile(bs.res[, var.names[var.index]], c(.5, .05, .95))
#   lines(q[2:3], c(-1, -1), lty = 1.5)
#   points(q[1], -1, pch = 19)
#   text(q[1], -1, "Cauchy Prior\n+ Bootstrap", pos = 3)
#   
#   ### Cauchy
#   m <- bayesglm(bs.form, 
#                 prior.df = 1, prior.scale = 2.5,
#                 family = binomial, data = d)
#   est <- coef(m)[var.names[var.index]]
#   se <- sqrt(diag(vcov(m)))[var.names[var.index]]
#   lines(c(est + 1.64*se, est - 1.64*se), c(-2, -2), lty = 1.5)
#   points(est, -2, pch = 19)
#   text(est, -2, "Cauchy Prior", pos = 3)
#   
#   ### Firth's
#   m <- brglm(bs.form, family = binomial, data = d)
#   est <- coef(m)[var.names[var.index]]
#   se <- sqrt(diag(vcov(m)))[var.names[var.index]]
#   lines(c(est + 1.64*se, est - 1.64*se), c(-3, -3), lty = 1.5)
#   points(est, -3, pch = 19)
#   text(est, -3, "Firth's Logit", pos = 3)
#   
#   ### MLE
#   m <- glm(bs.form, family = binomial, data = d)
#   est <- coef(m)[var.names[var.index]]
#   se <- sqrt(diag(vcov(m)))[var.names[var.index]]
#   lines(c(est + 1.64*se, est - 1.64*se), c(-4, -4), lty = 1.5)
#   points(est, -4, pch = 19)
#   text(est, -4, "Maximum Likelihood", pos = 3)
}
#dev.off()

d1 <- density(m4$BUGSoutput$sims.matrix[, "beta[2]"], from = -10, to = 50)
d2 <- density(m4b$BUGSoutput$sims.matrix[, "beta[2]"], from = -10, to = 50)
d3 <- density(m4c$BUGSoutput$sims.matrix[, "beta[2]"], from = -10, to = 50)

#mm(c(0, d1$x, d2$x, d3$x))
par(mfrow = c(1,1), oma = c(3,3,1,1))
eplot(xlim = c(-10, 50), ylim = mm(c(0, d1$y, d2$y, d3$y)))
lines(d1, lwd = 3, lty = 1)
curve(dcauchy(x, scale = 2.5), lty = 1, col = "red", add = TRUE)
lines(d2, lwd = 3, lty = 2)
curve(dcauchy(x, scale = 5), lty = 2, col = "red", add = TRUE)
lines(d3, lwd = 3, lty = 3)
curve(dcauchy(x, scale = 10), lty = 3, col = "red", add = TRUE)


