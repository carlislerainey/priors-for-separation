
# load and recode data
data(politics_and_need)  # load data set
d <- politics_and_need  # rename data set
d$dem_governor <- 1 - d$gop_governor  # create dem. gov. indicator
d$st_percent_uninsured <- rescale(d$percent_uninsured)  # standardize 

f <- oppose_expansion ~ dem_governor + percent_favorable_aca + gop_leg +
  st_percent_uninsured + bal2012 + multiplier + percent_nonwhite + percent_metro

m <- glm(f, data = d, family = binomial)


m.tol <- glm(f, data = d, family = binomial, control = list(epsilon = 10e-100))
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


###############################
## theorem 1 figures
###############################

shade <-  100
col1a <- rgb(shade, shade, shade, 150, maxColorValue = 255)
col2a <- rgb(shade, shade, shade, 150, maxColorValue = 255)
col1 <- rgb(shade, shade, shade, 255, maxColorValue = 255)
col2 <- rgb(shade, shade, shade, 255, maxColorValue = 255)

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


# normal + lik + post
pdf("doc/figs/thm-1-illustrated.pdf", height = 4, width = 6)
par(mfrow = c(2, 1), mar = c(.5, .5, .5, .5), oma = c(3, 1, 1, 1), xaxs = "r", yaxs = "r")
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Separating Variable",
      anny = FALSE)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
polygon(c(b1, rev(b1)), c(prior2/max(prior2), zeros), col = col1a, lty = 0)
text(3, .75, "Normal\nPrior", cex = 1, pos = 3, col = "black")
polygon(c(b1, rev(b1)), c(post2/max(post2), zeros), col = col2a, lty = 0)
text(b1[which(post2 == max(post2))], 1, "Posterior", cex = 1, pos = 3, col = "black")
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
#dev.off()

# cauchy + lik + post
eplot(xlim = 1.04*mm(b1), ylim = c(0, 1.2),
      xlab = "Coefficient for Separating Variable",
      anny = FALSE)
text(-17, 1, "Likelihood", cex = 1, pos = 3)
polygon(c(b1, rev(b1)), c(prior1/max(prior1), zeros), col = col1a, lty = 0)
text(3, .75, "Cauchy\nPrior", cex = 1, pos = 3, col = "black")
polygon(c(b1, rev(b1)), c(post1/max(post1), zeros), col = col2a, lty = 0)
text(b1[which(post1 == max(post1))], 1, "Posterior", cex = 1, pos = 3, col = "black")
lines(b1, lik/max(lik), lwd = 5, col = "black", xpd = NA)
dev.off()



