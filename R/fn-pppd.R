
pppd <- function(f, d, prior.sims, s, s.at, s.at.lo = TRUE) {
  # f - a logitistic regression model
  # d - a data frame
  # s - the name of the separating variable
  # s.at - value (typically one or zero) at which the separating variable
  #        should be set before computing the baseline chance of success.
  # s.at.lo - treat s.at as the low value when computing QIs? Defaults to
  #           TRUE. If FALSE, the s.at is treated as the high value.

  require(ggplot2)
  require(gridExtra)
  # estimate b.hat.mle
  mle <- glm(f, d, family = "binomial", x = TRUE)
  # create matrix at which to calculate the baseline
  X.pred <- apply(mle$x, 2, median)
  X.names <- names(X.pred)
  X.pred <- data.frame(matrix(X.pred, nrow = 1))
  names(X.pred) <- X.names
  X.pred <- X.pred[, -1]
  X.pred[, s] <- s.at
  # calculate the direction of the separation
  direction <- sign(coef(mle)[s])
  # choose only simulations in the direction of the separation
  prior.sims <- prior.sims[sign(prior.sims) == direction]
  # calculate the quantities of interest
  baseline <- predict(mle, newdata = X.pred)
  pr0 <- rep(plogis(baseline), length(prior.sims))
  pr1 <- plogis(baseline + prior.sims)
  if (s.at.lo == TRUE) {
    sims.fd <- pr1 - pr0
    sims.rr <- pr1/pr0
  }
  if (s.at.lo == FALSE) {
    sims.fd <- pr0 - pr1
    sims.rr <- pr0/pr1
  }
  # generate plots
  p1 <- qplot(pr1, 
        xlab = "probability") + ggtitle("PPPD of Probability")
  p2 <- qplot(pr1, log = "x",
        xlab = "probability") + ggtitle("PPPD of Probability")
  p3 <- qplot(sims.rr,
        xlab = "risk-ratio")
  p4 <- qplot(sims.rr, log = "x",
        xlab = "risk-ratio")
  p5 <- qplot(sims.fd,
        xlab = "first-difference")
  grid.arrange(p1, p2, p3, p4, p5, ncol=2)
  # define a qt function
  qt <- function(x) {
    quantile(x, 1:9/10)
  }
  q.pr <- qt(pr1)
  q.rr <- qt(sims.rr)
  q.fd <- qt(sims.fd)
  Q <- cbind(q.pr, q.rr, q.fd)
  colnames(Q) <- c("probability", "risk-ratio", "first-difference")
  print(Q, digits = 2)
  # return simulations
  ret <- list(mle, mle, 
              pr0 = pr0,
              pr1 = pr1, 
              sims.rr = sims.rr,
              sims.fd = sims.fd,
              Q = Q,
              direction = direction)
  return(ret)
}