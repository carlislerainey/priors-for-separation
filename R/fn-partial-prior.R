
ppd <- function(f, d, prior.sims, s, s.at, s.at.lo = TRUE) {
  # f - a logitistic regression model
  # d - a data frame
  # s - the name of the separating variable
  # s.at - value (typically one or zero) at which the separating variable
  #        should be set before computing the baseline chance of success.
  # s.at.lo - treat s.at as the low value when computing QIs? Defaults to
  #           TRUE. If FALSE, the s.at is treated as the high value.
  require(ggplot2)
  require(gridExtra)
  # set temporary options
  options.restore <- options()
  options("scipen" = 12)
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
    fd <- pr1 - pr0
    rr <- pr1/pr0
  }
  if (s.at.lo == FALSE) {
    fd <- pr0 - pr1
    rr <- pr0/pr1
  }
  # check for infinite risk ratios.
  if (max(rr) == Inf) {
    cat(paste("\n\n####### WARNING #######\n\n",
              "Of the ", 
              length(prior.sims), 
              " prior simulations, ",
              sum(rr == Inf), " (", round(100*sum(rr == Inf)/length(prior.sims), 2), "%)",
              " produced risk-ratios of infinity.\nThe largest ", 
              ceiling(300*sum(rr == Inf)/length(prior.sims)), "% of simulations are being dropped.\n\n", sep = ""))
    keep <- rr <= quantile(rr, 1 - ceiling(300*sum(rr == Inf)/length(prior.sims))/100)
    pr0 <- pr0[keep]
    pr1 <- pr1[keep]
    rr <- rr[keep]
    fd <- fd[keep]
  }
  # generate plots
  plots <- list()
  print(range(pr1))
  plots$pr1 <- qplot(pr1,
              xlab = "probability") + ggtitle("PPPD of Probability")
  print(range(log(pr1)))
  plots$log.pr1 <- qplot(pr1, log = "x",
              xlab = "probability") + ggtitle("PPPD of Probability (Log Scale)")
  print(range(rr))
  plots$rr <- qplot(rr,
              xlab = "risk-ratio") + ggtitle("PPPD of Risk-Ratio")
  print(range(log(rr)))
  plots$log.rr <- qplot(rr, log = "x",
              xlab = "risk-ratio") + ggtitle("PPPD of Risk-Ratio (Log Scale)")
  print(range(fd))
  plots$fd <- qplot(fd,
              xlab = "first-difference") + ggtitle("PPPD of First-Difference")
  grid.arrange(plots$pr1, 
               plots$log.pr1, 
               plots$rr,
               plots$log.rr,
               plots$fd, ncol=2)
  # define a qt function to compute the deciles
  qt <- function(x) {
    quantile(x, 1:9/10)
  }
  q.pr <- qt(pr1)
  q.rr <- qt(rr)
  q.fd <- qt(fd)
  Q <- cbind(q.pr, q.rr, q.fd)
  colnames(Q) <- c("probability", "risk-ratio", "first-difference")
  print(Q, digits = 2)
  # return simulations
  ret <- list(mle, mle, 
              pr0 = pr0,
              pr1 = pr1, 
              rr = rr,
              fd = fd,
              Q = Q,
              direction = direction,
              plots = plots)
  return(ret)
  # restore options
  options(options.restore)
}
