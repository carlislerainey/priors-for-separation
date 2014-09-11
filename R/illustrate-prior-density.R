# clear workspace
rm(list = ls())

# set working directory
setwd("~/Dropbox/projects/priors-for-separation/")

# load packages
library(compactr)

x <- seq(-10, 10, length.out = 1000)
d.norm <- dnorm(x, sd = 2.5)
d.cauchy <- dcauchy(x, scale = 2.5)

pdf("doc/figs/illustrate-prior-density.pdf", height = 3, width = 4)
par(mfrow = c(1,1), mar = c(3,4,1,1), oma = c(0,0,0,0))
eplot(xlim = mm(x), ylim = 1.1*mm(c(d.norm, d.cauchy)),
      xlab = expression(beta),
      ylab = expression(p(beta)),
      ylabpos = 2.3)
lines(x, d.norm, col = "black", lwd = 2)
lines(x, d.cauchy, col = "red", lwd = 2, lty = 10)
text(0, dnorm(0, sd = 2.5), "Normal(2.5)", 
     pos = 3, cex = .8)
lines(c(0, 2), y = c(0.132, 0.14), col = "red")
text(2, 0.14, "Cauchy(2.5)", 
     pos = 4, cex = .8, col = "red")
dev.off()