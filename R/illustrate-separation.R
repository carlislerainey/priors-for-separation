# set working director
setwd("~/Dropbox/projects/priors-for-separation")

# load packages
library(compactr)
library(texreg)

# set seed
set.seed(8742570)

pdf("doc/figs/illustrate-separation.pdf", height = 1.5, width = 6)
par(mfrow = c(1, 3), 
    mar = rep(.5, 4),
    oma = c(2,3,.5,.5))
col0 <- rgb(0,0,0,.25)

# overlap
n <- 100
x <- runif(n)
y <- rbinom(n, 1, plogis(-10 + 20*x))
m.overlap <- glm(y ~ x, family = "binomial")
m.overlap.precise <- glm(y ~ x, family = "binomial", control = list(epsilon = 1e-16))
eplot(xlim = c(-.025, 1.025), ylim = c(-.025, 1.025),
      xlab = "x", ylab = "y", ylabpos = 2.3,
      main = "Overlap")
points(x, y, col = col0)
curve(plogis(coef(m.overlap)[1] + coef(m.overlap)[2]*x), add = TRUE, lwd = 1)

# complete separation
y <- 1*(x > .5)
m.complete <- glm(y ~ x, family = "binomial")
m.complete.precise <- glm(y ~ x, family = "binomial", control = list(epsilon = 1e-16))

aplot("Complete Separation")
points(x, y, col = col0)
curve(plogis(coef(m.complete)[1] + coef(m.complete)[2]*x), add = TRUE, lwd = 1)

# quasicomplete separation
x <- round(x)
y[x == 0] <- 0
y[x == 1] <- rbinom(sum(x), 1, .5)
m.quasi <- glm(y ~ x, family = "binomial")
m.quasi.precise <- glm(y ~ x, family = "binomial", control = list(epsilon = 1e-16))
aplot("Quasicomplete Separation")
points(jitter(x, factor = .2), jitter(y, factor = .2), col = col0)
curve(plogis(coef(m.quasi)[1] + coef(m.quasi)[2]*x), add = TRUE, lwd = 1)
dev.off()

# create table
texreg(list(m.overlap, m.overlap.precise, 
               m.complete, m.complete.precise,
               m.quasi, m.quasi.precise))