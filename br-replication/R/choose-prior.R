
n.sims <- 1000
b0 <- rnorm(n.sims, 0, 3) #rep(-5, n.sims) #rcauchy(n, 0, 1)
b1 <- rnorm(n.sims, 0, 3)

p.lo <- plogis(b0 - .5*b1)
p.hi <- plogis(b0 + .5*b1)
dif <- p.hi - p.lo

par(mfrow = c(1,3))
hist(p.lo); hist(p.hi); hist(dif)

d.p.lo <- density(p.lo)
d.p.hi <- density(p.hi)
d.dif <- density(dif)




par(mfrow = c(1,3))
hist(p.hi - p.lo)
plot(p.lo, p.hi, col = rgb(1*(p.lo > p.hi), 0,0, .3))
eplot(xlim = c(0, 1), ylim = c(0, 1))
for (i in 1:n.sims) {
  if (p.lo[i] > p.hi[i] & p.lo[i] > 0.5) {
  lines(c(0, 1), c(p.lo[i], p.hi[i]), col = rgb(1*(p.lo[i] > p.hi[i]), 0,0, .1))
  }
}
