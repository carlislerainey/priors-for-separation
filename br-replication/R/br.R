
# load data
d <- read.csv("br-replication/data/need-rescale.csv")

d$dem_governor <- -1*d$gop_governor

# model formula
f <- oppose_expansion ~ dem_governor + obama_win + gop_leg + percent_uninsured + 
	income + percent_nonwhite + percent_metro

# simulate from posteriors
post_gelman <- sim_post_gelman(f, d, n_sims = n_sims, n_burnin = n_burnin, n_chains = n_chains, n_thin = 10)
post_jeffreys <- sim_post_jeffreys(f, d, n_sims = n_sims, n_burnin = n_burnin, n_chains = n_chains, n_thin = 10)
post <- combine_post(post_jeffreys, post_gelman)

# compare posterior locations
median(post_gelman$mcmc[, "dem_governor"])
median(post_jeffreys$mcmc[, "dem_governor"])
mean(post_gelman$mcmc[, "dem_governor"])
mean(post_jeffreys$mcmc[, "dem_governor"])
median(post_gelman$mcmc[, "dem_governor"])/median(post_jeffreys$mcmc[, "dem_governor"])
mean(post_gelman$mcmc[, "dem_governor"])/mean(post_jeffreys$mcmc[, "dem_governor"])

# compare hpd interval width
diff(HPDinterval(mcmc(post_gelman$mcmc[, "dem_governor"]), prob = .9)[1:2])/
	diff(HPDinterval(mcmc(post_jeffreys$mcmc[, "dem_governor"]), prob = .9)[1:2])

# functions
## find median and 90% hpd and add to plot
sumry <- function(sims, ht, lab = NA) {
	q <- quantile(sims, c(.05, .95))
	hpd <- HPDinterval(mcmc(sims), prob = .9)[1:2]
	lines(hpd, c(ht, ht), lwd = 2)
	#points(hpd, c(ht, ht), pch = 4)
	d <- density(sims, adjust = 1.5)
	i.max <- which.max(d$y)
	post.mode <- d$x[i.max]
	post.mean <- mean(sims)
	post.median <- median(sims)
	#points(post.mode, ht, pch = 22, bg = "white")
	#points(post.median, ht, pch = 24, bg = "grey")
	points(post.median, ht, pch = 21, bg = "black", cex = .7)
	text(post.median, ht, lab, pos = 3, cex = .7)
	#print(paste("HPD = ", round(hpd, 1)))
	#print(paste("CI = ", round(q, 1)))
	#print(paste("mode = ", round(post.mode, 1)))
	#print(paste("median = ", round(post.median, 1)))
	#print(paste("mean = ", round(post.mean, 1)))
	s <- round(c(hpd[1], post.median, hpd[2]), 1)
	return(s)
}
## plot posteior density and HPD for coefficients
plot.posterior.density <- function(sims) {
	dens <- density(sims)
	hpd <- HPDinterval(mcmc(sims), prob = .9)
	lo <- hpd[1]
	hi <- hpd[2]
	x <- dens$x[dens$x < hi & dens$x > lo]
	y <- dens$y[dens$x < hi & dens$x > lo]
	zeros <- rep(0, length(y))
	polygon(c(x, rev(x)), c(y, zeros), col = "grey80", lty = 0)
	abline(v = 0, col = "grey40")
	lines(dens, lwd = 2)
}

# coefficient plot 
pdf("doc/figs/br-coef-illustrate-importance.pdf", height = 3.5, width = 6)
par(mfrow = c(1,1), mar = c(4,1,1,1), oma = c(0,0,0,0))
eplot(xlim = c(-15, 0), ylim = c(0.55, 2.5),
			anny = FALSE,
			xlabpos = 2.5,
			xlab = "Posterior Median and 90% HPD for\nCoefficient of Democratic Governor Indicator")
abline(v = 0, col = "grey80")
sumry(post_jeffreys$mcmc[, "dem_governor"], 2, "Zorn's Default Jeffreys' Prior")
sumry(post_gelman$mcmc[, "dem_governor"], 1, "Gelman et al.'s Default Cauchy(0, 2.5) Prior")
dev.off()

# density plot
## calculate densities
d.jeffreys <- density(post_jeffreys$mcmc[, "dem_governor"])
d.gelman <- density(post_gelman$mcmc[, "dem_governor"])
pdf("doc/figs/br-posterior-density-illustrate-importance.pdf", height = 2.75, width = 8)
par(mfrow = c(1, 2), mar = c(1,1,1,1), oma = c(2,3,1,1))
eplot(xlim = c(-20, 0), 
			mm(c(d.jeffreys$y, d.gelman$y)),
			xlab = "Coefficient for Democratic Governor Indicator",
			ylab = "Posterior Density",
			ylabpos = 2.5,
			main = "Zorn's Default Jeffreys' Prior")
plot.posterior.density(post_jeffreys$mcmc[, "dem_governor"])
aplot("Gelman et al.'s Default Cauchy(0, 2.5) Prior")
plot.posterior.density(post_gelman$mcmc[, "dem_governor"])
dev.off()