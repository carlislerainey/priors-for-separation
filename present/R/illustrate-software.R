

# # install packages (if necessary)
# install.packages("devtools")
# devtools::install_github("carlislerainey/compactr")
# devtools::install_github("carlislerainey/separation")
# install.packages("arm")
# install.packages("ggplot")
# devtools::install_github("jrnold/ggthemes")

# load packages
library(separation)
library(arm)  # for rescale()

# load and recode data
data(politics_and_need) 
d <- politics_and_need
d$dem_governor <- 1 - d$gop_governor
d$st_percent_uninsured <- rescale(d$percent_uninsured)

# formula to use throughout
f <- oppose_expansion ~ dem_governor + percent_favorable_aca + gop_leg +
  st_percent_uninsured + bal2012 + multiplier + percent_nonwhite + percent_metro

# informative prior
prior_sims_4.5 <- rnorm(10000, 0, 4.5)
pppd <- calc_pppd(formula = f, 
                  data = d, 
                  prior_sims = prior_sims_4.5, 
                  sep_var_name = "dem_governor",
                  prior_label = "Normal(0, 4.5)")
print(pppd)
setwd("~/Dropbox/projects/priors-for-separation/present/figs")
png("ill-plot-pppd.png", height = 4, width = 6, units = "in", res = 300)
plot(pppd)
dev.off()

png("ill-plot-pppd-log.png", height = 4, width = 6, units = "in", res = 300)
plot(pppd, log_scale = TRUE)
dev.off()

# mcmc estimation
post <- sim_post_normal(f, d, sep_var = "dem_governor",
                        sd = 4.5, 
                        n_sims = 10000,
                        n_burnin = 1000,
                        n_chains = 4)
print(post)

# compute quantities of interest
## dem_governor
X_pred_list <- set_at_median(f, d)
x <- c(0, 1)
X_pred_list$dem_governor <- x
qi <- calc_qi(post, X_pred_list, qi_name = "fd")

png("ill-plot-qi-dem-gov.png", height = 4, width = 6, units = "in", res = 300)
plot(qi, xlim = c(-1, 1),
     xlab = "First Difference",
     ylab = "Posterior Density",
     main = "The Effect of Democratic Partisanship on Opposing the Expansion")
dev.off()
## st_percent_uninsured
X_pred_list <- set_at_median(f, d)
x <- seq(min(d$st_percent_uninsured),
         max(d$st_percent_uninsured),
         by = 0.1)
X_pred_list$st_percent_uninsured <- x
qi <- calc_qi(post, X_pred_list, qi_name = "pr")

png("ill-plot-qi-perc-unins.png", height = 4, width = 6, units = "in", res = 300)
plot(qi, x,
     xlab = "Percent Uninsured (Std.)",
     ylab = "Predicted Probability",
     main = "The Probability of Opposition as the Percent Uninsured (Std.) Varies")
dev.off()

# abandon separation for better things

## pull of coef. simulation for convenience
x <- post$mcmc[, "dem_governor"]
y <- post$mcmc[, "st_percent_uninsured"]

## base graphics
png("ill-abandon-base-graphics.png", height = 4, width = 6, units = "in", res = 300)
plot(x, y, 
     xlab = "Coefficient for Democratic Governor",
     ylab = "Coefficient for Percent Uninsured (Std.)",
     main = "Comparing the Effect of Partisanship and Need",
     col = rgb(0, 0, 0, .5))
abline(a = 0, b = 1, col = "red", lwd = 3)
dev.off()

## ggplot

### load a couple more packages
library(ggplot2)
library(ggthemes)

### default
q <- qplot(x, y, 
      xlab = "Coefficient for Democratic Governor",
      ylab = "Coefficient for Percent Uninsured (Std.)",
      main = "Comparing the Effect of Partisanship and Need") + 
  geom_abline(); 
png("ill-abandon-gg-default.png", height = 4, width = 6, units = "in", res = 300)
q
dev.off()

### Stata-like
png("ill-abandon-gg-stata.png", height = 4, width = 6, units = "in", res = 300)
q + theme_stata()
dev.off()

### Excel-like
png("ill-abandon-gg-excel.png", height = 4, width = 6, units = "in", res = 300)
q + theme_excel()
dev.off()



# mcmc estimation
n_sims <- 10000
n_burnin <- 1000
post_inf <- sim_post_normal(f, d, sep_var = "dem_governor",
                        sd = 4.5, 
                        n_sims = n_sims,
                        n_burnin = n_burnin,
                        n_chains = 8)
post_skep <- sim_post_normal(f, d, sep_var = "dem_governor",
                        sd = 2.25, 
                        n_sims = n_sims,
                        n_burnin = n_burnin,
                        n_chains = 4)
post_enth <- sim_post_normal(f, d, sep_var = "dem_governor",
                          sd = 9, 
                          n_sims = n_sims,
                          n_burnin = n_burnin,
                          n_chains = 4)
post_jeffreys <- sim_post_jeffreys(f, d,
                          n_sims = n_sims,
                          n_burnin = n_burnin,
                          n_chains = 4)
data(politics_and_need_rescaled)
politics_and_need_rescaled$dem_governor <- - politics_and_need_rescaled$gop_governor
post_gelman <- sim_post_jeffreys(oppose_expansion ~ dem_governor + percent_favorable_aca + gop_leg +
                                   percent_uninsured + bal2012 + multiplier + percent_nonwhite + percent_metro, 
                                 politics_and_need_rescaled,
                                   n_sims = n_sims,
                                   n_burnin = n_burnin,
                                   n_chains = 4)


col1a <- rgb(170, 86, 57, 150, maxColorValue = 255)
col2a <- rgb(39, 118, 80, 150, maxColorValue = 255)
col1 <- rgb(170, 86, 57, 255, maxColorValue = 255)
col2 <- rgb(39, 118, 80, 255, maxColorValue = 255)

library(compactr)
png("ill-conclusion.png", height = 4, width = 6.5, units = "in", res = 300)
par(mfrow = c(2, 3), oma = c(3, 3, 1, 1), mar = c(1, 1, 1, 1))
# skeptical prior
eplot(xlim = range(post_skep$mcmc[, "dem_governor"],
                   post_inf$mcmc[, "dem_governor"],
                   post_enth$mcmc[, "dem_governor"]),
      ylim = range(post_skep$mcmc[, "st_percent_uninsured"],
                   post_inf$mcmc[, "st_percent_uninsured"],
                   post_enth$mcmc[, "st_percent_uninsured"]),
      xlab = "Coefficient for Democratic Governor",
      ylab = "Coefficient for Percent Uninsured",
      main = "Normal(0, 2.25)")
post_skep$mcmc <- post_skep$mcmc[sample(1:n_sims), ]
cwh <- 1*(post_skep$mcmc[, "dem_governor"] < post_skep$mcmc[, "st_percent_uninsured"])
colr <- colg <- colb <- numeric(n_sims)
colr[cwh == 1] <- col2rgb(col2)[1]
colg[cwh == 1] <- col2rgb(col2)[2]
colb[cwh == 1] <- col2rgb(col2)[3]
colr[cwh == 0] <- col2rgb(col1)[1]
colg[cwh == 0] <- col2rgb(col1)[2]
colb[cwh == 0] <- col2rgb(col1)[3]

abline(a = 0, b = 1)
points(post_skep$mcmc[, "dem_governor"], 
       post_skep$mcmc[, "st_percent_uninsured"],
       col = rgb(colr, colg, colb, 50, maxColorValue = 250))
text(-12, 10, paste("Pr(Research Hyp.) = ", round(mean(cwh), 2), sep = ""), cex = 0.7)
# informative prior
aplot("Normal(0, 4.5)")
post_inf$mcmc <- post_inf$mcmc[sample(1:n_sims), ]
cwh <- 1*(post_inf$mcmc[, "dem_governor"] < post_inf$mcmc[, "st_percent_uninsured"])
colr <- colg <- colb <- numeric(n_sims)
colr[cwh == 1] <- col2rgb(col2)[1]
colg[cwh == 1] <- col2rgb(col2)[2]
colb[cwh == 1] <- col2rgb(col2)[3]
colr[cwh == 0] <- col2rgb(col1)[1]
colg[cwh == 0] <- col2rgb(col1)[2]
colb[cwh == 0] <- col2rgb(col1)[3]
abline(a = 0, b = 1)
points(post_inf$mcmc[, "dem_governor"], 
       post_inf$mcmc[, "st_percent_uninsured"],
       col = rgb(colr, colg, colb, 50, maxColorValue = 250))
text(-12, 10, paste("Pr(Research Hyp.) = ", round(mean(cwh), 2), sep = ""), cex = 0.7)
# enthusiastic prior
aplot("Normal(0, 9)")
post_enth$mcmc <- post_enth$mcmc[sample(1:n_sims), ]
cwh <- 1*(post_enth$mcmc[, "dem_governor"] < post_enth$mcmc[, "st_percent_uninsured"])
colr <- colg <- colb <- numeric(n_sims)
colr[cwh == 1] <- col2rgb(col2)[1]
colg[cwh == 1] <- col2rgb(col2)[2]
colb[cwh == 1] <- col2rgb(col2)[3]
colr[cwh == 0] <- col2rgb(col1)[1]
colg[cwh == 0] <- col2rgb(col1)[2]
colb[cwh == 0] <- col2rgb(col1)[3]
abline(a = 0, b = 1)
points(post_enth$mcmc[, "dem_governor"], 
       post_enth$mcmc[, "st_percent_uninsured"],
       col = rgb(colr, colg, colb, 50, maxColorValue = 250))
text(-12, 10, paste("Pr(Research Hyp.) = ", round(mean(cwh), 2), sep = ""), cex = 0.7)
# Jeffreys prior
aplot("Jeffreys' Prior")
post_jeffreys$mcmc <- post_jeffreys$mcmc[sample(1:n_sims), ]
cwh <- 1*(post_jeffreys$mcmc[, "dem_governor"] < post_jeffreys$mcmc[, "st_percent_uninsured"])
colr <- colg <- colb <- numeric(n_sims)
colr[cwh == 1] <- col2rgb(col2)[1]
colg[cwh == 1] <- col2rgb(col2)[2]
colb[cwh == 1] <- col2rgb(col2)[3]
colr[cwh == 0] <- col2rgb(col1)[1]
colg[cwh == 0] <- col2rgb(col1)[2]
colb[cwh == 0] <- col2rgb(col1)[3]
abline(a = 0, b = 1)
points(post_jeffreys$mcmc[, "dem_governor"], 
       post_jeffreys$mcmc[, "st_percent_uninsured"],
       col = rgb(colr, colg, colb, 50, maxColorValue = 250))
text(-12, 10, paste("Pr(Research Hyp.) = ", round(mean(cwh), 2), sep = ""), cex = 0.7)
# Gelman et al.'s prior
aplot("Gelman et al.'s Prior")
post_gelman$mcmc <- post_gelman$mcmc[sample(1:n_sims), ]
cwh <- 1*(post_gelman$mcmc[, "dem_governor"] < post_gelman$mcmc[, "percent_uninsured"])
colr <- colg <- colb <- numeric(n_sims)
colr[cwh == 1] <- col2rgb(col2)[1]
colg[cwh == 1] <- col2rgb(col2)[2]
colb[cwh == 1] <- col2rgb(col2)[3]
colr[cwh == 0] <- col2rgb(col1)[1]
colg[cwh == 0] <- col2rgb(col1)[2]
colb[cwh == 0] <- col2rgb(col1)[3]
abline(a = 0, b = 1)
points(post_gelman$mcmc[, "dem_governor"], 
       post_gelman$mcmc[, "percent_uninsured"],
       col = rgb(colr, colg, colb, 150, maxColorValue = 250))
text(-12, 10, paste("Pr(Research Hyp.) = ", round(mean(cwh), 2), sep = ""), cex = 0.7)
dev.off()
