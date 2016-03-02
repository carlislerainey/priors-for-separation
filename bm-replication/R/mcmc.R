

# load data set
d <- read.csv("bm-replication/data/bm.csv")
d <- d[, c("warl2", "onenukedyad", "twonukedyad", "logCapabilityRatio", "Ally",
           "SmlDemocracy", "SmlDependence", "logDistance", "Contiguity",
           "MajorPower", "NIGOs")]
d <- na.omit(d)

# center data
d$c.onenukedyad <- rescale(d$onenukedyad)
d$c.twonukedyad <- rescale(d$twonukedyad)
d$c.logCapabilityRatio <- rescale(d$logCapabilityRatio)
d$c.Ally <- rescale(d$Ally)
d$c.SmlDemocracy <- rescale(d$SmlDemocracy)
d$c.SmlDependence <- rescale(d$SmlDependence)
d$c.logDistance <- rescale(d$logDistance)
d$c.Contiguity<- rescale(d$Contiguity)
d$c.MajorPower <- rescale(d$MajorPower)
d$c.NIGOs <- rescale(d$NIGOs)

# set formula
f <- warl2 ~ c.onenukedyad + c.twonukedyad + c.logCapabilityRatio + 
  c.Ally + c.SmlDemocracy + c.SmlDependence + c.logDistance + 
  c.Contiguity + c.MajorPower + c.NIGOs

# my prior
m.me <- sim_post_normal(f, d, sep_var = "c.twonukedyad", sd = 4.5, 
               n_sims = 100, n_burnin = 50)
save(m.me, file = "output/m-me.RData")

# skeptical prior
m.skep <- sim_post_normal(f, d, sep_var = "c.twonukedyad", sd = 2, 
                 n_sims = 100, n_burnin = 50)
save(m.skep, file = "output/m-skep.RData")

# enthusiastic prior
m.enth <- sim_post_normal(f, d, sep_var = "c.twonukedyad", sd = 8, 
                 n_sims = 100, n_burnin = 50)
save(m.enth, file = "output/m-enth.RData")

# Zorn's default
m.zorn <- sim_post_jeffreys(f, d, n_sims = 100, n_burnin = 50, n_thin = 1, tune = .1)
save(m.zorn, file = "output/m-zorn.RData")

#G elman's default
m.gelman <- sim_post_gelman(f, d, n_sims = 100, n_burnin = 50, n_thin = 1)
save(m.gelman, file = "output/m-gelman.RData")

