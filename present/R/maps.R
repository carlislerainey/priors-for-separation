
setwd("~/Dropbox/projects/Need")

library(maps)
library(mapdata)
#library(extrafont)

d <- read.csv("~/Dropbox/projects/Need/Data/politics_and_need.csv")
dem.support <- as.character(d$state[d$gov_party == "Democrat"])
rep.support <- as.character(d$state[d$gov_position == "Supports"])
rep.oppose <- as.character(d$state[d$gov_position == "Opposes"])
rep.weighing <- as.character(d$state[d$gov_position == "Weighing Options"])


col0 <- "grey80"
col1 <- rgb(170, 86, 57, 255, maxColorValue = 255)
col2 <- rgb(39, 118, 80, 255, maxColorValue = 255)

setwd("~/Dropbox/projects/priors-for-separation/present/figs")
png("map.png", height = 4, width = 5, units = "in", res = 300)
par(mfrow = c(1, 1), oma = c(2.5,0,0,0), mar = c(0,0,0,0))
map("state", interior = FALSE, lty = 3)
map('state', 
    add = TRUE, 
    lty = 1, fill = TRUE, col = col0, 
    region = rep.weighing)
map('state', 
    add = TRUE, lty = 1, fill = TRUE, col = col2, 
    region = rep.support)
map('state', 
    add = TRUE, lty = 1, fill = TRUE, col = col1, 
    region = rep.oppose)
legend(par('usr')[1],par('usr')[3], xjust = 0, yjust = 1,
       legend = c("Supports Expansion", "Weighing Options", "Opposes Expansion"), 
       pt.bg = c(col2, col0, col1), pt.lwd = 1, pch = 22, pt.cex = 2, xpd = NA)
dev.off()

