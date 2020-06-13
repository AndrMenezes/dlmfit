rm(list = ls())
## Packages
library(dlm)

library(dplyr)
library(magrittr)
library(zoo)

library(ggplot2)
library(cowplot)

nms <- paste0("R/",c("models", "variance_law", "update_moments", "forecast", "fit",
                     "plots", "add_pi", "smooth"), ".R")
sapply(nms, source)

##################################################################
N     <- 100
mu    <- NULL
mu[1] <- 25
Wt    <- 0.5
Vt    <- 1.0
set.seed(12345)
for (i in 2:N) {
   mu[i] <- mu[i-1] + rnorm(1, 0, sqrt(Wt))
}
yt <- rnorm(N, mu, sqrt(Vt))
plot(yt, type = "l", lwd = 2)
lines(mu, lty = 2, col = 4, lwd = 2)

my_prior <- list(
  m0 = 25,
  C0 = diag(0.5, 1)
)
m <- poly.dm(order = 1, delta = 0.85) %>%
  fit.dm(y = yt, prior = my_prior)
plot.dm.predictive(m, interval = T)
plot.dm.state(m)

s <- poly.dm(order = 1, delta = 0.85) %>%
  smooth.dm(y = yt, prior = my_prior)
plot.dm.smooth(s, distr = "mean")
plot.dm.smooth(s, which.state = "theta_1", distr = "state")


gridExtra::grid.arrange(
  plot.dm.state(m) + ggtitle("Posterior"),
  plot.dm.smooth(s, which.state = "theta_1", distr = "state") + ggtitle("Smooth")
)

