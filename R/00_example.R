rm(list = ls())
## Packages
library(dlm)

library(dplyr)
library(magrittr)
library(zoo)

library(ggplot2)
library(cowplot)

nms <- paste0("R/",c("models", "variance_law", "update_moments", "forecast", "fit",
                "plots", "add_pi"), ".R")
sapply(nms, source)

## Simple model
y <- log(AirPassengers)

model <- poly.dm(order = 2, delta = 0.98) %>%
  superposition.dm(
    seas.dm(frequency = 12, delta=0.95)
    #trig.dm(s = 12, delta = 0.95)
  )
m <- model %>%
  fit.dm(y)
x11()
plot.dm.predictive(m, interval = T)

my_prior <- list(
  m0 = c(log(100), rep(0, 12)),
  C0 = diag(c(50, 10, 10, rep(1, 10)))
)
m1 <- model %>%
  fit.dm(y, prior = my_prior)

x11()
plot.dm.predictive(m1, interval = T)

x11()
(pf_1 <- plot.dm.forecast(m1, horizon = 12, distr = "predictive"))
fit_ets <- ets(y, model = "ANA")
x11();plot(forecast(fit_ets, h = 12))

## Variance law model (poison)
m2 <- model %>%
  fit.dm(y, var_law = "poisson")
m2$lpl

plot.dm.predictive(m2)
plot.dm.state(m2, which = "theta_1")
(pf_2 <- plot.dm.forecast(m2, horizon = 12, distr = "predictive"))

## Variance law and evolution

m3 <- model %>%
  fit.dm(y, var_law = "identity", delta_phi = 0.96)
m3$lpl

plot.dm.predictive(m3)
plot.dm.state(m3, which = "theta_1")
(pf_3 <- plot.dm.forecast(m3, horizon = 12, distr = "predictive"))

## Comparison in terms of LPL
tibble(
  model = paste0("m_", 1:3),
  lpl = c(m1$lpl, m2$lpl, m3$lpl)
)

## Multiple forecast plots
plot_grid(
  pf_1+ggtitle("m1"),
  pf_2+ggtitle("m2"),
  pf_3+ggtitle("m3")
)

library(forecast)
fit_ets <- ets(y, model = "MAM")
x11();plot(forecast(fit_ets, h = 12))
x11();pf_3
m1 %>%
  forecast.dm(horizon = 12, which = "predictive")


## Example West and Harisson (1993) pag. 318
y  <- c(8.48, 8.7, 8.09, 8.58, 8.94, 8.86, 8.45, 9, 9.2, 9.11, 8.69, 8.87, 9.13,
        9.23, 8.65, 8.84, 9.23, 9.21, 8.68, 9.2, 9.49, 9.54, 9.06, 9.35, 9.37,
        9.66, 9.03, 9.44, 9.56, 9.98, 9.19, 9.5, 9.71, 9.6, 9.18, 9.53, 9.72,
        9.88, 9.11, 9.49, 9.82, 9.9, 8.87, 9.38, 10.11, 9.9, 9.47, 9.47)
y <- ts(y, start = 73, frequency = 4)
plot(y, type = "o", pch = 19)

model <- poly.dm(order = 2, delta = 0.85) %>%
  superposition.dm(
    trig.dm(s = 4, delta = 0.97)
  )
m1 <- model %>%
  fit.dm(
    y = y,
    prior = list(
      m0 = c(8, rep(0, 4)),
      C0 = diag(c(100, 10, 10, 1, 1))
    )
  )

x11()
plot.dm.predictive(m1)
plot.dm.state(m1, which = "theta_2", interval.level = 0.1)
x11()
plot.dm.forecast(m1, horizon = 12, distr = "predictive")

fit_ets <- ets(y, model = "MAM")
x11();plot(forecast(fit_ets, h = 12))
