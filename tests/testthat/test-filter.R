test_that("filterDist works", {

  library(dlmfit)
  library(magrittr)
  model <- StatePoly(order = 2, delta = 0.985) %>%
    superposition(
      StateTrig(s = 12, delta = 0.95)
    )
  my_prior <- list(
    m0 = c(100, rep(0, 12)),
    C0 = diag(c(100, 10, 10, rep(1, 10)))
  )
  m1 <- filterDist(
      model = model,
      y = AirPassengers,
      prior = my_prior
    )
  m1
  expect_equal(length(m1), 13)
})
