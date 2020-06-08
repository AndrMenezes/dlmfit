variance_law <- function(type, mu, p) {
  switch (type,
          "identity" = 1,
          "poisson" = mu,
          "binomial" = mu * (1 - mu),
          "power" = mu ^ p
  )
}
