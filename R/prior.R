state.prior.dm <- function(y, model) {

  dim_p <- model$dim_p
  nm <- names(dim_p)
  p <- sum(dim_p)
  if (is.null(model$superposition)) {

  } else {



  }


}

poly.dm.prior <- function(y, order) {
  t   <- 1:length(y)
  tmp <- lm(y ~ poly(t, degree = order - 1))
  unname(coef(tmp))
}

y = as.numeric(AirPassengers)
frequency = 12
seas.dm.prior <- function(y, frequency) {
  y <- ts(data = y, frequency = frequency)
  X <- forecast::seasonaldummy(y)
  tmp <- lm(y ~ X)
  m0 <- unname(coef(tmp))
  C0 <- diag()
}

y = as.numeric(AirPassengers)
s = 12
trig.dm.prior <- function(y, s) {
  y <- ts(data = y, frequency = s)
  X <- forecast::fourier(y, K = s/2)
  tmp <- lm(y ~ X)
  unname(coef(tmp))
}



