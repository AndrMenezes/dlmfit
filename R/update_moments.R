#' @importFrom dplyr tibble
#'
#' @name update_moments
#'
#' @title Compute the posterior and predictive moments at time \eqn{t}.
#'
#' @description The function is used in \code{fit.dm} to calculate recursively the expected value and variance of posterior
#' and predictive distributions at time \eqn{t}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @param y observed value at time \eqn{t}.
#' @param t observed time.
#' @param FF vector \eqn{\bold{FF}_t} of the model.
#' @param GG matrix \eqn{\bold{GG}_t} of the model.
#' @param a vector \eqn{\bold{a}_t} of the model.
#' @param R matrix \eqn{\bold{R} = \bold{P} + \bold{W}} of the model.
#' @param C0 prior variance at time \eqn{t} for \eqn{\bold{\theta}_t}.
#' @param n0 prior shape parameter of the Gamma distribuition at time \eqn{t} for \eqn{\phi=V^{-1}}.
#' @param d0 prior degree of freedom at time \eqn{t} for \eqn{\phi=V^{-1}}.
#' @param alpha discount factor of variance evolution.
#' @param law specify the type of variance law.
#' @param power parameter for power variance law.


update_moments <- function(y, t, FF, GG, a, R, n0, d0, alpha=1, law="identity", power=1) {

  s0 <- d0 / n0

  ## one-step ahead predictive distribution
  f <- crossprod(FF, a)
  k <- variance_law(type=law, mu=f, p=power)
  Q <- c(crossprod(FF, R %*% FF) + k * s0)

  ## posterior distribution on D_t
  e <- y - f
  A <- (R %*% FF) %*% solve(Q)

  n <- alpha * n0 + 1
  d <- c(alpha * d0 + s0 * e^2 / Q)
  s <- d / n
  m <- c(a + A %*% e)
  C <- (R - tcrossprod(A, A) * Q) * s / s0

  state <- tibble(
    t = t,
    a = m,
    dR = diag(C),
    df =  n - 1,
    parameter = paste0("theta_",1:nrow(FF))
  )
  pred <- tibble(
    t = t,
    f = f,
    q = Q,
    df = n - 1
  )

  list(
    m = m,
    C = C,
    d = d,
    n = n,
    state = state,
    pred = pred
  )
}




