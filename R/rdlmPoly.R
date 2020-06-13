#' @importFrom dlm dlmModPoly
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#'
#'
#' @name rdlmPoly
#'
#' @title Simulate univariate DLM with polynomial structure
#'
#' @description The function simulate values (\eqn{y} and \eqn{\bold{\theta}})
#' from a univariate DLM with polynomial structure.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @param n number of observations.
#' @param order order of the polynomial model.
#' @param m0 vector for prior mean (given \eqn{D_0}).
#' @param C0 matrix for prior variance (given \eqn{D_0}).
#' @param V numeric value for the observational variance \eqn{V}.
#' @param delta discount factor \eqn{0< \delta \leq 1} to compute the evolution variance matrices \eqn{\bold{W}_t}.
#'
#'
#' @return A list with two entries: \code{y} the simulated series and \code{theta} the simulated parameters.
#'
#'
#' @details To simulate we use the hierarchical represention of DLM. For further details see the package vignette.
#'
#' @seealso \code{\link{dlmfit}}
#'
#' @examples
#'
#' set.seed(1212)
#' out = rdlmPoly(n = 100, order = 1, m0 = 100, C0 = 10, V=1.5, delta=0.9)
#' plot(out[['y']], type = 'l')
#' plot(out[['theta']][, 1], type = 'l')
#'
#' set.seed(1212)
#' out = rdlmPoly(n = 50, order = 2, m0 = c(100, 10), C0 = diag(10, 2), V=2, delta=0.8)
#' plot(out[['y']], type = 'l')
#' plot(out[['theta']][, 1], type = 'l')
#' plot(out[['theta']][, 2], type = 'l')
#'
#'
#'
#' @rdname rdlmPoly
#' @export

rdlmPoly <- function(n, order, m0, C0, V, delta=0.95)
{

  theta = matrix(nrow = n, ncol = order)
  y = numeric(n)

  m = dlmModPoly(order)
  GG = m$GG
  FF = t(m$FF)

  for(t in 1:n)
  {

    ## theta_t | D_{t-1} (priori)
    P = tcrossprod(GG %*% C0, GG)
    W = P * (1/delta - 1)
    R = P + W
    a = GG %*% m0
    theta[t, ] = MASS::mvrnorm(1, mu = a, Sigma = R * V)

    ## Y_t (preditiva)
    f = crossprod(FF, a)
    Q = as.numeric(1 + crossprod(FF,R %*% FF ))
    y[t] = rnorm(1, mean = f, sd = sqrt(Q))

    ## Atualiza parâmetros (posteriori em D_t)
    e = y[t] - f
    A = (R %*% FF) / Q
    m0 = a + A %*% e
    C0 = R - tcrossprod(A, A) * Q
  }
  list(y=y, theta=theta)
}


# n <- 50
# omega_1 <- rnorm(n, mean = 0, sd = 3)
# omega_2 <- rnorm(n, mean = 0, sd = 1)
# epsilon <- rnorm(n, mean = 0, sd = 35)
# theta_1 <- c(100, rep(NA, n-1))
# theta_2 <- c(10, rep(NA, n-1))
# y <- c()
# for (j in 2:n) {
#   theta_1[j] <- theta_1[j-1] + theta_2[j-1] + omega_1[j]
#   theta_2[j] <- theta_2[j-1] + omega_2[j]
#   y[j] <- theta_1[j] + epsilon[j]
# }
# plot(y, type = "l")
#
# model <- poly.dm(order = 2, delta = 0.95)
# fit <- fit.dm(model, y[-1],
#               prior = list(m0 = c(100, 0),
#                            C0 = diag(1000, 2)))






























