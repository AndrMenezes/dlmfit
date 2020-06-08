#' @importFrom magrittr %>%
#' @importFrom dplyr tibble rename select bind_cols bind_rows
#' @importFrom tsibble as_tsibble
#' @importFrom dlm dlmModPoly dlmModSeas dlmModTrig dlmSum bdiag
#' @importFrom tidyr pivot_wider
#' @importFrom stats is.ts is.mts
#'
#' @name dlmfit
#'
#' @title Fit univariate Dynamic Linear Models
#'
#' @description The function fits univariate Dynamic Linear Models with polynomial, season and/or harmonic
#' structures by using discount factors for \eqn{\bold{W}_t}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @param y a numeric vector or a \code{ts} class object specify the time series.
#' @param index Optionally, a vector for the time index variable.
#' @param xreg Optionally, a vector or matrix of external regressors, which must have the same number of rows as x.
#' @param order order of the polynomial model. The default corresponds to second order model.
#' @param freq number of seasons to create a DLM with seasonal factors or period to create Fourier representation of a periodic DLM component.
#' @param seasonality could be \code{none}, \code{factor} or \code{fourier}.
#' @param delta a vector indicating discount factors \eqn{0< \delta_i \leq 1} for \eqn{i=1,\ldots,k}
#' components of DLM.
#' @param level cumulative probability (\eqn{\alpha}) to compute confidence limits. The default is \code{qn = 0.05}.
#' @param prior a list to specify the prior parameters \eqn{\bold{m}_0},  \eqn{\bold{C}_0}, \eqn{n_0} and \eqn{S_0}.
#'
#'
#' @return A tidy \code{tsibble} with the columns:
#' \describe{
#'  \item{\code{index}}{time index variable}
#'  \item{\code{value}}{the observed value at time t}
#'  \item{\code{theta_k}}{the posterior mean at time \eqn{t} for the \code{k}th component of \eqn{\bold{\theta}_t}}
#'  \item{\code{dC_k}}{the posterior variance at time \eqn{t} for the \code{k}th component of \eqn{\bold{\theta}_t}}
#'  \item{\code{Lw_k}; \code{Up_k}}{the \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles, respectively, of
#'    the posterior distribution at time \eqn{t} for the \code{k}th component of \eqn{\bold{\theta}_t}}
#'  \item{\code{f}}{the mean of the one-step ahead predictive distribution at time \eqn{t}}
#'  \item{\code{Q}}{the variance of the one-step ahead predictive distribution at time \eqn{t}}
#'  \item{\code{Lw}; \code{Up}}{the  \eqn{\alpha/2} and \eqn{1-\alpha/2}
#'    quantiles, respectively, of the predictive distribution at time \eqn{t}.}
#' }
#'
#'
#' @details Inference procedure is implemented following the unified theoretical
#' of the univariate DLM. For further details see West and Harrison (1997).
#'
#' @examples
#' y <- AirPassengers
#' fit_poly <- dlmfit(y = y, order = 2, delta = 0.98)
#' fit_poly
#'
#' fit_seas1 <- dlmfit(y = y, order = 2, freq = 3, seasonality = "fourier", delta = c(0.98, 0.95))
#' fit_seas1
#'
#' fit_seas2 <- dlmfit(y = y, order = 2, freq = 12, seasonality = "factor", delta = c(0.90, 0.95))
#' fit_seas2
#'
#'
#' @rdname dlmfit
#' @export
#'
dlmfit <- function(y, index=NULL, xreg=NULL, order, freq, seasonality="none", delta, level=0.05, prior=NULL) {

  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  nt <- length(y)
  if (is.ts(y)) {
    tsb <- tsibble::as_tsibble(y)
    y   <- c(y)
  } else if (!is.ts(y)) {
    if (is.null(index)) {
      tsb <- tsibble::as_tsibble(dplyr::tibble(index=1:nt, value=y), index=index)
    } else {
      tsb <- tsibble::as_tsibble(dplyr::tibble(index=index, value=y), index=index)
    }
  }
  if (!is.null(xreg)) {
    if (!is.numeric(xreg)) {
      stop("xreg must be numeric")
    } else {
      tsb_x <- dplyr::as_tibble(xreg)
      xreg  <- as.matrix(xreg)
    }
    # tsb <- dplyr::bind_cols(tsb, tsb_x)
  }

  ## Modelo polinomial
  if (!missing(order) & seasonality == "none") {
    m0 <- rep(0, order)
    C0 <- diag(1e4, order)
    mod <- dlm::dlmModPoly(order = order)[c("FF", "GG")]
  }
  ## Modelo de sazonalidade via fatores
  else if (missing(order) & seasonality == "factor") {
    m0 <- rep(0, freq - 1)
    C0 <- diag(1e4, freq - 1)
    mod <- dlm::dlmModSeas(frequency = freq)[c("FF", "GG")]
  }
  ## Modelo de sazonalidade via fourier
  else if (missing(order) & seasonality == "fourier") {
    m0 <- rep(0, freq - 1)
    C0 <- diag(1e4, freq - 1)
    mod <- dlm::dlmModTrig(s = freq)[c("FF", "GG")]
  }
  ## Modelo polinomial + sazonalidade (via fatores)
  else if (!missing(order) & seasonality == "factor") {
    m0 <- rep(0, order + freq - 1 )
    C0 <- diag(1e4, order + freq - 1)
    mod <- dlm::dlmSum(dlm::dlmModPoly(order = order),
                       dlm::dlmModSeas(frequency = freq))[c("FF", "GG")]
    mod$super <- TRUE
  }
  ## Modelo polinomial + sazonalidade (via dourier)
  else if (!missing(order) & seasonality == "fourier") {
    m0 <- rep(0, order + freq - 1 )
    C0 <- diag(1e4, order + freq - 1)
    mod <- dlm::dlmSum(dlm::dlmModPoly(order = order),
                       dlm::dlmModTrig(s = freq))[c("FF", "GG")]
    mod$super <- TRUE
  }
  ## Modelo de regressão
  else if (!is.null(xreg)) {
    k   <- ncol(xreg)
    mod <- list(GG=diag(k+1), FF=matrix(c(1), nrow = 1))
    # ## Incluí componente sazonal se for especificada
    if (seasonality == "factor") {
      seas <- dlm::dlmModSeas(frequency = freq)[c("FF", "GG")]
      mod[["GG"]] <- dlm::bdiag(mod[["GG"]][, 1], seas[["GG"]], mod[["GG"]][, -1])
      mod[["FF"]] <- cbind(mod[["FF"]], seas[["FF"]])
    } else if (seasonality == "fourier") {
      seas <- dlm::dlmModTrig(s = freq)[c("FF", "GG")]
      mod[["GG"]] <- dlm::bdiag(mod[["GG"]][, 1], seas[["GG"]], mod[["GG"]][, -1])
      mod[["FF"]] <- cbind(mod[["FF"]], seas[["FF"]])
    }
    p  <- nrow(mod[["GG"]])
    m0 <- rep(0, p)
    C0 <- diag(1e4, p)
  }

  n0 <- 1
  S0 <- 10
  if(!is.null(prior)) {
    nm <- names(prior)
    if (!(nm %in% c("n0", "S0", "m0", "C0"))) {
      stop("the prior names is incorrectly specified")
    }
    else {
      if (nm == "n0") n0 <- prior[["n0"]]
      if (nm == "S0") S0 <- prior[["S0"]]
      if (nm == "m0") m0 <- prior[["m0"]]
      if (nm == "C0") C0 <- prior[["C0"]]
    }
  }

  aux <- list()

  if (is.null(xreg)) {
    GG  <- mod$GG
    FF  <- t(mod$FF)
    if (is.null(mod$super)) {
      for (t in 1:nt) {
        P <- tcrossprod(GG %*% C0, GG)
        W <- P * (1 / delta - 1)
        R <- P + W

        out <- update_moments(y[t], t, FF, GG, R, m0, C0, S0, n0, level)

        m0 <- out[["m"]]
        C0 <- out[["C"]]
        S0 <- out[["S"]]
        n0 <- out[["n"]]

        aux[[t]] <- out[["tb"]]
      }
    }
    else {
      if(length(delta) < 2) delta <- rep(delta, 2)

      G1 <- GG[1:order, 1:order]
      G2 <- GG[-(1:order), -(1:order)]
      FF <- as.matrix(FF[,1] + FF[,2])

      for (t in 1:nt) {
        ## Priori em t ( theta_t | D_{t-1} )
        P1 <- tcrossprod(G1 %*% C0[1:order, 1:order], G1)
        W1 <- P1 * (1 / delta[1] - 1)
        P2 <- tcrossprod(G2 %*% C0[-(1:order), -(1:order)], G2)
        W2 <- P2 * (1 / delta[2] - 1)
        P  <- dlm::bdiag(P1, P2)
        W  <- dlm::bdiag(W1, W2)
        R  <- P + W

        out <- update_moments(y[t], t, FF, GG, R, m0, C0, S0, n0, level)

        m0 <- out[["m"]]
        C0 <- out[["C"]]
        S0 <- out[["S"]]
        n0 <- out[["n"]]

        aux[[t]] <- out[["tb"]]
      }
    }
  } else if (!is.null(xreg))  {
    if (seasonality == "none") {
      if(length(delta) < 2) delta <- rep(delta, 2)
      GG <- mod$GG

      G1 <- GG[1, 1]
      G2 <- GG[-1, -1]

      for (t in 1:nt) {
        P1 <- tcrossprod(G1 %*% C0[1, 1], G1)
        W1 <- P1 * (1 / delta[1] - 1)
        P2 <- tcrossprod(G2 %*% C0[-1, -1], G2)
        W2 <- P2 * (1 / delta[2] - 1)
        P  <- dlm::bdiag(P1, P2)
        W  <- dlm::bdiag(W1, W2)
        R  <- P + W
        FF <- rbind(1, xreg[t,])

        out <- update_moments(y[t], t, FF, GG, R, m0, C0, S0, n0, level)

        m0 <- out[["m"]]
        C0 <- out[["C"]]
        S0 <- out[["S"]]
        n0 <- out[["n"]]

        aux[[t]] <- out[["tb"]]
      }
    }
    else if (seasonality != "none") {
     stop("Regression with seasonality is not available yet")
    }
  }

  pred <- dplyr::bind_rows(aux) %>%
    tidyr::pivot_wider(names_from = parms, values_from = value)

  dplyr::bind_cols(tsb, dplyr::select(pred, -c(t, y)))

}
