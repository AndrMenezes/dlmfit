#' @title Forecast function for Dynamic Models
#'
#' @description Forecast univariate Dynamic Models.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @param object \code{dm} class object from \code{fit.dm} function.
#' @param horizon horizon of the forecast
#' @param distr forecast the predictive, state or both distributions?
#'
#' @importFrom utils tail
#' @importFrom zoo index
#'
#' @rdname forecast
#' @export

forecast <- function(object, horizon, distr=c("state", "predictive", "both")) {

  ## specifying the changes in time
  t  <- tail(index(object$y), 2)
  dt <- diff(tail(t, 2))
  t  <- max(t)

  FF <- t(object$model$FF)
  GG <- object$model$GG

  ## Static values
  Wt <- object$final_state$Wt
  st <- object$final_state$st
  nt <- object$final_state$nt

  ## initial values
  mt <- object$final_state$mt
  Rt_h <- object$final_state$Ct

  fore_state <- tibble()
  fore_pred <- tibble()

  GG_h <- GG
  for (h in seq_len(horizon)) {

    ## state distribution
    at_h <- GG_h %*% mt
    Rt_h <- tcrossprod(GG_h %*% Rt_h, GG_h) + Wt
    ## predictive distribution
    ft_h <- c(crossprod(FF, at_h))
    qt_h <- c(crossprod(FF, Rt_h %*% FF) + st)

    GG_h <- GG_h %*% GG

    fore_state <- rbind(
      fore_state,
      tibble(
        t = t + h * dt,
        a = c(at_h),
        dR = diag(Rt_h),
        df = nt,
        parameter = paste0("theta_",1:nrow(FF))
      )
    )

    fore_pred <- rbind(
      fore_pred,
      tibble(
        t = t + h * dt,
        f = ft_h,
        q = qt_h,
        df = nt
      )
    )
  }

  out <- switch (distr,
    "state" = fore_state,
    "predictive" = fore_pred,
    "both" = list(
      fore_state = fore_state,
      fore_pred = fore_pred
    )
  )
  out
}
