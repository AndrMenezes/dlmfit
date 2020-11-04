#' @name kalmanDist
#'
#' @title Fitting Dynamic Models
#'
#' @description Fitting univariate Dynamic Models via discount factor principle.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @param model \code{dm} class object.
#' @param y observed data.
#' @param prior list specifying the prior parameters.
#' @param var_law name of the variance law.
#' @param pow parameter power for the variance law.
#' @param delta_phi discount factor for variance evolution.
#'
#'
#' @importFrom zoo as.zoo
#' @importFrom dplyr bind_rows mutate summarise pull
#' @importFrom magrittr %>%
#' @importFrom stats dt sd quantile
NULL

#' @rdname kalmanDist
#' @export

filterDist <- function(model, y, prior=NULL, var_law="identity", pow=1, delta_phi=1) {

  FF <- t(model$FF)
  GG <- model$GG
  D  <- model$D

  p  <- nrow(FF)
  ## Define prior distribution
  m0 <- rep(0, p)
  C0 <- diag(100, p)
  n0 <- 0.01
  d0 <- 0.01
  if(!is.null(prior)) {
    nm <- names(prior)
    if (!(any(nm %in% c("n0", "d0", "m0", "C0")))) {
      stop("the prior names is incorrectly specified")
    }
    else {
      if (any(nm %in% "n0")) n0 <- prior[["n0"]]
      if (any(nm %in% "d0")) d0 <- prior[["d0"]]
      if (any(nm %in% "m0")) m0 <- prior[["m0"]]
      if (any(nm %in% "C0")) C0 <- prior[["C0"]]
    }
  }

  prior <- list(
    m0 = m0,
    C0 = C0,
    n0 = n0,
    d0 = d0
  )

  ## response variable
  y_zoo <- as.zoo(y)
  y <- as.numeric(y)
  nt <- length(y)

  ## empty tibbles to save results
  tb_state <- tibble()
  tb_pred <- tibble()

  for (t in 1:nt) {
    ## Prior (theta_t | D_{t-1})
    P <- tcrossprod(GG %*% C0, GG)
    R <- D * P
    a <- GG %*% m0

    out <- update_moments(y[t], t,
                          FF, GG, a, R,
                          n0, d0)
    m0 <- out[["m"]]
    C0 <- out[["C"]]
    d0 <- out[["d"]]
    n0 <- out[["n"]]

    tb_state <- bind_rows(tb_state, out[["state"]])
    tb_pred <- bind_rows(tb_pred, out[["pred"]])
  }


  final_state <- list(
    Wt = R - P,
    mt = m0,
    Ct = C0,
    st = d0 / n0,
    nt = n0
  )

  meas <- tb_pred %>%
    mutate(
      y = y,
      pl = dt((y - f) / sqrt(q), df) / sqrt(q)
    ) %>%
    summarise(
      lpl = sum(log(pl), na.rm = TRUE),
      rmse = Metrics::rmse(y, f),
      mase = Metrics::mase(y, f),
      mape = Metrics::mape(y, f),
      mdae = Metrics::mdae(y, f),
      mae  = Metrics::mae(y, f)
    )

  measures <- c(rmse = meas$rmse,
                mase = meas$mase,
                mape = meas$mape,
                mdae = meas$mdae,
                mae = meas$mae)

  res <- (y - tb_pred$f) / sd(y - tb_pred$f)
  res < zoo::zoo(res, index = zoo::index(y_zoo))

  out <- list(
    call = match.call(),
    model = model,
    lpl = meas$lpl,
    residuals = res,
    measures = measures,
    var_law = var_law,
    pow = pow,
    delta_phi = delta_phi,
    prior = prior,
    state = tb_state,
    pred = tb_pred,
    final_state = final_state,
    y = y_zoo
  )

  class(out) <- "filterDist"
  out
}

#' @export

print.filterDist <- function(x, ...) {

  cat(paste0("Call:\n", deparse(x$call), "\n\n"))

  cat("Filter distribution", "\n\n")

  aux <- x$model$dim_p
  if (length(aux) > 1) cat(paste0(aux[1], " " ,names(aux)[1], " with ", aux[2], " components of ", names(aux)[2]), "\n \n")
  else cat(paste0(aux[1], names(aux)[1]))

  cat(paste0("Variance law: ", x$var_law), "\n")

  cat(paste0("Variance discount factor: ", x$delta_phi), "\n\n")

  cat(paste0("Predictive log-likelihod: ", round(x$lpl, 4)), "\n \n")

  cat("Measures of one-step ahead predictions errors:\n")
  print(round(x$measures, 4))
  cat("\n\n")

  qts <- quantile(x$residuals, probs = c(0.25, 0.5, 0.75))
  names(qts) <- c("1Q", "Median", "3Q")
  sts <- round(c(Min = min(x$residuals), qts, Max = max(x$residuals)), 4)
  cat("Residuals:\n")
  print(qts)

  invisible()
}
