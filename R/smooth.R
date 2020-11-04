#' @title Smooth distribution
#'
#' @description Get the smooth distribution for univariate Dynamic Models.
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
#' @rdname kalmanDist
#' @export

smoothDist <- function(model, y, prior=NULL, var_law="identity", pow=1, delta_phi=1) {

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

  ## save results
  l_m <- list()
  l_a <- list()
  l_C <- list()
  l_R <- list()
  st <- c()

  # Get the filter distribution
  for ( t in 1:nt ) {
    ## prior distribution in t-1
    P <- tcrossprod(GG %*% C0, GG)
    R <- D * P
    a <- GG %*% m0

    l_a[[t]] <- a
    l_R[[t]] <- R

    out <- update_moments(y[t], t,
                          FF, GG,
                          a, R,
                          n0, d0,
                          alpha=delta_phi, law=var_law)
    m0 <- out[["m"]]
    C0 <- out[["C"]]
    d0 <- out[["d"]]
    n0 <- out[["n"]]
    l_m[[t]] <- m0
    l_C[[t]] <- C0
    st[t] <- d0 / n0
  }

  t <- nt
  ## a_t(0)
  at_k <- l_m[[t]]
  ## R_t(0)
  Rt_k <- l_C[[t]]

  tb_state <- tibble()
  tb_mean  <- tibble()
  # Get the smooth distribution
  for (k in 1:(t-1)) {
    ## R_{t-k}
    Rt <- D * tcrossprod(GG %*% l_C[[t-k]], GG)
    ## B_{t-k}
    Bt_k <- tcrossprod(l_C[[t-k]],  GG) %*% solve(Rt)
    ## a_t(-k) e R_t(-k)        a_t(-1) - a_{t}
    at_k <- l_m[[t-k]] + Bt_k %*% (at_k - l_a[[t-k+1]])
    Rt_k <- l_C[[t-k]] + tcrossprod(Bt_k %*% (Rt_k - Rt), Bt_k)
    ## f_t(-k) e q_t(-k)
    ft_k <- crossprod(FF, at_k)
    qt_k <- crossprod(FF, Rt_k %*% FF)

    tb_state <- bind_rows(
      tb_state,
      data.frame(
        t = t-k,
        at_k = c(at_k),
        Rt_k = diag(Rt_k),
        df =  nt,
        parameter = paste0("theta_",1:p)
      )
    )
    tb_mean <- bind_rows(
      tb_mean,
      data.frame(
        t = t-k,
        ft_k = c(ft_k),
        qt_k = c(qt_k),
        df =  nt
      )
    )
  }
  tb_state <- dplyr::arrange(tb_state, t)
  tb_mean <- dplyr::arrange(tb_mean, t)
  out <- list(
    call = match.call(),
    y = y_zoo,
    model = model,
    var_law = var_law,
    pow = pow,
    delta_phi = delta_phi,
    prior = prior,
    state = tb_state,
    mean = tb_mean
  )
  class(out) <- "smoothDist"
  out
}

#' @export

print.smoothDist <- function(x, ...) {

  cat(paste0("Call:\n", deparse(x$call), "\n\n"))

  cat("Smooth distribution", "\n\n")
  aux <- x$model$dim_p
  if (length(aux) > 1) cat(paste0(aux[1], " " ,names(aux)[1], " with ", aux[2], " components of ", names(aux)[2]), "\n \n")
  else cat(paste0(aux[1], names(aux)[1]))

  cat(paste0("Variance law: ", x$var_law), "\n")

  cat(paste0("Variance discount factor: ", x$delta_phi), "\n\n")

  invisible()
}


