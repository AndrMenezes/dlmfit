fit.dm <- function(model, y, prior=NULL, var_law="identity", pow=1, delta_phi=1) {

  ## Define some objects
  if (is.null(model$superposition)) {
    FF <- t(model$FF)
    GG <- model$GG
    dim_p <- nrow(GG)
    D <- matrix(1/model$delta,nrow=dim_p, ncol=dim_p)
  } else {
    FF <- t(do.call("cbind", model$FF))
    GG <- do.call("bdiag", model$GG)
    dim_p <- sapply(model$FF, ncol)
    delta <- model$delta
    aux <- sapply(1:length(delta), function(i) {
      matrix(1/delta[i], nrow=dim_p[i],ncol=dim_p[i])
      #diag(1/delta[i], nrow = dim_p[i], ncol=dim_p[i])
    })
    D <- do.call("bdiag", aux)
  }


  ## Define prior distribution
  p  <- nrow(FF)
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
                          n0, d0,
                          alpha=delta_phi, law=var_law)
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

  lpl <- tb_pred %>%
    mutate(
      y = y,
      pl = dt((y - f) / sqrt(q), df) / sqrt(q)
    ) %>%
    summarise(sum(log(pl), na.rm = TRUE)) %>%
    pull()


  out <- list(
    y = y_zoo,
    model = model,
    lpl = lpl,
    var_law = var_law,
    pow = pow,
    delta_phi = delta_phi,
    prior = prior,
    state = tb_state,
    pred = tb_pred,
    final_state = final_state
  )
  class(out) <- "dm"
  out
}







# smooth.dm <- function(formula, y, prior) {
#
# }
