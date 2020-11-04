#' @name models
#' @aliases poly.dm seas.dm trig.dm superposition.dm
#'
#' @title Dynamic models
#'
#' @description Available Dynamic Models.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @param delta discount factor.
#' @param order order of the polynomial model.
#' @param frequency number of seasons to create seasonal factors representation.
#' @param s the period to create Fourier representation.
#' @param ... pass the models separate by commma to make the superposition.
#' @param lt list of models to make the superposition.
#'
#' @importFrom dlm dlmModPoly dlmModSeas dlmModTrig bdiag

#' @rdname models
#' @export
StatePoly <- function(order, delta=0.98) {
  out <- dlm::dlmModPoly(order = order)[c("FF", "GG")]
  out$dim_p <- ncol(out$GG)
  out$D <- matrix(1/delta, nrow = out$dim_p, ncol = out$dim_p)
  out$mod <- "polynomial"
  class(out) <- "dm"
  out
}
#' @rdname models
#' @export
StateSeas <- function(frequency, delta=0.96) {
  out <- dlm::dlmModSeas(frequency = frequency)[c("FF", "GG")]
  out$dim_p <- ncol(out$GG)
  out$D <- matrix(1/delta, nrow = out$dim_p, ncol = out$dim_p)
  out$mod <- "seasonal_free"
  out
}
#' @rdname models
#' @export
StateTrig <- function(s, delta=0.96) {
  out   <- dlm::dlmModTrig(s = s)[c("FF", "GG")]
  out$dim_p <- ncol(out$GG)
  out$D <- matrix(1/delta, nrow = out$dim_p, ncol = out$dim_p)
  out$mod <- "seasonal_fourier"
  out
}
#' @rdname models
#' @export
superposition <- function(..., lt) {
  if (missing(lt)) lt <- list(...)

  FF  <- do.call("cbind", purrr::map(lt, "FF"))
  GG  <- do.call("bdiag", purrr::map(lt, "GG"))
  D   <- do.call("bdiag", purrr::map(lt, "D"))
  mod <- paste(purrr::map(lt, "mod"))
  dim_p <- do.call("c", purrr::map(lt, "dim_p"))
  names(dim_p) <- mod

  out <- list(
    FF = FF,
    GG = GG,
    D = D,
    dim_p = dim_p,
    mod = paste(mod, collapse = " + "),
    superposition = TRUE
  )

  out
}
