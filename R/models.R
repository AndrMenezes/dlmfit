poly.dm <- function(order, delta=0.98) {
  out <- dlm::dlmModPoly(order = order)[c("FF", "GG")]
  out$dim_p <- ncol(out$GG)
  out$D <- matrix(1/delta, nrow = out$dim_p, ncol = out$dim_p)
  out$mod <- "polynomial"
  class(out) <- "dm"
  out
}
seas.dm <- function(frequency, delta=0.96) {
  out <- dlm::dlmModSeas(frequency = frequency)[c("FF", "GG")]
  out$dim_p <- ncol(out$GG)
  out$D <- matrix(1/delta, nrow = out$dim_p, ncol = out$dim_p)
  out$mod <- "seasonal_free"
  class(out) <- "dm"
  out
}
trig.dm <- function(s, delta=0.96) {
  out   <- dlm::dlmModTrig(s = s)[c("FF", "GG")]
  out$dim_p <- ncol(out$GG)
  out$D <- matrix(1/delta, nrow = out$dim_p, ncol = out$dim_p)
  out$mod <- "seasonal_fourier"
  class(out) <- "dm"
  out
}
superposition.dm <- function(..., lt) {
  if(missing(lt)) lt <- list(...)

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

  class(out) <- "dm"
  out
}
