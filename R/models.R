poly.dm <- function(order, delta=0.98) {
  out <- dlm::dlmModPoly(order = order)[c("FF", "GG")]
  out$delta <- delta
  out$mod <- "polynomial"
  class(out) <- "dm"
  out
}
seas.dm <- function(frequency, delta=0.96) {
  out <- dlm::dlmModSeas(frequency = frequency)[c("FF", "GG")]
  out$delta <- delta
  out$mod <- "seasonal_free"
  class(out) <- "dm"
  out
}
trig.dm <- function(s, delta=0.96) {
  out   <- dlm::dlmModTrig(s = s)[c("FF", "GG")]
  out$delta <- delta
  out$mod <- "seasonal_fourier"
  class(out) <- "dm"
  out
}
superposition.dm <- function(..., lt) {
  if(missing(lt)) lt <- list(...)

  FF  <- purrr::map(lt, "FF")
  GG  <- purrr::map(lt, "GG")
  mod <- purrr::map(lt, "mod")
  delta <- do.call("c", purrr::map(lt, "delta"))
  names(delta) <- mod
  mod <- paste(mod)

  out <- list(
    FF = FF,
    GG = GG,
    delta = delta,
    mod = mod,
    superposition = TRUE
  )

  class(out) <- "dm"
  out
}
