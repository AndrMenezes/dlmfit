add_pi.dm <- function(tb, type = "exact", level = 0.05) {

  if (!any(names(tb) == "parameter")) {
    tb <- tb %>%
      rename(a = f, dR = q)
  }

  if (type == "exact") {
    tb <- tb %>%
      mutate(
        lw = a + qt(level / 2, df = df) * sqrt(dR),
        up = a + qt(1-level / 2, df = df) * sqrt(dR)
      )
  }
  else if (type == "aprox") {
    tb <- tb %>%
      mutate(
        lw = a + qnorm(level / 2) * sqrt(dR),
        up = a + qnorm(1-level / 2) * sqrt(dR)
      )
  }

  if (!any(names(tb) == "parameter")) {
    tb <- tb %>%
      rename(f = a, q = dR)
  }

  tb
}
