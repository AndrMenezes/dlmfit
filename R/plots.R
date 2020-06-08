# fit <- model_fit

# Plot of the one-step ahead predictive distribution ----------------------

plot.dm.predictive <- function(fit, interval = TRUE,
                               interval.type = "exact", interval.level = 0.05) {

  ## Number of parameters
  p <- length(fit$final_state$mt)

  df <- bind_cols(
    y = as.numeric(fit$y),
    date = index(fit$y),
    fit$pred
  ) %>%
    slice(-(1:(p + 1)))

  out <- ggplot(df, aes(x = date, y = y)) +
    geom_line() +
    geom_line(aes(y = f), col = "red") +
    theme_cowplot() +
    background_grid()

  if (interval) {
    df <- add_pi.dm(df, type = interval.type, level = interval.level)
    out <- out +
      geom_ribbon(data = df, aes(ymin = lw, ymax = up), alpha = 0.6, fill = "grey69")
  }

  out
}

# Plot of the parameter state distribution --------------------------------

plot.dm.state <- function(fit, which = "all",
                          interval.type = "exact", interval.level = 0.05) {

  p <- length(fit$final_state$mt)
  aux  <- tibble(t = fit$pred$t, date = index(fit$y))
  df <- fit$state %>%
    left_join(aux, by = "t") %>%
    filter(t > p + 1) %>%
    add_pi.dm(type = interval.type, level = interval.level)

  if(any(which == "all")) {
    out <- ggplot(df, aes(x = date, y = a)) +
      facet_wrap(~parameter, scales = "free_y") +
      geom_line() +
      geom_ribbon(aes(ymin = lw, ymax = up), alpha = 0.6, fill = "grey69") +
      theme_cowplot() +
      background_grid()
  } else {
    out <- df %>%
      filter(parameter %in% which) %>%
      ggplot(aes(x = date, y = a)) +
      facet_wrap(~parameter, scales = "free_y") +
      geom_line() +
      geom_ribbon(aes(ymin = lw, ymax = up), alpha = 0.6, fill = "grey69") +
      theme_cowplot() +
      background_grid()
  }
  out
}

# Plot of the forecast distribution ---------------------------------------

plot.dm.forecast <- function(fit, horizon, distr = "predictive", which.state = "all",
                             interval.type = "exact", interval.level = 0.05) {
  ## number of parameters
  p <- length(fit$final_state$mt)

  if (distr == "predictive") {
    df <- bind_cols(
      y = as.numeric(fit$y),
      date = index(fit$y)
    )
    fore <- fit %>%
      forecast.dm(horizon = horizon, which = "predictive") %>%
      add_pi.dm(type = interval.type, level = interval.level) %>%
      rename(date = t)

    out <- ggplot() +
      geom_line(data = df, aes(x = date, y = y)) +
      geom_line(data = fore, aes(x = date, y = a), col = "#0000AA") +
      geom_ribbon(data = fore, aes(x = date, ymin = lw, ymax = up),
                  fill = "grey69", alpha = 0.6) +
      theme_cowplot() +
      background_grid()
  } else {
      fore <- fit %>%
        forecast.dm(horizon = horizon, which = "state") %>%
        add_pi.dm(type = interval.type, level = interval.level) %>%
        rename(date = t)
      if(which.state != "all") {
        fore %>%
          filter(parameter %in% which.state)
      }

      out <- ggplot(fore, aes(x = date, y = a)) +
        facet_wrap(~parameter, scales = "free_y") +
        geom_line(col = "#0000AA") +
        geom_ribbon(aes(ymin = lw, ymax = up), fill = "grey69", alpha = 0.6) +
        theme_cowplot() +
        background_grid()
    }
  out
}


