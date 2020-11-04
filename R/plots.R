#' @name plots
#' @aliases PlotPredictive PlotState PlotForecast PlotSmooth
#'
#' @title Plottting functions
#'
#' @description Functions to plot the results of a model fit or smooth.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object \code{dm} class object from \code{fit} or \code{smooth} functions.
#' @param horizon horizon of forecast
#' @param which.distr which distribution to plot
#' @param which.state which state to plot
#' @param interval logical indicating if the interval will be plot
#' @param interval.type which type of interval to use
#' @param interval.level which significance level to use
#' @param geom_data if the observed data will be represent by point or lines
#' @param x não sei
#' @param ... currently not used.
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom utils tail

plot.dm <- function(x, ...) {

  mod <- x$model$mod
  p   <- x$model$dim_p
  if (mod == "polynomial") {
    # Quais plots fazer para modelos só polinomiais?
    p1 <- PlotPredictive(x)
    p2 <- PlotState(x, which.state = "theta_1")
    #p3 <- PlotResiduals(x)
    if (p == 1) {


    }

  }
  else if (mod %in% c("seasonal_free", "seasonal_fourier")) {

  }
  else if (grepl("+", mod, fixed = T) ) {
    if (any(mod == "regression")) {

    } else {



    }

  }

}


# Plot of the one-step ahead predictive distribution ----------------------
#' @export
#' @rdname plots
PlotPredictive <- function(object, interval = TRUE,
                               geom_data = "line",
                               interval.type = "exact",
                               interval.level = 0.05) {

  ## Number of parameters
  p <- length(object$final_state$mt)
  mark_t <- min(index(object$y)) + diff(tail(index(object$y), 2)) * p

  df <- dplyr::bind_cols(
    y = as.numeric(object$y),
    date = index(object$y),
    object$pred
  ) %>%
    dplyr::mutate(
      f = ifelse(t <= p +1, NA_real_, f)
    )
  if(geom_data == "line") {
    out <- ggplot(df, aes(x = date, y = y)) +
      geom_line()
  } else {
    out <- ggplot(df, aes(x = date, y = y)) +
      geom_point()
  }

  out <- out +
    geom_line(aes(y = f), col = "#0000AA", size = 1.0) +
    geom_vline(xintercept = mark_t, linetype = "dashed") +
    theme_cowplot() +
    background_grid()

  if (interval) {
    df <- add_pi(df, type = interval.type, level = interval.level)
    out <- out +
      geom_ribbon(data = df, aes(ymin = lw, ymax = up), alpha = 0.6, fill = "grey69")
  }

  out
}

# Plot of the parameter state distribution --------------------------------
#' @export
#' @rdname plots
PlotState <- function(object, which.state = "all", interval.type = "exact", interval.level = 0.05) {

  p <- length(object$final_state$mt)
  aux <- tibble(t = object$pred$t, date = index(object$y))
  df <- object$state %>%
    dplyr::left_join(aux, by = "t") %>%
    dplyr::filter(t > p + 1) %>%
    add_pi(type = interval.type, level = interval.level)

  if (any(which.state == "all")) {
    out <- ggplot(df, aes(x = date, y = a)) +
      facet_wrap(~parameter, scales = "free_y") +
      geom_line() +
      geom_ribbon(aes(ymin = lw, ymax = up), alpha = 0.6, fill = "grey69") +
      theme_cowplot() +
      background_grid()
  } else {
    out <- df %>%
      dplyr::filter(parameter %in% which.state) %>%
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
#' @export
#' @rdname plots
PlotForecast <- function(object, horizon, which.distr = "predictive", which.state = "all",
                         interval.type = "exact", interval.level = 0.05) {
  ## number of parameters
  p <- length(object$final_state$mt)

  if (which.distr == "predictive") {
    df <- dplyr::bind_cols(
      y = as.numeric(object$y),
      date = index(object$y)
    )
    fore <- object %>%
      forecast(horizon = horizon, distr = "predictive") %>%
      add_pi(type = interval.type, level = interval.level) %>%
      dplyr::rename(date = t)
    mark_t <- min(fore$date)

    out <- ggplot() +
      geom_line(data=df, aes(x = date, y = y)) +
      geom_line(data = fore, aes(x = date, y = f), col = "#0000AA", size = 1) +
      geom_ribbon(data = fore, aes(x = date, ymin = lw, ymax = up),
                  fill = "grey69", alpha = 0.6) +
      geom_vline(xintercept = mark_t, linetype = "dashed") +
      theme_cowplot() +
      background_grid()
  } else {
      fore <- object %>%
        forecast(horizon = horizon, distr = "state") %>%
        add_pi(type = interval.type, level = interval.level) %>%
        dplyr::rename(date = t)
      if (any(which.state != "all")) {
        fore <- fore %>%
          dplyr::filter(parameter %in% which.state)
      }

      out <- ggplot(fore, aes(x = date, y = a)) +
        facet_wrap(~parameter, scales = "free_y") +
        geom_line() +
        geom_ribbon(aes(ymin = lw, ymax = up), fill = "grey69", alpha = 0.6) +
        theme_cowplot() +
        background_grid()
    }
  out
}


# Plot of the standardized residuals --------------------------------------
#' @export
#' @rdname plots
PlotResiduals <- function(object) {



}


# Plot of the smooth distribution -----------------------------------------
#' @export
#' @rdname plots
PlotSmooth <- function(object, which.distr = "mean", which.state = "all",
                       interval.type = "exact", interval.level = 0.05) {

  p <- length(object$prior$m0)

  if (which.distr == "mean") {
    df <- dplyr::bind_cols(
      y = as.numeric(object$y),
      date = index(object$y)
    ) %>%
      dplyr::slice(-n()) %>%
      dplyr::bind_cols(
        object$mean %>%
          dplyr::rename(f = ft_k, q = qt_k) %>%
          add_pi(type = interval.type, level = interval.level)
      )
    out <- ggplot(data=df, aes(x = date, y = y)) +
      geom_line() +
      geom_line(aes(y = f), col = "#0000AA", size = 1) +
      geom_ribbon(aes(ymin = lw, ymax = up), fill = "grey69", alpha = 0.6) +
      theme_cowplot() +
      background_grid()
  } else {
    n <- length(object$y)
    aux  <- dplyr::tibble(t = 1:(n-1), date = index(object$y)[-n])
    df <- object$state %>%
      dplyr::arrange(t) %>%
      dplyr::left_join(aux, by = "t") %>%
      dplyr::rename(a = at_k, dR = Rt_k) %>%
      add_pi(type = interval.type, level = interval.level)

    if (any(which.state != "all")) {
      df <- df %>%
        dplyr::filter(parameter %in% which.state)
    }

    out <- ggplot(df, aes(x = date, y = a)) +
      facet_wrap(~parameter, scales = "free_y") +
      geom_line() +
      geom_ribbon(aes(ymin = lw, ymax = up), fill = "grey69", alpha = 0.6) +
      theme_cowplot() +
      background_grid()
  }
  out
}

