#' @import ggplot2 magrittr dplyr
#'
#' @importFrom tidyr separate pivot_wider
#'
#' @name plotdlm
#'
#' @title Plots for fitted DLM from \code{dlmfit} function
#'
#' @description The function provides several plots of posteriori and predictive distributions.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @param mod the output from \code{dlmfit} function, a \code{tsibble} object.
#'
#' @param which a character indicating the type of plot. The options are
#' (i) \code{level} for the posterior mean level with confidence limits and the observed data,
#' (ii) \code{pred} for the mean of predictive distribution with confidence limits and the observed data,
#' (iii) \code{theta} for the mean posterior paramters with confidence limits and
#' (vi) \code{all} for all plots.
#'
#'
#' @return The plot(s) according to \code{which} option.
#'
#' @seealso \code{\link{dlmfit}}
#'
#' @examples
#'
#' y <- AirPassengers
#' fit_poly = dlmfit(y = y, order = 2, delta = 0.95)
#' plotdlm(fit_poly, which="level")
#' plotdlm(fit_poly, which="pred")
#' plotdlm(fit_poly, which="thetas")
#'
#' fit_seas <- dlmfit(y = y, order = 2, freq = 3, seasonality = "fourier", delta = c(0.95, 0.95))
#' plotdlm(fit_seas, which="level")
#' plotdlm(fit_seas, which="pred")
#' plotdlm(fit_seas, which="thetas")
#'
#' @rdname plotdlm
#' @export

plotdlm <- function(mod, which = c("level", "thetas", "pred", "all"))
{

  ## Nível da série (mu_t)
  p_mu <- mod %>%
    # dplyr::filter(dC < 1e5) %>%
    ggplot(aes(x = index, y = value)) +
    geom_line() +
    geom_line(aes(y = theta_1), col = "blue") +
    geom_line(aes(y = Up_1), col = "red", linetype = "dashed") +
    geom_line(aes(y = Lw_1), col = "red", linetype = "dashed") +
    theme_bw()

  ## Outros parâmetros
  p_thetas <- mod %>%
    dplyr::select(-c(f, Q, Lw, Up)) %>%
    tidyr::gather(parms, m, -c(index, value)) %>%
    tidyr::separate(parms, c("est", "parms")) %>%
    dplyr::mutate(parms = paste0("theta_", parms)) %>%
    tidyr::pivot_wider(names_from = est, values_from = m) %>%
    # dplyr::filter(dC < 1e4) %>%
    ggplot(aes(x = index, y = theta)) +
    facet_wrap(.~ parms, scales = "free") +
    geom_line() +
    geom_line(aes(y = Up), col = "red", linetype = "dashed") +
    geom_line(aes(y = Lw), col = "red", linetype = "dashed") +
    theme_bw() +
    labs(y = "Posterior mean")

  p_pred <- mod %>%
    # dplyr::mutate(f = c(f[-1], NA_real_),
    #               Lw = c(Lw[-1], NA_real_),
    #               Up = c(Up[-1], NA_real_)) %>%
    # dplyr::filter(Q < 1e4) %>%
    ggplot(aes(x = index, y = value)) +
    geom_line() +
    geom_line(aes(y = f), col = "blue") +
    geom_line(aes(y = Up), col = "red", linetype = "dashed") +
    geom_line(aes(y = Lw), col = "red", linetype = "dashed") +
    theme_bw()

  out <- switch(which,
               "all" = {list(p_mu = p_mu, p_thetas = p_thetas, p_pred = p_pred)},
               "level" = {p_mu},
               "thetas"= {p_thetas},
               "pred" = {p_pred})
  out
}
