% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast.R
\name{forecast}
\alias{forecast}
\title{Forecast function for Dynamic Models}
\usage{
forecast(object, horizon, distr = c("state", "predictive", "both"))
}
\arguments{
\item{object}{\code{dm} class object from \code{fit.dm} function.}

\item{horizon}{horizon of the forecast}

\item{distr}{forecast the predictive, state or both distributions?}
}
\description{
Forecast univariate Dynamic Models.
}
\references{
West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
}
\author{
André F. B. Menezes \email{andrefelipemaringa@gmail.com}
}
