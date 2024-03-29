% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.bayeshaz.R
\name{plot.bayeshaz}
\alias{plot.bayeshaz}
\title{Plot Average Baseline Hazard Rates}
\usage{
\method{plot}{bayeshaz}(
  bayeshaz_object,
  col_hazard = "black",
  col_CI = rgb(0.5, 0.5, 0.5, 0.5),
  level_CI = 0.95,
  type = "s",
  pch = 20,
  xlim = NULL,
  ylim = NULL,
  xlab = "Time",
  ylab = NULL,
  main = NULL,
  cex = 0.5,
  lwd = 1.5,
  ...
)
}
\arguments{
\item{bayeshaz_object}{an object of the class \code{bayeshaz} created by the \code{bayeshaz()} function}

\item{col_hazard}{the color parameter for the baseline hazard points, default is \code{black}}

\item{col_CI}{the color parameter for the confidence intervals of the baseline hazard, default is semitransparent-grey}

\item{level_CI}{Credible interval level, for a specified value an equal-tailed,
level_CI\% credible interval will be plotted which has ((1-level_CI*100)/2)\%
posterior probability below and above the interval.
E.g. level_CI=.95 (the default) plots a 95\% credible interval.}

\item{...}{other graphical parameters for the plot function. Default ones will be used if not provided.}
}
\description{
This function plots the baseline hazard rate at the midpoint for each interval,
and it also marks the number at risk at the bottom corresponding to the ticks of the time axis.
To enable the plot with two axes to show both time and the number at risk but also with a clean display,
we limit the freedom of changing any graphical parameters, i.e. some will take the default values overriding \code{NULL}.
If the users don't like the output format, they can extract the parameters from \code{bayeshaz} object directly
to generate plot(s).
}
\details{
Plot the average baseline hazard rate from all posterior draws on all intervals from the given \code{bayeshaz} object.
}
\examples{
# example demo
## Continued from ?bayeshaz
set.seed(1)
post_draws_ind = bayeshaz(
  d = data, ## data set
  reg_formula = Surv(y, delta) ~ A,
  num_partitions = 100, 
  model = 'independent',
  sigma = 3,
  A = 'A',
  warmup = 1000,
  post_iter = 1000)
plot(post_draws_ind, ylim=c(0,.11),
  xlim=c(0, 900),
  type='p',
  main='Independent Prior Process',
  ylab = 'Baseline Hazard Rate', 
  xlab = 'Time (days)')
}
