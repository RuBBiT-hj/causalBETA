#' Plot Average Baseline Hazard Rates
#' 
#' Plot the average baseline hazard rate from all posterior draws on all intervals from the given `bayeshaz` object.
#' 
#' @param bayeshaz_object an object of the class `bayeshaz` created by the `bayeshaz()` function
#' @param col_hazard the color parameter for the baseline hazard points, default is `black`
#' @param col_CI the color parameter for the confidence intervals of the baseline hazard, default is semitransparent-grey
#' @param level_CI Credible interval level, for a specified value an equal-tailed, 
#' level_CI% credible interval will be plotted which has ((1-level_CI*100)/2)% 
#' posterior probability below and above the interval. 
#' E.g. level_CI=.95 (the default) plots a 95% credible interval.
#' @param ... other graphical parameters for the plot function. Default ones will be used if not provided.
#' 
#' @description
#' This function plots the baseline hazard rate at the midpoint for each interval,
#' and it also marks the number at risk at the bottom corresponding to the ticks of the time axis.
#' To enable the plot with two axes to show both time and the number at risk but also with a clean display,
#' we limit the freedom of changing any graphical parameters, i.e. some will take the default values overriding `NULL`.
#' If the users don't like the output format, they can extract the parameters from `bayeshaz` object directly
#' to generate plot(s).
#' 
#' 
#' @examples
#' # example demo
#' ## Continued from ?bayeshaz
#' set.seed(1)
#' post_draws_ind = bayeshaz(
#'   d = data, ## data set
#'   reg_formula = Surv(y, delta) ~ A,
#'   num_partitions = 100, 
#'   model = 'independent',
#'   sigma = 3,
#'   A = 'A',
#'   warmup = 1000,
#'   post_iter = 1000)
#' plot(post_draws_ind, ylim=c(0,.11),
#'   xlim=c(0, 900),
#'   type='p',
#'   main='Independent Prior Process',
#'   ylab = 'Baseline Hazard Rate', 
#'   xlab = 'Time (days)')
## usethis namespace: start
## usethis namespace: end
#' @export

plot.bayeshaz = function(bayeshaz_object, col_hazard = "black", 
                         col_CI = rgb(0.5, 0.5, 0.5, 0.5), level_CI=.95,
                         type = 's', pch = 20, xlim = NULL, ylim = NULL,
                         xlab = "Time", ylab = NULL, main = NULL, cex = 0.5,
                         lwd = 1.5,
                         ...){
 # extract the key components from the object
  hazard <- do.call(rbind, bayeshaz_object$haz_draws)
  
  partitions <- bayeshaz_object$partition
  xv <- bayeshaz_object$midpoint
  obs_time <- bayeshaz_object$data[, bayeshaz_object$time]
  
  # Generate the vector for time axis
  endpoint <- max(partitions)
  time_unit <- 10^(floor(log10(endpoint) - 1))
  time_coef <- round(endpoint / time_unit / 10)
  time_unit <- time_coef * time_unit
  time_axis <- seq(0, endpoint, by = time_unit)
  
  # Generate the vector for number of subjects axis
  n_obs_axis <- c()
  for (i in 1:length(time_axis)){
    n_obs_axis[i] <- sum(obs_time > time_axis[i])
  }
  
  ## posterior mean/ 95% interval
  bslhaz_mean = colMeans(hazard)
  
  tail_prob = (1-level_CI) / 2
  bslhaz_lwr = apply(hazard, 2, quantile, probs= tail_prob )
  bslhaz_upr = apply(hazard, 2, quantile, probs= 1 - tail_prob )
  
  # graphic parameters checking for necessary ones
  if (is.null(xlim)) xlim = c(0, endpoint)
  if (is.null(ylim)) ylim = c(0, 1.1*max(bslhaz_upr))
  if (is.null(ylab)) ylab = "Average Baseline Hazard Rates"
  if (is.null(main)) main = "Average Baseline Hazard Rates Over Time"
  
  # Plot
  par(mar = c(6.1, 4.1, 3.1, 2.1),
      mgp = c(3,1.5,0))
  
  plot( xv, bslhaz_mean, pch=pch, 
        xlim = xlim, ylim=ylim,
        type=type, cex = cex, col = col_hazard,
        xlab = "", xaxt = "n",
        main = main, ylab = ylab, ...) 
  segments(x0 = xv, y0 = bslhaz_lwr,
           x1 = xv, y1 = bslhaz_upr,
           lwd = lwd, col = col_CI)
  axis(1, time_axis, paste0(time_axis, "\n (n=", n_obs_axis, ")"),
       line = 1, tck = -0.03)
  title(xlab = xlab, line = 4)
  
  # reset
  par(mar = c(5.1, 4.1, 4.1, 2.1),
      mgp = c(3,1,0))
}