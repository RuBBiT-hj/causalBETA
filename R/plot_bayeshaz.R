#' Plot Average Baseline Hazard Rates
#' 
#' Plot the average baseline hazard rate from all posterior draws on all intervals generated from the bayeshaz() function
#' @param hazard a vector containing baseline hazard rates from all posterior draws for each interval
#' @param partitions a vector of time values for partitions
#' @param xv a vector of middle points for intervals. It is an returned output from bayeshaz()
#' @param obs_time a vector of the observed time of all subjects before events or censoring

#' @examples
#' # example demo
## usethis namespace: start
## usethis namespace: end
#' @export

plot_bayeshaz = function(hazard, partitions, xv, obs_time){
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
  bslhaz_lwr = apply(hazard, 2, quantile, probs=.025)
  bslhaz_upr = apply(hazard, 2, quantile, probs=.975)
  
  # Plot
  par(mar = c(6.1, 4.1, 3.1, 2.1),
      mgp = c(3,1.5,0))
  plot( xv, bslhaz_mean, pch=20, ylim=c(0, 1.1*max(bslhaz_upr)),
        type = "o", cex = 0.5,
        xlab = "", xaxt = "n",
        main = "Average Baseline Hazard Rates Over Time",
        ylab = "Average Baseline Hazard Rates") 
  segments(x0 = xv, y0 = bslhaz_lwr,
           x1 = xv, y1 = bslhaz_upr,
           lwd = 1.5, col = rgb(0.5, 0.5, 0.5, 0.5))
  axis(1, time_axis, paste0(time_axis, "\n (n=", n_obs_axis, ")"),
       line = 1, tck = -0.03)
  title(xlab = "Time", line = 4)
  
  # reset
  par(mar = c(5.1, 4.1, 4.1, 2.1),
      mgp = c(3,1,0))
}

