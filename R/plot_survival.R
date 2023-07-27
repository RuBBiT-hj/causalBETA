#' Plot Survival Curve Predicted
#' 
#' Plot the survival curve based on the predictions from the posterior draw
#' @param x the individual used for prediction, with the variables specified in the model
#' @param beta_draws the beta coefficients from the posterior draws, with dimension number of posterior draws * the number of variables 
#' @param haz_draws the baseline hazard rates from the posterior draws, with dimension number of posterior draws * the number of intervals
#' @param partition the time values for each partition
#' @param n the number of prediction for each posterior draw; the default is 1000
#' @examples
#' # example demo
## usethis namespace: start
#' @import survival
#' @importFrom mets rpch
## usethis namespace: end
#' @export


plot_survival = function(x, beta_draws, haz_draws, partition, n=1000){
  
  # call predict.haz as a helper function
  # the number of the list - the number of posterior draws
  # the length of each list - the number of predictions (n)
  all_surv_time = predict.haz(x, beta_draws, haz_draws, partition, n, func = list)
  
  # let t be the middle points of each partition (the number is the same as the intervals)
  t <- diff(partition)/2 + partition[-length(partition)]
  
  surv_prob <- matrix(nrow = length(all_surv_time),
                      ncol = length(t))
  
  # For each t, we calculate the proportion
  for (i in 1:length(t)){
    surv_prob[ ,i] <- sapply(all_surv_time, function(x) mean(x > t[i]))
  }
  
  # Upper and lower quantiles
  lwr <- apply(surv_prob, 2, quantile, probs=.025)
  upr <- apply(surv_prob, 2, quantile, probs=.975)
  
  plot(c(0, t), c(1, colMeans(surv_prob)), pch=20, ylim=c(0, 1),
       type = "o", cex = 0.5,
       xlab = "Time", ylab = "Predicted Survival Probability")
  segments(x0 = t, y0 = lwr,
           x1 = t, y1 = upr,
           lwd = 1.5, col = rgb(0.5, 0.5, 0.5, 0.5))
  
  return(surv_prob)
}
