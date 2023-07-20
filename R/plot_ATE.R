#' Plot Average Treatment Effect Estimated
#' 
#' Plot the average treatment effect (ATE) estimated from `ATE_estimation` function 
#' @param ATE_estimated the list with the format: t, time points; ATE, ATE estimated from many posterior draws

#' @examples
#' # example demo
## usethis namespace: start
## usethis namespace: end
#' @export

plot_ATE = function(ATE_estimated){
  t <- ATE_estimated[[1]]
  ATE <- ATE_estimated[[2]]
  # Upper and lower quantiles
  lwr <- apply(ATE, 2, quantile, probs=.025)
  upr <- apply(ATE, 2, quantile, probs=.975)
  
  plot(t, colMeans(ATE), pch=20, ylim=c(-1, 1),
       type = "o", cex = 0.5,
       xlab = "Time", ylab = "Predicted ATE")
  segments(x0 = t, y0 = lwr,
           x1 = t, y1 = upr,
           lwd = 1.5, col = rgb(0.5, 0.5, 0.5, 0.5))
}
