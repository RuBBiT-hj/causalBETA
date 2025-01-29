#' Plot the Posterior Predictive Survival Curve for a Single individual
#' 
#' Plot the survival curve based on the predictions from the posterior draws for the specified individual
#' @param x an object of the class `bayeshaz` created by the `bayeshaz()` function
#' @param data the data frame for an individual used for prediction, with the variables specified in the model
#' @param n the number of prediction for each posterior draw; the default is 1000
#' @param col the color parameter for the baseline hazard points, default is `black`
#' @param col_CI the color parameter for the confidence intervals of the baseline hazard, default is semitransparent grey
#' @param xlim Limits for x-axis
#' @param ylim Limits for y-axis
#' @param type Plot type
#' @param pch Plot character
#' @param cex Character expansion factor
#' @param lwd Line width
#' @param ... other graphical parameters for the plot function. Default ones will be used if not provided.
#' 
#' @details
#' This function uses `predict.bayeshaz` to obtain all predictions for an individual,
#' and it generates a matrix contains the posterior predictive probability for each time interval.
#' 
#' It also makes a plot of survival curve, and the user can modify by providing plot parameters.
#' The plot has the the mean survival probability and the 95% CI from bootstrapping.
#' If the users want to make their own new plot, they can utilizes the returned matrix.
#' 
#' @returns
#' This function returns a matrix contains the posterior predictive probability for each time interval.
#' Each row represents a posterior draw, and each column corresponds to a time interval.
#' 
#' @examples
#' \dontrun{
#' # example demo
#' # after getting the posterior draw object
#' pred_survival <- plot_survival(post_draws_ar1_adj, df_veteran[1, ])
#' }
## usethis namespace: start
#' @import survival
#' @import graphics
#' @import grDevices
#' @importFrom mets rpch
## usethis namespace: end
#' @export


plot_survival = function(x, data, n=1000,
                         col = "black", col_CI = rgb(0.5, 0.5, 0.5, 0.5),
                         type = 'o', pch = 20, cex = 0.5, lwd = 1.5,
                         xlim = NULL, ylim = c(0,1),
                         ...){
  
  # extract
  beta_draws = x$beta_draws
  haz_draws = x$haz_draws
  partition = x$partition
  
  # call predict.bayeshaz as a helper function
  # the number of the list - the number of posterior draws
  # the length of each list - the number of predictions (n)
  all_surv_time = predict.bayeshaz(x, data, n, func = list)
  all_surv_time = all_surv_time[[1]] # as only a single individual
  
  # let t be the middle points of each partition (the number is the same as the intervals)
  t <- x$midpoint
  
  surv_prob <- matrix(nrow = length(all_surv_time),
                      ncol = length(t))
  
  # For each t, we calculate the proportion
  for (i in 1:length(t)){
    surv_prob[ ,i] <- sapply(all_surv_time, function(surv_time) mean(surv_time > t[i]))
  }
  
  # Upper and lower quantiles
  lwr <- apply(surv_prob, 2, quantile, probs=.025)
  upr <- apply(surv_prob, 2, quantile, probs=.975)
  
  # graphic parameters checking for necessary ones
  if (is.null(xlim)) xlim = c(0, max(partition))
  if (is.null(ylim)) ylim = c(0, 1) # should be between 0 and 1
  
  plot(c(0, t), c(1, colMeans(surv_prob)), pch=pch, ylim=ylim,
       type = type, cex = cex, ...)
  segments(x0 = t, y0 = lwr,
           x1 = t, y1 = upr,
           lwd = lwd, col = col_CI)
  
  return(surv_prob)
}
