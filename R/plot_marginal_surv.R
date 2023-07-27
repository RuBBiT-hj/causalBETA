#' Marginal Survival Curves Predicted
#' 
#' Plot the marginal survival curves that are based on the results from ATE_estimation() function.
#' If both treatments are chosen, the difference between two curves would be the average treatment effect (ATE) estimated. It could be plotted with plot_ATE() function
#' @param ATE_estimated the list with the format: t, surv_prob_1, surv_prob_2, and ATE.
#' @param trt_name the name of the treatment variable
#' @param trt_values the value of the treatment variable, and the order should be the same as the ATE_estimation() function 
#' @param trt_option select treatment 1 or 2 or both. Default is c(1,2)
#' @examples
#' # example demo
## usethis namespace: start
## usethis namespace: end
#' @export

plot_marginal_surv = function(ATE_estimated, trt_name, trt_values, trt_option = c(1,2)){
  # read from the input
  t <- ATE_estimated[[1]]
  surv_1 <- ATE_estimated[[2]]
  surv_2 <- ATE_estimated[[3]]
  
  title_trt <- paste0("Marginal Survival Curves for ", trt_name)
  
  lwr_1 <- apply(surv_1, 2, quantile, probs=.025)
  upr_1 <- apply(surv_1, 2, quantile, probs=.975)
  
  lwr_2 <- apply(surv_2, 2, quantile, probs=.025)
  upr_2 <- apply(surv_2, 2, quantile, probs=.975)
  
  # different options from trt_option
  if ( !all(trt_option %in% c(1,2)) ) stop("trt_option given in not valid. 
                                        Valid values are 1, 2, or c(1,2).")
  if (all(c(1,2) %in% trt_option)){
    plot(t, colMeans(surv_1), pch=20, ylim=c(0, 1),
         type = "o", cex = 0.5,  col = rgb(0.9, 0.1, 0.1, 0.9),
         xlab = "", ylab = "")
    segments(x0 = t, y0 = lwr_1,
             x1 = t, y1 = upr_1,
             lwd = 1.5, col = rgb(0.9, 0.5, 0.5, 0.4))
    lines(t, colMeans(surv_2), pch=20,
          type = "o", cex = 0.5,  col = rgb(0.1, 0.1, 0.9, 0.9))
    segments(x0 = t, y0 = lwr_2,
             x1 = t, y1 = upr_2,
             lwd = 1.5, col = rgb(0.5, 0.5, 0.9, 0.4))
    title(xlab = "Time", ylab = "Predicted Marginal Survival Prob",
          main = title_trt)
    legend("topright", 
           legend = c(paste0(trt_values[1]), paste0(trt_values[2])),
           lty = c(1,1),
           col = c(rgb(0.9, 0.1, 0.1, 0.9),  rgb(0.1, 0.1, 0.9, 0.9)))
  }
  else if (trt_option == 1){
    plot(t, colMeans(surv_1), pch=20, ylim=c(0, 1),
         type =  "o", cex = 0.5,
         col = rgb(0.9, 0.1, 0.1, 0.9),
         xlab = "", ylab = "")
    segments(x0 = t, y0 = lwr_1,
             x1 = t, y1 = upr_1,
             lwd = 1.5, col = rgb(0.9, 0.5, 0.5, 0.4))
    title(xlab = "Time", ylab = "Predicted Marginal Survival Prob",
          main = paste0(title_trt, " ", trt_values[1]))
  } 
  else if (trt_option == 2){
    plot(t, colMeans(surv_2), pch=20, ylim=c(0, 1),
         type = "o", cex = 0.5,
         col = rgb(0.1, 0.1, 0.9, 0.9),
         xlab = "", ylab = "")
    segments(x0 = t, y0 = lwr_2,
             x1 = t, y1 = upr_2,
             lwd = 1.5, col = rgb(0.5, 0.5, 0.9, 0.4))
    title(xlab = "Time", ylab = "Predicted Marginal Survival Prob",
          main = paste0(title_trt, " ", trt_values[2]))
  }
  else {
    stop("use for debug") # need to delete in the future
  }
  
  # I need to modify the if-else statements of thinking about the length of option?

}
