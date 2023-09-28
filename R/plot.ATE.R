#' Plot Average Treatment Effect or Marginal Survival Probability Estimated for an `ATE` Object
#' 
#' This function displays the information stored in the `ATE` object: it can plot the marginal survival probability of
#' a treatment option or overlay both, and it can plot the ATE estimated corresponding to the reference value,
#' depending on the plotting choices user provided. It has some default plotting parameters if users don't specify otherwise.
#' 
#' @param ATE_object ATE_object an object of the class `ATE` created by the `ATE_estimation()` function
#' @param mode Mode of the plot, the default is `"ATE"` to plot the ATE in the `ATE` object. The other choices are
#' `0`, `1`, and `c(0,1)`, which plot the control, treatment, and both respectively.  
#' @param col_ATE,col_0,col_1 the colors for the ATE, control and treatment points and lines. The default are
#' black, red, and blue respectively
#' @param col_CI_ATE,col_CI_0,col_CI_1 the colors for the confidence interval of ATE, control and treatment.
#' The default are semitransparent grey, red, and blue
#' @param ... other other graphical parameters for the plot function. Default ones will be used if not provided.
#' @details
#' When the mode is `ATE`, this function plots the ATE and its 95% CI, and the default color is black and grey.
#' The default range of y-axis is from -1 to 1. The default title is `ATE Estimated Over Time`.
#' 
#' When the mode is `0` or `1`, this function plots the marginal survival probability of the control or treatment solely with
#' the 95% CI. The default colors are red and blue, and the default title is `Marginal Survival Curve for ...`.
#' 
#' When the mode is `c(0,1)` (or `c(1,0)`), this function plots the control and treatment together with a default legend.
#' The default title is `Marginal Survival Curves for Both Treatments`.
#' 
#' @examples
#' # example demo
## usethis namespace: start
## usethis namespace: end
#' @export

plot.ATE = function(ATE_object, mode = "ATE",
                    col_ATE = "black", col_0 = rgb(0.9, 0.1, 0.1, 0.9), col_1 = rgb(0.1, 0.1, 0.9, 0.9),
                    col_CI_ATE = rgb(0.5, 0.5, 0.5, 0.5), col_CI_0 = rgb(0.9, 0.5, 0.5, 0.4),
                    col_CI_1 = rgb(0.5, 0.5, 0.9, 0.4),
                    xlim = NULL, ylim = NULL, type = "o",
                    pch = 20, cex = 0.5, lwd = 1.5,
                    main = NULL, xlab = NULL, ylab = NULL,
                    ...){
  # extract
  t <- ATE_object$t
  surv_1 <- ATE_object$surv_ref
  surv_2 <- ATE_object$surv_trt
  ref <- ATE_object$ref
  ATE <- ATE_object$ATE
  trt_values <- ATE_object$trt_values
  
  # if ATE is selected
  if (!is.vector(mode)) stop("Mode provided is not valid")
  
  # Upper and lower quantiles
  lwr <- apply(ATE, 2, quantile, probs=.025)
  upr <- apply(ATE, 2, quantile, probs=.975)
  
  lwr_1 <- apply(surv_1, 2, quantile, probs=.025)
  upr_1 <- apply(surv_1, 2, quantile, probs=.975)
  
  lwr_2 <- apply(surv_2, 2, quantile, probs=.025)
  upr_2 <- apply(surv_2, 2, quantile, probs=.975)
  
  if (length(mode) == 2) {
    if (all(mode %in% c(0, 1))) { # if two treatment are selected

      # set necessary parameters if not specified
      if (is.null(xlim)) xlim = c(0, max(t))
      if (is.null(ylim)) ylim = c(0, 1)
      if (is.null(xlab)) xlab = "Time"
      if (is.null(ylab)) ylab = "Marginal Survival Probability" 
      if (is.null(main)) main = paste0("Marginal Survival Curves for Both Treatments")
      
      plot(t, colMeans(surv_1), pch = pch, 
           xlim = xlim, ylim = ylim,
           type = type, cex = cex, col = col_0,
           xlab = xlab, ylab = ylab, main = main,
           ...)
      segments(x0 = t, y0 = lwr_1,
               x1 = t, y1 = upr_1,
               lwd = lwd, col = col_CI_0)
      lines(t, colMeans(surv_2), pch = pch, 
            xlim = xlim, ylim = ylim,
            type = type, cex = cex, col = col_1,
            ...)
      segments(x0 = t, y0 = lwr_2,
               x1 = t, y1 = upr_2,
               lwd = lwd, col = col_CI_1)
      legend("topright", 
             legend = c(paste0(trt_values[1]), paste0(trt_values[2])),
             lty = c(1,1),
             col = c(col_0, col_1))
    } else {
      stop("Mode provided is not valid.")
    }
  } else if (length(mode) == 1) {
    if (mode == "ATE") {
      # set necessary parameters if not specified
      if (is.null(xlim)) xlim = c(0, max(t))
      if (is.null(ylim)) ylim = c(-1, 1)
      if (is.null(main)) main = "ATE Estimated Over Time"
      if (is.null(xlab)) xlab = "Time"
      if (is.null(ylab)) ylab = "ATE Estimated" 
      
      message(paste0("The current reference is ", ref))
      
      plot(t, colMeans(ATE), pch = pch, col = col_ATE, 
           xlim = xlim, ylim=ylim,
           type = type, cex = cex,
           xlab = xlab, ylab = ylab, main = main,
           ...)
      segments(x0 = t, y0 = lwr,
               x1 = t, y1 = upr,
               lwd = lwd, col = col_CI_ATE)
    } else if (mode == 0){
      if (is.null(main)) main = paste0("Marginal Survival Curve for ", trt_values[1])
      plot(t, colMeans(surv_1), pch = pch, 
           xlim = xlim, ylim = ylim,
           type = type, cex = cex,  col = col_0,
           xlab = xlab, ylab = ylab, main = main)
      segments(x0 = t, y0 = lwr_1,
               x1 = t, y1 = upr_1,
               lwd = lwd, col = col_CI_0)
    } else if (mode == 1){
      if (is.null(main)) main = paste0("Marginal Survival Curve for ", trt_values[2])
      plot(t, colMeans(surv_2), pch = pch, 
           xlim = xlim, ylim = ylim,
           type = type, cex = cex,  col = col_1,
           xlab = xlab, ylab = ylab, main = main)
      segments(x0 = t, y0 = lwr_2,
               x1 = t, y1 = upr_2,
               lwd = lwd, col = col_CI_1)
    } else {
      stop("Mode provided is not valid.")
    } 
  } else {
    stop("Mode provided is not valid.")
  }
}
