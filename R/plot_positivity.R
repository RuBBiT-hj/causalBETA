#' Positivity Assumption Checking
#' 
#' Plot the histogram of the proprensity score of treatment assignment conditional on the observed treatment. 
#' This plot can be used to check the positivity assumption for casual inference.
#' If the assumption is met, the mirrored histogram should look similar.
#' If there is an obvious shift in the distribution of histograms, the assumption is likely violated.
#' 
#' @param d the data frame in survival format, i.e. processed by Surv function
#' @param formula the formula for the logistic regression for propensity score (a formula object in this version)
#' @param treatment the treatment variable name

#' @examples
#' # example demo
## usethis namespace: start
## usethis namespace: end
#' @export

plot_positivity = function(d, formula, treatment){
  # Fit the logistic regression
  pos_check <- glm(formula, family = binomial(link = "logit"),
                   data = d)
  d_trt <- d[, treatment]
  trt_1 <- unique(d_trt)[1]
  trt_2 <- unique(d_trt)[2]
  
  
  score_1 <- predict(pos_check, newdata = d[d_trt == trt_1, ], type = "response")
  score_2 <- predict(pos_check, newdata = d[d_trt == trt_2, ], type = "response")
  
  # Get parameters for histogram
  h_1 <- hist(score_1, breaks = "Scott", probability = T)
  h_2 <- hist(score_2, breaks = "Scott", probability = T)
  
  x_min = min(c(score_1, score_2))*0.98
  x_max = min(1, max(c(score_1, score_2))*1.02) # not exceed 1
  y_max = 1.02 * max(c(h_1$density, h_2$density))
  
  # Plot
  par(mfrow=c(2,1), mar=c(0,5,3,3))
  hist(score_1 , main="Positivity Assumption Checking" , ylab="Density", xlab="", 
       ylim=c(0, y_max), xlim = c(x_min, x_max), 
       xaxt="n", las=1 , col="slateblue1", breaks="Scott", probability = T)
  par(mar=c(5,5,0,3))
  hist(score_2 , main="" , ylab="Density", xlab="Propensity Score", 
       ylim=c(y_max, 0),  xlim = c(x_min, x_max),
       las=1 , col="tomato3"  , breaks = "Scott", probability = T)
  # reset
  par(mfrow=c(1,1),
      mar = c(5.1, 4.1, 4.1, 2.1),
      mgp = c(3,1,0))
  
}

