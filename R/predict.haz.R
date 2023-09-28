#' Prediction helper function
#' 
#' Predict quantities of interest using the posterior draws from the piece-wise exponential model
#' @param x the individual used for prediction, with the variables specified in the model
#' @param beta_draws the beta coefficients from the posterior draws, with dimension number of posterior draws * the number of variables 
#' @param haz_draws the baseline hazard rates from the posterior draws, with dimension number of posterior draws * the number of intervals
#' @param partition the time values for each partition
#' @param n the number of prediction for each posterior draw; the default is 1000
#' @param func the function specified to get the quantities of interest; a common choice is `mean` to get the posterior distribution of the expectation of survival time
#' @return The posterior predictive quantity in interest


predict.haz = function(x, beta_draws, haz_draws, partition, n = 1000, func){
  
  post_iter = nrow(haz_draws)
  res = numeric(length = post_iter)
  
  for(j in 1:post_iter){
    hazv = haz_draws[j, ]*exp( sum( beta_draws[j,] * x  ) )
    simv = mets::rpch(n, hazv, partition)
    res[j] = func(simv) # get the result from the function specified
  }
  
  return(res)
}
