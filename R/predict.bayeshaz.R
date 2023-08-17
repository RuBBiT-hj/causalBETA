#' Piece-wise Constant Hazard Distribution Prediction
#' 
#' Predict quantities of interest using the posterior draws from the piecewise exponential model
#' @param bayeshaz_object an object of the class `bayeshaz` created by the `bayeshaz()` function
#' @param x the individual used for prediction, with the variables specified in the model
#' @param n the number of prediction for each posterior draw; the default is 1000
#' @param func the function specified to get the quantities of interest; a common choice is `mean` to get the posterior distribution of the expectation of survival time
#' @examples
#' # example demo
## usethis namespace: start
#' @import survival
#' @importFrom mets rpch
## usethis namespace: end
#' @export


predict.bayeshaz = function(bayeshaz_object, x, n=1000, func){
  
  beta_draws = bayeshaz_object$beta_draws
  haz_draws = bayeshaz_object$haz_draws
  partition = bayeshaz_object$partition
  
  post_iter = nrow(haz_draws)
  res = numeric(length = post_iter)
  
  for(j in 1:post_iter){
    hazv = haz_draws[j, ]*exp( sum( beta_draws[j,] * x  ) )
    simv = mets::rpch(n, hazv, partition)
    res[j] = func(simv) # get the result from the function specified
  }
  
  return(res)
}
