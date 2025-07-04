#' Piece-wise Constant Hazard Distribution Prediction
#' 
#' Predict quantities of interest using the posterior draws from the piecewise exponential model
#' @param object an object of the class `bayeshaz` created by the `bayeshaz()` function
#' @param x a data frame for the individual(s) used for prediction, with the variables specified in the model
#' @param n a numeric value as the number of predictions for each posterior draw; the default is 1000
#' @param func a function used to obtain the quantities of interest; 
#' a common choice is `mean` to get the posterior distribution of the expectation of survival time (more details in the demo).
#' @param ... additional arguments affecting the predictions produced.
#' 
#' @details
#' The data frame used for prediction could be for a single individual or multiple individuals.
#' It is required that the number and the order of the variables matches the model for generating the
#' `bayeshaz` object exactly.
#' 
#' @returns
#' This function returns a matrix that contains the predicted quantities of interest.
#' Each column represents an individual, and each row represents a posterior draw.
#' 
#' This function use bootstrapping method to get the quantities of interest.
#' It uses the parameters from each posterior draw to make a set of predictions for the individual(s) provided.
#' Within each posterior draw, the function eventually obtains the quantity of interest based on `n` predictions and the `func` given.
#' 
#' @examples
#' \dontrun{
#' # example demo
#' ## Continued from ?bayeshaz
#' predict(post_draws_ar1_adj,
#'         x = data[1, c("A", "age", "karno", 
#'                       "celltypesquamous",
#'                       "celltypesmallcell",
#'                       "celltypeadeno")],
#'         func = mean)
#' }
## usethis namespace: start
#' @import survival
#' @importFrom mets rpch
## usethis namespace: end
#' @export


predict.bayeshaz = function(object, x, n = 1000, func, ...){
  
  beta_draws = do.call(rbind, object$beta_draws)
  haz_draws = do.call(rbind, object$haz_draws)
  partition = object$partition
  
  reg_formula = object$formula
  
  # one hot encoding for the x

  # # model matrix (one hot encoding)
  x =  model.matrix(formula(
    paste(as.character(reg_formula[c(1,3)]), collapse = "")
  ), data = x)
  x = data.frame(x[, -1, drop = F])
  
  
  
  # check that the dimension matches
  if (is.null(dim(beta_draws))) { # only trt
    if (!is.null(dim(x))) stop("The number of variables in the data frame doesn't match the model")
  } else {
    if (dim(x)[2] != dim(beta_draws)[2]) stop("The number of variables in the data frame doesn't match the model")
  }
  
  post_iter = nrow(haz_draws)
  
  res = apply(x, MARGIN = 1, function(y) {
    res_single = predict_haz(x = y,
                             beta_draws = beta_draws,
                             haz_draws = haz_draws,
                             partition = partition,
                             n = n, func = func)
    return(res_single)
  })
  
  return(res)
}
