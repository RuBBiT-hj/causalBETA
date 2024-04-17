#' causaBETA pipeline
#' 
#' Conduct the causalBETA pipeline in a row: bayesian piece-wise exponential modeling and the following average treatment effect estimation
#' @param d data, a data frame in survival format. Categorical variables should be
#' transformed into dummy variables
#' @param reg_formula a formula object that specifies the formula for the poisson regression.
#' This also decides the formula will be used the function to check positivity overlap.
#' @param A a character variable that specifies the name of the treatment
#' @param model a character variable that tells the stan model used to implement the Bayesian piece-wise exponential model, 
#' default is "AR1" and the other option is "independent"
#' @param sigma a numeric variable as the user-defined standard deviation for beta coefficients prior, the default is 3
#' @param num_partitions a numeric variable as the number of partitions of the study time, the default is 100
#' @param warmup a numeric variable as the number of warmup in MCMC, the default is 1000
#' @param post_iter a numeric variable as the number of iterations to draw from the posterior, the default is 1000
#' @param ref the reference value of the treatment, so it should be one of the treatment values
#' @param t optional; a numeric vector of time points at which users want to compute marginal
#' survival probabilities
#' @param V the number of prediction for each posterior draw; the default is 1000
#' @param func the function to calculate the statistics in interest; the default is the posterior survival
#' probability at each value of t, `function(x){mean(x > t[i])}`
#' @param ... additional arguments required for func
#' 
#' @details
#' The function is the wrapper function of the main function `bayeshaz()` and `bayesgcomp` to complete the computation
#' part of the analysis. For details, check the documentation for these two functions by `help(bayeshaz)` and `help(bayesgcomp)`.
#' 
#' @returns
#' This function returns a list of one object of class `bayesahz` and `ATE` storing the modeling results
#' and the ATE result respectively. For details regarding these two objects, check `help(bayeshaz)` and `help(bayesgcomp)`.
#' 
#' It also plots the result of ATE by calling `plot()` on the ATE object.
#'
#' @examples
#' # example demo
## usethis namespace: start
#' @import survival
#' @importFrom mets rpch
#' @importFrom LaplacesDemon rdirichlet
## usethis namespace: end
#' @export

bayespipeline <- function(d, reg_formula, A, model = "AR1", sigma = 3, 
                          num_partitions=100, warmup=1000, post_iter=1000,
                          ref, t = NULL, V = 1000,
                          func = function(x){mean(x > t[i])}, ...){
  bayeshaz_object = bayeshaz(d = d,
                             reg_formula = reg_formula, 
                             A = A,
                             model = model,
                             sigma = sigma,
                             num_partitions = num_partitions,
                             warmup = warmup,
                             post_iter = post_iter)
  ATE_object = bayesgcomp(bayeshaz_object,
                          ref = ref,
                          t = t,
                          V = V,
                          func = func, ...)
  plot(ATE_object)
  return(list(bayeshaz_object, ATE_object))
}
