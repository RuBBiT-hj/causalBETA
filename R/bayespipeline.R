#' causaBETA pipeline
#' 
#' Conduct the causalBETA pipeline in a row: bayesian piece-wise exponential modeling and the following average treatment effect estimation
#' @param data data, a data frame in survival format. Categorical variables should be
#' transformed into dummy variables
#' @param reg_formula a formula object that specifies the formula for the poisson regression.
#' This also decides the formula will be used the function to check positivity overlap.
#' @param A a character variable that specifies the name of the treatment
#' @param model a character variable that tells the stan model used to implement the Bayesian piece-wise exponential model, 
#' default is "AR1" and the other option is "independent"
#' @param priorSD a numeric variable as the user-defined standard deviation for beta coefficients prior, the default is 3
#' @param num_partition a numeric variable as the number of intervals in the partition, the default is 100
#' @param warmup a numeric variable as the number of warmup in MCMC, the default is 1000
#' @param post_iter a numeric variable as the number of iterations to draw from the posterior, the default is 1000
#' @param chains the number of chains for sampling, the default is 1
#' @param ref the reference value of the treatment, so it should be one of the treatment values
#' @param B the number of prediction for each posterior draw; the default is 1000
#' @param estimand the statistics in interest; the default is the posterior survival
#' probability at each value of t, `prob`. The other options are the median survival time, `median`, and
#' the Restricted mean survival time, `rmean`.
#' @param t optional; a numeric vector of time points at which users want to compute marginal
#' survival probabilities
#' @param threshold optional; the threshold used for Restricted mean survival time.
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
#' \dontrun{
#' # example demo
#' data = survival :: veteran
#' data$A = 1*(data$trt==2)
#' ## rename variables
#' var_names = colnames(data)
#' colnames(data)[var_names=='status'] = 'delta'
#' colnames(data)[var_names=='time'] = 'y'
#' results <- bayespipeline(data,
#'   reg_formula = Surv(y, delta) ~ A,
#'   A = 'A',
#'   ref = 0
#' )}
## usethis namespace: start
#' @import survival
#' @import cmdstanr
#' @import coda
#' @importFrom mets rpch
#' @importFrom LaplacesDemon rdirichlet
## usethis namespace: end
#' @export

bayespipeline <- function(data, reg_formula, A, model = "AR1", priorSD = 3, 
                          num_partition=100, warmup=1000, post_iter=1000,
                          chains = 1,
                          ref, B = 1000,
                          estimand = "prob", 
                          t = NULL, 
                          threshold = NULL, ...){
  bayeshaz_object = bayeshaz(data = data,
                             reg_formula = reg_formula, 
                             A = A,
                             model = model,
                             priorSD = priorSD,
                             num_partition = num_partition,
                             warmup = warmup,
                             post_iter = post_iter,
                             chains = chains)
  
  ATE_object = bayesgcomp(bayeshaz_object,
                          ref = ref,
                          B = B,
                          estimand = estimand,
                          t = t, 
                          threshold = threshold,
                          ...)
  
  plot(ATE_object)
  return(list(bayeshaz_object, ATE_object))
}

