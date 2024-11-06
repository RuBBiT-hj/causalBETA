#' Construct the Object for Bayesian Piece-wise Exponential Model
#' 
#' Construct a `bayeshaz` object to store the information of the bayesian piece-wise
#' exponential model and the posterior parameters. This constructor function is hidden
#' from the users
#' 
#' @param data data, a data frame in survival format
#' @param formula a formula object that specifies the formula for the poisson regression.
#' This also decides the formula will be used the function to check positivity overlap.
#' @param treatment a character variable that specifies the name of the treatment
#' @param covariates a character variable storing the name of the covariate(s)
#' @param time a character variables that specifies the name of the time variable
#' @param outcome a character variables that specifies the name of the outcome variable
#' @param model a character variable recording which mode of model was used
#' @param priorSD a numeric variable as the user-defined standard deviation for beta coefficients prior
#' @param chains the number of chains for sampling, the default is 1
#' @param partition a numeric vector for the partition
#' @param midpoint a numeric vector for the midpoints of intervals
#' @param haz_draws a `mcmc.list` for the baseline hazard rates from posterior draws
#' @param beta_draws a `mcmc.list` for the beta coefficients estimated from posterior draws
#' 
#' @return It returns an object of class `bayeshaz`

create_bayeshaz <- function(data, formula, treatment, covariates, time, outcome,
                            model, priorSD, chains, partition,
                            midpoint, haz_draws, beta_draws) {

  my_object <- structure(list(
    data = data, formula = formula, treatment = treatment, covariates = covariates,
    time = time, outcome = outcome,
    model = model, priorSD = priorSD, 
    chains = chains,
    partition = partition, midpoint = midpoint,
    haz_draws = haz_draws, beta_draws = beta_draws
  ), class = "bayeshaz")
  return(my_object)
}