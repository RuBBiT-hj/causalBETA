#' Bayesian Piece-Wise Exponential Model
#' 
#' Perform a bayesian piece-wise exponential model on the given survival data, 
#' and this function implements it by an equivalent poisson regression in MCMC.
#' 
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
#' 
#' @details
#' A typical model has the form `Surv(time, outcome) ~ covariates`. The function will capture
#' the outcome and covariates based on this formula object.
#' 
#' The Bayesian piece-wise exponential model uses normal prior for baseline hazard rate,
#' beta coefficients, and the error term for baseline hazard rate to form the first-order Gaussian process.
#' It uses a beta prior to generate the correlation from -1 to 1 for hazard rate.
#' For more details, users can check the full model listed in the reference.
#' 
#' Under the `AR1` model, user can specify the standard deviation for the normal prior of the beta coefficients.
#' The default is 3, and only values between 0 and 3 are accepted since 3 is already a relatively weak prior. 
#' 
#' @return It returns an object of class `bayeshaz` that contains the information about the data, model, etc.
#' This serves as the basis for the extended functions in this package.
#' 
#' An object of class `bayeshaz` is a list containing at least the following components:
#' 
#' * `data`, a data frame for the data
#' * `formula`, the regression formula
#' * `treatment`, the name of the treatment variable
#' * `covariates`, the name of the covariate(s)
#' * `time`, the name of the time variable
#' * `outcome`, the name of the outcome variable
#' * `model`, the type of the model used
#' * `sigma`, the sigma specified
#' * `partition`, the partition vector
#' * `midpoint`, the midpoints of intervals
#' * `haz_draws`, the baseline hazard rate from each posterior draws
#' * `beta_draws`, the beta coefficients estimated from each posterior draws
#' 
#' @references
#' Oganisian, Arman, Anthony Girard, Jon A. Steingrimsson, and Patience Moyo.
#' "A Bayesian Framework for Causal Analysis of Recurrent Events in Presence of Immortal Risk."
#' arXiv preprint arXiv:2304.03247 (2023).
#' 
#' @examples
#' # example demo
#' df_veteran <- survival::veteran
#' df_veteran$trt <- ifelse(df_veteran$trt == 2, 1, 0)
#' post_draws <- bayeshaz(d = df_veteran,
#'    reg_formula = Surv(time, status) ~ trt,
#'    A = 'trt')
## usethis namespace: start
#' @import cmdstanr
#' @importFrom survival survSplit
## usethis namespace: end
#' @export

bayeshaz = function(d, reg_formula, A, model = "AR1", sigma = 3, 
                    num_partitions=100, warmup=1000, post_iter=1000){
  ## dependency checkings
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package \"cmdstanr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(
      "Package \"survival\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # the address of the stan files
  path_stan <- paste0(.libPaths(), "/causalBETA/data/")
  
  # check for sigma value
  if (sigma > 3 | sigma <= 0) {
    warning("sigma must be not less than 0 and less than 3. Forced to be 3\n")
    ## force sigma to be 3
    sigma <- 3
  }
  
  
  ## user-specified intervention variable
  trt_names = A
  
  outcome_char = gsub(" ", "", as.character(reg_formula[2]))
  
  ### collect user inputs 
  
  y_name =substr(outcome_char, 6, gregexpr(",", outcome_char)[[1]][1] - 1  )
  delta_name = substr(outcome_char, gregexpr(",", outcome_char)[[1]][1]+1 ,  gregexpr(")", outcome_char)[[1]][1]-1 )
  
  covar_char = gsub(" ", "", as.character(reg_formula[3]))
  
  # consider * situations only
  covar_names = as.character(attr(terms(reg_formula), "variables"))
  covariates = covar_names[3:length(covar_names)]
  covariates = covariates[covariates != trt_names]
  
  y = d[, y_name]
  delta = d[,delta_name]
  
  ## partition time interval
  partition = seq(0, max(y)+.01,length.out=num_partitions)
  
  #outcome = paste0("Surv(", y_name,", ", delta_name,") ~ ")
  #reg_formula = as.formula(paste0(outcome, paste0(covar_names, collapse = "+"), "+", paste0(trt_names, collapse = "+") ) )
  
  ## create long-form data set
  dsplit = survival::survSplit(data = d, formula = reg_formula, cut = partition, id='id')
  
  dsplit$offset = dsplit[,y_name] - dsplit$tstart ## amount of survival time in each interval
  dsplit$interval = as.factor(dsplit$tstart) ## factor interval
  dsplit$interval_num = as.numeric(as.factor(dsplit$tstart)) ## interval number
  
  ## create covariate model matrix
  xmat = model.matrix(data=dsplit, object = as.formula(paste0(" ~ -1 + ", covar_char  )) )
  
  ## create list of data to pass to Stan model 
  
  ## use different stan files for different model input
  
  if (model == "independent"){ # independent instead of being auto-regressive
    dlist = list(N=nrow(dsplit),
                 P = ncol(xmat),
                 n_pieces = length( unique(dsplit$interval_num) ),
                 delta = dsplit[, delta_name], 
                 off_set  = dsplit$offset, 
                 interval_num = dsplit$interval_num,
                 xmat = xmat)
    mod = cmdstan_model(paste0(path_stan, "hazard_mod_v1.stan"))
  } else if (model == "AR1"){ # a different variance for beta coefficients

    dlist = list(N=nrow(dsplit),
                 P = ncol(xmat),
                 n_pieces = length( unique(dsplit$interval_num) ),
                 delta = dsplit[, delta_name], 
                 off_set  = dsplit$offset, 
                 interval_num = dsplit$interval_num,
                 xmat = xmat,
                 sigma_beta = sigma)
    mod = cmdstan_model(paste0(path_stan, "hazard_mod_v2.stan"))
  } else { # the model input is not correct
    stop("The model input is not valid")
  }
  
  
  res = mod$sample(data= dlist,
                   chains = 1,  iter_warmup = warmup, iter_sampling = post_iter)
  
  haz_draws = exp(res$draws("log_haz", format = 'matrix') )
  beta_draws = res$draws("beta", format = 'matrix')
  
  xv = (partition[-1] - .5*mean(diff(partition)) ) ## midpoint of each interval
  
  draws = create_bayeshaz(data = d, formula = reg_formula, treatment = A,
                          covariates = covariates,
                          time = y_name, outcome = delta_name,
                          model = model, sigma = sigma, 
                          partition = partition, midpoint = xv,
                          haz_draws = haz_draws, beta_draws=beta_draws)
  
  return(draws)
}

