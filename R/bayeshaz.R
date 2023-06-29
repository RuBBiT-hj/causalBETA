#' Bayesian Piece-Wise Exponential Model
#' 
#' Perform a bayesian piece-wise exponential model on the given survival data
#' Use an equivalent poisson regression in MCMC
#' @param d The data in survival format, i.e. processed by Surv function
#' @param reg_formula The formula for the poisson regression
#' @param A The name of the treatment
#' @param model The stan model used, default is "AR1", other options are "independent" and "beta"
#' @param sigma The user-defined sd for odds ratio prior, the default is 3, the same as the default model
#' @param num_intervals The number of intervals to partition the study time, the default is 100
#' @param warmup The number of warmup in MCMC, the default is 1000
#' @param post_iter The number of iterations to draw from the posterior, the default is 1000
#' @examples
#' # example demo
## usethis namespace: start
#' @import cmdstanr
#' @importFrom survival survSplit
## usethis namespace: end
#' @export
bayeshaz = function(d, reg_formula, A, model = "AR1", sigma = 3, 
                    num_intervals=100, warmup=1000, post_iter=1000){
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
  path_stan <- paste0(.libPaths(), "/BayesSurvival/data/")
  
  ## user-specified intervention variable
  trt_names = A
  
  outcome_char = gsub(" ", "", as.character(reg_formula[2]))
  
  ### collect user inputs 
  
  y_name =substr(outcome_char, 6, gregexpr(",", outcome_char)[[1]][1] - 1  )
  delta_name = substr(outcome_char, gregexpr(",", outcome_char)[[1]][1]+1 ,  gregexpr(")", outcome_char)[[1]][1]-1 )
  
  covar_char = gsub(" ", "", as.character(reg_formula[3]))
  
  # can use strsplit(covar_char, "[\\+]+")
  # covar_names = c('age', 'ph.ecog', 'sex') ## need to get this from covar_char somehow.
  covar_names = strsplit(covar_char, "[\\+]+")[[1]]
  
  y = d[, y_name]
  delta = d[,delta_name]
  
  ## partition time interval
  partition = seq(0, max(y)+.01,length.out=num_intervals)
  
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
                 offset  = dsplit$offset, 
                 interval_num = dsplit$interval_num,
                 xmat = xmat)
    mod = cmdstan_model(paste0(path_stan, "hazard_mod_v1.stan"))
  } else if (model == "beta"){ # a different variance for beta coefficients
    if (sigma > 3 | sigma <= 0) {
      warning("B must be not less than 0 and less than 3. Forced to be 3")
      ## force B to be 3
      sigma <- 3
    }
    dlist = list(N=nrow(dsplit),
                 P = ncol(xmat),
                 n_pieces = length( unique(dsplit$interval_num) ),
                 delta = dsplit[, delta_name], 
                 offset  = dsplit$offset, 
                 interval_num = dsplit$interval_num,
                 xmat = xmat,
                 sigma = sigma)
    mod = cmdstan_model(paste0(path_stan, "hazard_mod_v2.stan"))
    
  } else if (model == "AR1"){ # the original version
    dlist = list(N=nrow(dsplit),
                 P = ncol(xmat),
                 n_pieces = length( unique(dsplit$interval_num) ),
                 delta = dsplit[, delta_name], 
                 offset  = dsplit$offset, 
                 interval_num = dsplit$interval_num,
                 xmat = xmat)
    
    mod = cmdstan_model(paste0(path_stan, "hazard_mod.stan"))
  } else { # the model input is not correct
    stop("The model input is not valid")
  }
  
  
  res = mod$sample(data= dlist,
                   chains = 1,  iter_warmup = warmup, iter_sampling = post_iter)
  
  haz_draws = exp(res$draws("log_haz", format = 'matrix') )
  beta_draws = res$draws("beta", format = 'matrix')
  
  xv = (partition[-1] - .5*mean(diff(partition)) ) ## midpoint of each interval
  
  draws = list(haz_draws = haz_draws, beta_draws=beta_draws,
               xv = xv, partition = partition )
  
  return(draws)
}

