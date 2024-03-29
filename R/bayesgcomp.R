#' Average Treatment Effect Estimation
#' 
#' Takes as input results from bayeshaz() and runs g-computation and returns posteriors draws of the average difference in survival probabilities between the two treatments.
#' By default, survival probabilities and differences are calculated at the midpoints of each of the partitions specified in bayeshaz().
#' @param bayeshaz_object an object of the class `bayeshaz` created by the `bayeshaz()` function
#' @param ref the reference value of the treatment, so it should be one of the treatment values
#' @param t optional; a numeric vector of time points at which users want to compute marginal
#' survival probabilities
#' @param B the number of prediction for each posterior draw; the default is 1000
#' 
#' @details
#' Takes as input a `bayeshaz` object generated by bayeshaz() and runs g-computation and returns posteriors draws of the average difference in survival probabilities between the two treatments.
#' Averaging is done over a Bayesian bootstrap draw of the baseline confounder distribution.
#' 
#' It is important to specify the correct reference level for this function and also other ATE-related
#' ones. Users can specify or change reference by `relevel` function or the `ref` argument in these
#' ATE-related functions. The reference value needs to match the type of treatment variable, and it
#' needs to belong to one of the treatment values.
#' 
#' In addition, users can use argument `t` to specify the time points where they want to estimate 
#' the posterior draws of the average difference of the survival probabilities. If not specified,
#' the function will use the default time points - the midpoints stored in the `bayeshaz` object.
#' 
#' These ATE-related functions only considers binary treatments so it is not applicable for more than two treatment values, or the multiple levels need to be coerced into a
#' binary variable before creating the `bayeshaz` object.
#' 
#' @returns
#' This function returns an object of class `ATE` that stores the following information:
#' 
#' * `surv_ref`, an `mcmc` object storing the marginal survival probability for the reference
#' * `surv_trt`, an `mcmc` object storing the marginal survival probability for the treatment
#' * `ref`, the value of the reference treatment
#' * `trt_values`, the possible values of the treatment
#' * `ATE`, the difference between the marginal survival probability of the treatment and the reference
#' 
#' @references
#' Ji, Han, and Oganisian, Arman. 2023. 
#' "causalBETA: An R Package for Bayesian Semiparametric Causal Inference with 
#' Event-Time Outcomes. *arXiv:2310.12358 \[Stat\]*, October. 
#' \url{http://arxiv.org/abs/2310.12358}.
#' 
#' @examples
#' # example demo
#' ## Continued from ?bayeshaz
#' gcomp_res = bayesgcomp(post_draws_ar1_adj, ## bayeshaz output 
#'                        ref = 0, ## treatment reference group
#'                        B = 1000) ## monte carlo iterations in g-comp
## usethis namespace: start
#' @import survival
#' @importFrom mets rpch
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom coda mcmc
## usethis namespace: end
#' @export


bayesgcomp = function(bayeshaz_object, ref, t = NULL, B = 1000){
  
  # calculate the survival probability for each subject at each time
  
  # extract values
  d <- bayeshaz_object$data
  trt <- bayeshaz_object$treatment
  trt_vector <- d[, trt]
  
  beta_draws <- bayeshaz_object$beta_draws
  haz_draws <- bayeshaz_object$haz_draws
  partition <- bayeshaz_object$partition
  covariates <- bayeshaz_object$covariates
  n_draws = nrow(beta_draws)
  
  # check if users provided valid t values
  if (is.null(t)) {
    t <- bayeshaz_object$midpoint # default
  } else {
    # numeric and vector
    if ( !(is.numeric(t) & is.vector(t)) ) {
      stop("t provided should be a numeric vector")
    }
    # not exceed the limit
    if (max(t) > max(partition)) {
      stop("t exceeds the limit of partitions")
    }
  }
  
  # the possible values of the treatment
  trt_values <- unique(trt_vector)
  
  # check that the treatment is binary
  
  if (length(trt_values) !=2 ) stop("Treatment variable needs to be binary")
  
  # check the ref input
  # if the value is one of the treatment variable
  if (!(ref %in% trt_values)) stop("Reference provided is not a valid treatment value")

  # the number of subjects
  n_subject <- nrow(d)
  
  # make two datasets with two treatment values - ref and treatment
  # variables to include
  var_included <- c(trt, covariates)
  
  d_1 <- subset(d, select = var_included)
  d_1[, trt] <- ref
  
  d_2 <- subset(d, select = var_included)
  d_2[, trt] <- trt_values[trt_values != ref]
  
  # the survival time for the first and the second treatment value
  #   then calculate the difference
  # each matrix is for one subject, and n_draws x n_partitions
  cat(paste0('Running g-comp procedure...','\n'))
  cat(paste0('Calculating Posterior Draws of Survival Probability under reference value...',Sys.time(),'\n'))
  surv_prob_1 <- lapply(1:n_subject, function(x){
    # call predict.haz as a helper function
    all_surv_time = predict.haz(d_1[x, ], beta_draws, haz_draws, partition, B, func = list)
    # the number of the list - the number of posterior draws
    # the length of each list - the number of predictions (B)
    surv_prob <- matrix(nrow = length(all_surv_time),
                        ncol = length(t))
    # For each t, we calculate the proportion
    for (i in 1:length(t)){
      surv_prob[ ,i] <- sapply(all_surv_time, function(x) mean(x > t[i]))
    }
    return(surv_prob)
  })
  
  cat(paste0('Calculating Posterior Draws of Survival Probability under treatment...',Sys.time(),'\n'))
  surv_prob_2 <- lapply(1:n_subject, function(x){
    # call predict.haz as a helper function
    all_surv_time = predict.haz(d_2[x, ], beta_draws, haz_draws, partition, B, func = list)
    # the number of the list - the number of posterior draws
    # the length of each list - the number of predictions (B)
    surv_prob <- matrix(nrow = length(all_surv_time),
                        ncol = length(t))
    # For each t, we calculate the proportion
    for (i in 1:length(t)){
      surv_prob[ ,i] <- sapply(all_surv_time, function(x) mean(x > t[i]))
    }
    return(surv_prob)
  })
  
  cat(paste0('Calculating Posterior Draws of Surv. Prob. Difference...',Sys.time(),'\n'))
  surv_prob_diff <- lapply(1:n_subject, function(i){
    return(surv_prob_2[[i]] - surv_prob_1[[i]])
  })
  
  # calculate the posterior surv_prob for each treatment and posterior ATE
  surv_1_post <- matrix(nrow = n_draws,
                        ncol = length(t))
  surv_2_post <- matrix(nrow = n_draws,
                        ncol = length(t))
  ATE <- matrix(nrow = n_draws,
                ncol = length(t))
  
  # alpha vector for rdirichlet
  alpha <- rep(1, n_subject)
  
  cat(paste0('Computing Bayesian Bootstrap Weighted Average...',Sys.time(),'\n'))
  for (i in 1:n_draws){

    # draw Dirichlet for weighted sum
    weights_dir <- matrix(LaplacesDemon::rdirichlet(1, alpha = alpha), ncol = 1)
    
    surv_1_matrix <- sapply(surv_prob_1, function(x) x[i,])
    surv_2_matrix <- sapply(surv_prob_2, function(x) x[i,])
    diff_matrix <- sapply(surv_prob_diff, function(x) x[i,]) # n_parition x n_subject
    
    surv_1_post[i, ] <- surv_1_matrix %*% weights_dir
    surv_2_post[i, ] <- surv_2_matrix %*% weights_dir
    ATE[i, ] <- diff_matrix %*% weights_dir
  }
  
  surv_1_post = mcmc(surv_1_post, start = 1, end = n_draws, thin = 1)
  surv_2_post = mcmc(surv_2_post, start = 1, end = n_draws, thin = 1)
  ATE = mcmc(ATE, start = 1, end = n_draws, thin = 1)
  
  ATE_object = create_ATE(surv_ref = surv_1_post, surv_trt = surv_2_post, 
                          trt_values = c(ref, trt_values[trt_values != ref]), 
                          ref = ref, ATE = ATE, t = t)
  
  cat(paste0('g-comp complete...',Sys.time(),'\n'))
  
  return(ATE_object)
}
