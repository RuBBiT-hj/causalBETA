#' Average Treatment Effect Estimation
#' 
#' Estimate the average treatment effect (ATE) based on the predicted survival probability
#' @param bayeshaz_object an object of the class `bayeshaz` created by the `bayeshaz()` function
#' @param ref the reference value of the treatment, so it should be one of the treatment values
#' @param n the number of prediction for each posterior draw; the default is 1000
#' 
#' @details
#' This function estimates the average treatment effect (ATE) for the data set used to generate
#' the `bayeshaz` object. It calculates the marginal survival probability for all individuals
#' while assuming they all receive the same treatment, like all in control or all in treatment group.
#' The ATE estimated would be the weighted sum of the difference between the survival probabilities,
#' and we use a Dirichlet prior for the weights.
#' 
#' It is important to specify the correct reference level for this function and also other ATE-related
#' ones. Users can specify or change reference by `relevel` function or the `ref` argument in these
#' ATE-related functions. The reference value needs to match the type of treatment variable, and it
#' needs to belong to one of the treatment values.
#' 
#' These ATE-related functions only considers treatment as a binary variable (i.e. control and treatment),
#' so it is not applicable for more than two treatment values, or they needs to be coerced into a
#' binary variable before creating the `bayeshaz` object.
#' 
#' @returns
#' This function returns an object of class `ATE` that stores the following information:
#' 
#' * `surv_ref`, the marginal survival probability for the reference
#' * `surv_trt`, the marginal survival probability for the treatment
#' * `ref`, the value of the reference treatment
#' * `trt_values`, the possible values of the treatment
#' * `ATE`, the difference between the marginal survival probability of the treatment and the reference
#' @examples
#' # example demo
## usethis namespace: start
#' @import survival
#' @importFrom mets rpch
#' @importFrom LaplacesDemon rdirichlet
## usethis namespace: end
#' @export


ATE_estimation = function(bayeshaz_object, ref, n = 1000){
  
  # calculate the survival probability for each subject at each time
  
  # extract values
  d <- bayeshaz_object$data
  trt <- bayeshaz_object$treatment
  trt_vector <- d[, trt]
  
  t <- bayeshaz_object$midpoint
  beta_draws <- bayeshaz_object$beta_draws
  haz_draws <- bayeshaz_object$haz_draws
  partition <- bayeshaz_object$partition
  covariates <- bayeshaz_object$covariates
  
  
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
  surv_prob_1 <- lapply(1:n_subject, function(x){
    # call predict.haz as a helper function
    all_surv_time = predict.haz(d_1[x, ], beta_draws, haz_draws, partition, n, func = list)
    # the number of the list - the number of posterior draws
    # the length of each list - the number of predictions (n)
    surv_prob <- matrix(nrow = length(all_surv_time),
                        ncol = length(t))
    # For each t, we calculate the proportion
    for (i in 1:length(t)){
      surv_prob[ ,i] <- sapply(all_surv_time, function(x) mean(x > t[i]))
    }
    return(surv_prob)
  })
  
  surv_prob_2 <- lapply(1:n_subject, function(x){
    # call predict.haz as a helper function
    all_surv_time = predict.haz(d_2[x, ], beta_draws, haz_draws, partition, n, func = list)
    # the number of the list - the number of posterior draws
    # the length of each list - the number of predictions (n)
    surv_prob <- matrix(nrow = length(all_surv_time),
                        ncol = length(t))
    # For each t, we calculate the proportion
    for (i in 1:length(t)){
      surv_prob[ ,i] <- sapply(all_surv_time, function(x) mean(x > t[i]))
    }
    return(surv_prob)
  })
  
  surv_prob_diff <- lapply(1:n_subject, function(i){
    return(surv_prob_2[[i]] - surv_prob_1[[i]])
  })
  
  # calculate the posterior surv_prob for each treatment and posterior ATE
  surv_1_post <- matrix(nrow = nrow(beta_draws),
                        ncol = length(t))
  surv_2_post <- matrix(nrow = nrow(beta_draws),
                        ncol = length(t))
  ATE <- matrix(nrow = nrow(beta_draws),
                ncol = length(t))
  
  # alpha vector for rdirichlet
  alpha <- rep(1, n_subject)
  
  for (i in 1:nrow(beta_draws)){
    
    cat(paste0('Running G-comp. iteration for each posterior draw...',i,'\n'))
    if( i %% 50 == 0 ){
      cat(paste0('Iteration ',i,'\n'))
    }
    
    # draw Dirichlet for weighted sum
    weights_dir <- matrix(LaplacesDemon::rdirichlet(1, alpha = alpha), ncol = 1)
    
    surv_1_matrix <- sapply(surv_prob_1, function(x) x[i,])
    surv_2_matrix <- sapply(surv_prob_2, function(x) x[i,])
    diff_matrix <- sapply(surv_prob_diff, function(x) x[i,]) # n_parition x n_subject
    
    surv_1_post[i, ] <- surv_1_matrix %*% weights_dir
    surv_2_post[i, ] <- surv_2_matrix %*% weights_dir
    ATE[i, ] <- diff_matrix %*% weights_dir
  }
  
  # construct the ATE object
  # the constructor function - hidden from user as it is embedded in ATE_estimation function
  create_ATE <- function(surv_ref, surv_trt, ref, trt_values, ATE, t) {
    #  surv_ref, the marginal survival probability for the reference
    #  surv_trt, the marginal survival probability for the treatment
    #  ref, the value of the reference treatment
    #  trt_values, the possible values of the treatment variable
    #  ATE, the difference between the marginal survival probability of the treatment and the reference
    #  return an object of class 'ATE'
    my_object <- structure(list(
      surv_ref = surv_ref, surv_trt = surv_trt, ref = ref, trt_values = trt_values, ATE = ATE, t = t
    ), class = "ATE")
    return(my_object)
  }
  
  ATE_object = create_ATE(surv_ref = surv_1_post, surv_trt = surv_2_post, 
                          trt_values = c(ref, trt_values[trt_values != ref]), 
                          ref = ref, ATE = ATE, t = t)
  
  return(ATE_object)
}
