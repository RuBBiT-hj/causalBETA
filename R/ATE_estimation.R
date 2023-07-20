#' Average Treatment Effect Estimation
#' 
#' Estimate the average treatment effect (ATE) based on the predicted survival times
#' @param d the data
#' @param beta_draws the beta coefficients from the posterior draws, with dimension number of posterior draws * the number of variables 
#' @param haz_draws the baseline hazard rates from the posterior draws, with dimension number of posterior draws * the number of intervals
#' @param partition the time values for each partition
#' @param covariates the names of the covariates used in the model
#' @param trt the treatment variable name
#' @param trt_values the possible values of the treatment, and the length should be two: reference and treatment
#' @param n the number of prediction for each posterior draw; the default is 1000
#' @examples
#' # example demo
## usethis namespace: start
#' @import survival
#' @importFrom mets rpch
#' @importFrom LaplacesDemon rdirichlet
## usethis namespace: end
#' @export


ATE_estimation = function(d, beta_draws, haz_draws, partition, covariates, trt, trt_values, 
                          n = 1000){
  
  # calculate the survival probability for each subject at each time
  # let t be the middle points of each partition (the number is the same as the intervals)
  t <- diff(partition)/2 + partition[-length(partition)]
  
  # the number of subjects
  n_subject <- nrow(d)
  
  # make two datasets with two treatment values
  d_1 <- d[, covariates]
  d_1[, trt] <- trt_values[1]
  
  d_2 <- d[, covariates]
  d_2[, trt] <- trt_values[2]
  
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
  
  ATE <- matrix(nrow = nrow(beta_draws),
                ncol = length(t))
  
  # alpha vector for rdirichlet
  alpha <- rep(1, n_subject)
  
  for (i in 1:nrow(beta_draws)){
    
    # draw Dirichlet for weighted sum
    weights_dir <- matrix(LaplacesDemon::rdirichlet(1, alpha = alpha), ncol = 1)
    
    diff_matrix <- sapply(surv_prob_diff, function(x) x[i,]) # n_parition x n_subject
    ATE[i, ] <- diff_matrix %*% weights_dir
  }
  
  return(list(t, ATE))
}
