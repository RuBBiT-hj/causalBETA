#' Construct an `ATE` object to store the results from g-computation
#' 
#' @param surv_ref a numeric matrix for the posterior estimates for the reference
#' @param surv_trt a numeric matrix for the posterior estimates for the treatment
#' @param ref the value of the reference treatment
#' @param trt_values a vector for the possible values of the treatment variable (binary)
#' @param estimand the estimand for this ATE
#' @param ATE a numeric matrix for the difference between the marginal survival probability 
#' of the treatment and the reference
#' @param t timepoints used in this ATE object, specifically for posterior survival probability
#' @param threshold optional; the threshold used for restrictive mean survival time. If not provided,
#' it will use the maximum time of the study as the default
#' survival time as a special case of restrictive mean survival time 
#' 
#' @return It returns an object of class 'ATE'

create_ATE <- function(surv_ref, surv_trt, ref, trt_values, ATE, estimand, t, threshold) {
  my_object <- structure(list(
    surv_ref = surv_ref, surv_trt = surv_trt, ref = ref, trt_values = trt_values, 
    estimand = estimand, ATE = ATE, t = t, threshold = threshold
  ), class = "ATE")
  return(my_object)
}