#' Construct an `ATE` object to store the results from g-computation
#' 
#' @param surv_ref a numeric matrix for the marginal survival probability for the reference
#' @param surv_trt a numeric matrix for the marginal survival probability for the treatment
#' @param ref the value of the reference treatment
#' @param trt_values a vector for the possible values of the treatment variable (binary)
#' @param ATE a numeric matrix for the difference between the marginal survival probability 
#' of the treatment and the reference
#' 
#' @return It returns an object of class 'ATE'

create_ATE <- function(surv_ref, surv_trt, ref, trt_values, ATE, t) {
  my_object <- structure(list(
    surv_ref = surv_ref, surv_trt = surv_trt, ref = ref, trt_values = trt_values, ATE = ATE, t = t
  ), class = "ATE")
  return(my_object)
}