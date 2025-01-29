#' Modify the Reference Value
#' 
#' The reference value of the `ATE` class object is changed to be the argument `ref`. If the new value provided is the
#' other treatment, then this function will exchange the marginal survival probability between reference and treatment.
#' @param x an object of the class `ATE` created by the `bayesgcomp()` function
#' @param ref the new reference value, so it should be one of the treatment values
#' @param ... additional arguments for future methods
#' @details
#' When the reference provided is the other treatment than the current `ref` in the `ATE` object,
#' this function swaps the marginal survival probabilities, but it doesn't recalculate them for computational
#' efficiency. For calculation with a different set of weights, a new `ATE` object needs to be generated from
#' the `ATE_estimation()` function.
#' @returns An object of `ATE` with the new reference
#' @examples
#' \dontrun{
#' # example demo
#' ## Continued from ?bayesgcomp
#' ## original reference is 0
#' gcomp_res_relevel <- relevel(gcomp_res, ref = 1)
#' plot(gcomp_res_relevel)
#' plot(gcomp_res) ## comparison
#' }
## usethis namespace: start
#' @importFrom stats relevel
## usethis namespace: end
#' @export

relevel.ATE = function(x, ref, ...){
  # extract
  trt_values <- x$trt_values
  ref_current <- x$ref
  
  # if the ref is valid
  if ( !(ref %in% trt_values) ) stop("Reference provided is not a valid treatment value")
  # if the ref is the other treatment
  if ( ref != ref_current){
    temp <- x$surv_ref
    x$surv_ref <- x$surv_trt
    x$surv_trt <- temp
    x$ATE <- -x$ATE
  }
  
  x$ref <- ref
  return(x)
}
