#' Positivity Assumption Checking
#' 
#' Plot the histogram of the estimated propensity score by observed treatment. 
#' This plot can be used to check the positivity assumption for casual inference.
#' If the assumption is met, there should be good overlap in the histograms. Histograms that do not overlap may indicate near violations of positivity.
#' 
#' @param bayeshaz_object an object of the class `bayeshaz` created by the `bayeshaz()` function
#' @param formula an optional formula variable used in fitting the propensity score model
#' @param breaks the parameter set for breaks, default is `"Scott"`
#' 
#' @details
#' If the formula is not given, this function will constructs the formula 
#' for propensity score model by removing all treatment related terms from
#' the survival model formula, while the outcome should be the treatment variable.
#' If users provide a formula, it should not contain treatment variable nor the variables
#' that are not used in the survival model.
#' 
#' @examples
#' # example demo
## usethis namespace: start
## usethis namespace: end
#' @export

plot_positivity = function(bayeshaz_object, formula = NULL, breaks="Scott"){
  
  # extract
  d = bayeshaz_object$data
  
  reg_formula = bayeshaz_object$formula
  covariates = bayeshaz_object$covariates
  treatment = bayeshaz_object$treatment
  
  # if the formula is not given, we construct from the model formula
  if (is.null(formula)) {
    # remove the terms and interaction terms related to treatment
    variables <- attr(terms(reg_formula), "term.labels")
    variables <- variables[!grepl(treatment, variables, fixed = TRUE)]
    # construct the formula
    formula = formula(
    paste0(treatment, "~", paste0(variables, collapse = "+"))
    )
  }
  else{
    # check if the formula is valid in terms of the covariate scope
    variables <- as.character(attr(terms(formula), "variables"))
    # remove list and outcome
    variables <- variables[3:length(variables)]
      
    if (treatment %in% variables) warning("Treatment varibale is in the formula given")
    if ( !all(variables %in% covariates) ) warning("Formula has variables that are not
                                                   listed in the survival model")
  }
  
  # output the formula
  formula_text <- as.character(formula)
  formula_text <- paste0(formula_text[2], formula_text[1],  formula_text[3], 
                         collapse = " ")
  message(formula_text)
  
  # Fit the logistic regression
  pos_check <- glm(formula, family = binomial(link = "logit"),
                   data = d)
  d_trt <- d[, treatment]
  trt_1 <- unique(d_trt)[1]
  trt_2 <- unique(d_trt)[2]
  
  
  score_1 <- predict(pos_check, newdata = d[d_trt == trt_1, ], type = "response")
  score_2 <- predict(pos_check, newdata = d[d_trt == trt_2, ], type = "response")
  
  # Get parameters for histogram
  h_1 <- hist(score_1, breaks = breaks, probability = T)
  h_2 <- hist(score_2, breaks = breaks, probability = T)
  
  x_min = min(c(score_1, score_2))*0.9
  x_max = min(1, max(c(score_1, score_2))*1.1) # not exceed 1
  y_max = 1.02 * max(c(h_1$density, h_2$density))
  
  message(paste0("The above is ", trt_1, " and the below is ", trt_2))
  
  # Plot
  par(mfrow=c(2,1), mar=c(0,5,3,3))
  hist(score_1 , main="Distribution of Estimated P-Scores" , ylab="Density", xlab="", 
       ylim=c(0, y_max), xlim = c(x_min, x_max), 
       xaxt="n", las=1 , col="slateblue1", breaks=breaks, probability = T)
  par(mar=c(5,5,0,3))
  hist(score_2 , main="" , ylab="Density", xlab="Estimated Propensity Score", 
       ylim=c(y_max, 0),  xlim = c(x_min, x_max),
       las=1 , col="tomato3", breaks = breaks, probability = T)
  # reset
  par(mfrow=c(1,1),
      mar = c(5.1, 4.1, 4.1, 2.1),
      mgp = c(3,1,0))
  
}

