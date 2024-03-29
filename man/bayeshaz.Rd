% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayeshaz.R
\name{bayeshaz}
\alias{bayeshaz}
\title{Bayesian Piece-Wise Exponential Model}
\usage{
bayeshaz(
  d,
  reg_formula,
  A,
  model = "AR1",
  sigma = 3,
  num_partitions = 100,
  warmup = 1000,
  post_iter = 1000
)
}
\arguments{
\item{d}{data, a data frame in survival format. Categorical variables should be
transformed into dummy variables}

\item{reg_formula}{a formula object that specifies the formula for the poisson regression.
This also decides the formula will be used the function to check positivity overlap.}

\item{A}{a character variable that specifies the name of the treatment}

\item{model}{a character variable that tells the stan model used to implement the Bayesian piece-wise exponential model,
default is "AR1" and the other option is "independent"}

\item{sigma}{a numeric variable as the user-defined standard deviation for beta coefficients prior, the default is 3}

\item{num_partitions}{a numeric variable as the number of partitions of the study time, the default is 100}

\item{warmup}{a numeric variable as the number of warmup in MCMC, the default is 1000}

\item{post_iter}{a numeric variable as the number of iterations to draw from the posterior, the default is 1000}
}
\value{
It returns an object of class \code{bayeshaz} that contains the information about the data, model, etc.
This serves as the basis for the extended functions in this package.

An object of class \code{bayeshaz} is a list containing at least the following components:
\itemize{
\item \code{data}, a data frame for the data
\item \code{formula}, the regression formula
\item \code{treatment}, the name of the treatment variable
\item \code{covariates}, the name of the covariate(s)
\item \code{time}, the name of the time variable
\item \code{outcome}, the name of the outcome variable
\item \code{model}, the type of the model used
\item \code{sigma}, the sigma specified
\item \code{partition}, the partition vector
\item \code{midpoint}, the midpoints of intervals
\item \code{haz_draws}, the baseline hazard rate from each posterior draws (as an \code{mcmc} object)
\item \code{beta_draws}, the beta coefficients estimated from each posterior draws (as an \code{mcmc} object)
}
}
\description{
Perform a Bayesian piece-wise exponential model on the given survival data,
and this function implements it by an equivalent poisson regression in MCMC.
}
\details{
A typical model has the form \code{Surv(time, outcome) ~ covariates}. The function will capture
the outcome and covariates based on this formula object.

The Bayesian piece-wise exponential model uses normal prior for baseline hazard rate,
beta coefficients, and the error term for baseline hazard rate to form the first-order Gaussian process.
It uses a beta prior to generate the correlation from -1 to 1 for hazard rate.
For more details, users can check the full model listed in the reference.

Under the \code{AR1} model, user can specify the standard deviation for the normal prior of the beta coefficients.
The default is 3, and only values between 0 and 3 are accepted since 3 is already a relatively weak prior.
}
\examples{
# example demo
data = survival :: veteran
data$A = 1*(data$trt==2)
## rename variables
var_names = colnames(data)
colnames(data)[var_names=='status'] = 'delta'
colnames(data)[var_names=='time'] = 'y'
## append one-hot encoded celltypes
data = cbind(data, model.matrix(data=data, ~ -1 + celltype))
formula1 = Surv(y, delta) ~ A + age + karno + 
  celltypesquamous + celltypesmallcell + celltypeadeno
post_draws_ar1_adj = bayeshaz(
  d = data,
  reg_formula = formula1 ,
  model = 'AR1',
  A = 'A',
  warmup = 1000, post_iter = 1000)
}
\references{
Ji, Han, and Oganisian, Arman. 2023.
"causalBETA: An R Package for Bayesian Semiparametric Causal Inference with
Event-Time Outcomes. \emph{arXiv:2310.12358 [Stat]}, October.
\url{http://arxiv.org/abs/2310.12358}.
}
