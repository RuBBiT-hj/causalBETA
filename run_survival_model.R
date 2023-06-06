library(survival)
library(cmdstanr)
library(eha)

setwd("~/Desktop/BayesianSurvivalReg/")

#### ----------- prepare data ---------- ####
d = survival::cancer
d = d[complete.cases(d),] ## drop missing data
d$id = 1:nrow(d) ## patient id
d$status[d$status==2] = 0
## recode binary variables to 0/1
d$sex[d$sex==2] = 0
d$age = ( d$age - mean(d$age) )/ sd(d$age) ## standardize age


#### ----------- run Bayesian piecewise constant hazard model ------------- ####

post_draws = bayeshaz(d = d, 
                      reg_formula = Surv(time, status) ~ age + ph.ecog + sex,
                      A = 'sex', num_intervals = 100, warmup = 1000, post_iter = 1000 )

#### ----------- Plot baseline hazard model ------------- ####

### plot posterior mean of baseline hazard estimate

## posterior mean/ 95% interval
bslhaz_mean = colMeans(post_draws$haz_draws)
bslhaz_lwr = apply(post_draws$haz_draws, 2, quantile, probs=.025)
bslhaz_upr = apply(post_draws$haz_draws, 2, quantile, probs=.975)

plot( post_draws$xv, bslhaz_mean, pch=20, ylim=c(0, .05) ) 
segments(x0 = post_draws$xv, y0 = bslhaz_lwr,
         x1 = post_draws$xv, y1 = bslhaz_upr)

## estimate frequentist piecewise exponential hazard function
freq_res = eha::pchreg(data=d,cuts = seq( 0 , max(d$time) + .01, length.out = 100 ),
                       formula = Surv(time, status) ~  age + ph.ecog + sex)
## overlay frequentist point estimates 
points(post_draws$xv, freq_res$hazards, col='red')


## compare coefficients with frequentist estimates
colMeans(post_draws$beta_draws)
apply(post_draws$beta_draws, 2, quantile, probs=c(.025, .975))
freq_res$coefficients





