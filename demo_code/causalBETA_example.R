## install and load causalBETA
#devtools::install_github("RuBBiT-hj/causalBETA") 
library(causalBETA)
## load other packages
library(survival)

## load data and re-code variables
data = survival::veteran
data$A = 1*(data$trt==2)

## rename variables
var_names = colnames(data)
colnames(data)[var_names=='status'] = 'delta'
colnames(data)[var_names=='time'] = 'y'

#------------------------------------------------------------------------------#
###       Run Unadjusted Analyses                                            ###
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
###       Run Model with independent hazards across increments               ###
#------------------------------------------------------------------------------#
set.seed(1) ## set seed so MCMC draws are reproducible
post_draws_ind = bayeshaz(d = data, ## data set
                          reg_formula = Surv(y, delta) ~ A, ## hazard regression formula
                          num_partition = 100, ## number of partitions, K
                          model = 'independent', ## prior on baseline hazard
                          sigma = 3, ## prior standard deviation for coefficients 
                          A = 'A', ## column name of treatment variable in d
                          warmup = 1000, ## number of warmup/burnin iterations
                          post_iter = 1000, ## number of post-warmup iterations to output
                          chains = 1)

#------------------------------------------------------------------------------#
###          Run Model with Smoothed hazards across increments               ###
#------------------------------------------------------------------------------#
set.seed(122132)
post_draws_ar1 = bayeshaz(d = data, ## data set
                          reg_formula = Surv(y, delta) ~ A,
                          model = 'AR1', ## choice of smoothing prior - independent or AR1 smoothing
                          A = 'A', ## which variable is treatment
                          warmup = 1000, post_iter = 1000)

png(filename = "hazard_plots_unadj.png", width = 800, height = 400)
par(mfrow=c(1,2))
## ----- plot results with independence prior ----- ##
plot(post_draws_ind, 
     ylim=c(0,.11), xlim=c(0, 900),
     type='p', level_CI = .95,
     main='Independent Prior Process', ylab = 'Baseline Hazard Rate', 
     xlab = 'Time (days)')

abline(h=mean(colMeans(post_draws_ind$haz_draws)), lty=2, col='black')

legend(x=0, y=.111, 
       legend = c('Posterior Point/Interval Estimate', 'Frequentist Estimate'), 
       col=c('black', 'red'), pch=c(20,20), bty='n')

## overlay frequentist point estimates 
freq_res = eha::pchreg(data=data,
                       cuts = post_draws_ind$partition,
                       formula = Surv(y, delta) ~  A)
points(post_draws_ind$midpoint, freq_res$hazards, col='red',pch=20, cex=.5)

## ----- plot results with AR1 prior ----- ##
plot(post_draws_ar1, 
     ylim=c(0,.11), xlim=c(0, 900), 
     type='p', level_CI = .95,
     main='AR1 Prior Process',
     ylab = 'Baseline Hazard Rate', 
     xlab = 'Time (days)')
points(post_draws_ar1$midpoint, freq_res$hazards, col='red',pch=20, cex=.5)
dev.off()

#------------------------------------------------------------------------------#
###         Adjusted Hazard Models Curves                                    ###
#------------------------------------------------------------------------------#
set.seed(123) ## set seed so MCMC draws are reproducible

formula1 = Surv(y, delta) ~ A + age + karno + celltype

post_draws_ar1_adj = bayeshaz(d = data, 
                               reg_formula = formula1,
                               model = 'AR1',
                               A = 'A', 
                               warmup = 1000, post_iter = 1000) 

summary(post_draws_ar1_adj$beta_draws, quantiles = c(.025, .975))
summary(eha::pchreg(data=data,cuts = post_draws_ar1_adj$partition,formula = formula1))

set.seed(32123)
gcomp_res = bayesgcomp(post_draws_ar1_adj, ## posterior draws of hazard model
                       ref = 0, ## treatment reference group
                       B = 1000, ## monte carlo iterations in g-comp
                       estimand = "prob") ## Posterior Survival Prob 

summary( gcomp_res$ATE, quantiles = c(.025, .975) )

#------------------------------------------------------------------------------#
###         Adjusted Survival Curves                                         ###
#------------------------------------------------------------------------------#

png(filename = "survival_plots_adj.png", width = 700, height = 300)
par(mfrow=c(1,3))
## survival curve under A=0
plot(gcomp_res, mode = 0, type ='p', 
     ylim=c(0,1), xlim=c(0, 1000), 
     xlab='Time (days)', ylab='Survival Probability',
     main='Marginal Survival Curve under Standard Chemo')
## overlay kaplan-meier
d0 = data[data$A==0, ]
lines(survfit(data=d0, Surv(y, delta)~ 1))

## survival curve under A=1
plot(gcomp_res, mode = 1, type ='p', 
     ylim=c(0,1), xlim=c(0, 1000), 
     xlab='Time (days)', ylab='Survival Probability',
     main='Marginal Survival Curve under Novel Chemo')
## overlay kaplan-meier
d1 = data[data$A==1, ]
lines(survfit(data=d1, Surv(y, delta)~ 1))

## difference in survival rate
plot(gcomp_res, mode = 'ATE', ylim=c(-1,1), 
     main = 'Difference in Marginal Survival Probability',
     xlab = 'Time (days)', 
     ylab = 'Difference')
abline(h=0, lty=2, col='red')
dev.off()


#------------------------------------------------------------------------------#
###         Convergence Checks                                               ###
#------------------------------------------------------------------------------#

## run three chains:
n_chains = 3

set.seed(1)
post_draws = bayeshaz(d = data, 
                      reg_formula = formula1,
                      model = 'AR1',
                      A = 'A', 
                      warmup = 2000, 
                      post_iter = 1000, ## output 1000 draws
                      chains = n_chains)

ATE_chains = bayesgcomp(post_draws,
                        ref = 0,
                        B = 1000,
                        estimand = "prob",
                        t = c(365, 2*365))

## summarize posterior draws of the three chains;
## compute posterior mean and 95% credible interval endpoints
## for 1-yr and 2-yr survival rate differences
summary(ATE_chains$ATE, quantiles = c(.025, .975))

## plot traceplots and density

png(filename = "traceplots.png", width = 700, height = 500)
plot(ATE_chains$ATE, ask = F)
dev.off()

## at each iteration, plot median, .025, and .975 percentiled 
## of previous iterations - these quantiles should stabilize over iterations.
coda::cumuplot(ATE_chains$ATE)

## if converged, gelman-rubin diagnostic Upper C.I. should be near 1
coda::gelman.diag(ATE_chains$ATE)

#------------------------------------------------------------------------------#
###         Monte Carlo Iterations Checks                                    ###
#------------------------------------------------------------------------------#
set.seed(32123)

B_vec = c(1, 100, 500, 1000)

ATE_list = list()

for(B in B_vec){
        gcomp_res = bayesgcomp(post_draws_ar1_adj, ## posterior draws of hazard model
                               ref = 0, ## treatment reference group
                               B = B, ## monte carlo iterations in g-comp
                               t = 365)
        ATE_list[[length(ATE_list) + 1 ]] = gcomp_res$ATE
}

png(filename = 'density.png', width = 700, height = 500)
par(mfrow=c(1,1))
plot(density(ATE_list[[1]][[1]]), ylim=c(0, 20), col='red', 
     main = latex2exp::TeX(paste0('Posterior Density of $\\Psi(t)$ for $t=365$') ) )
lines(density(ATE_list[[2]][[1]]), col='black')
lines(density(ATE_list[[3]][[1]]), col='blue')
lines(density(ATE_list[[4]][[1]]), col='green')
legend('topleft', 
       legend = c('B=1', 'B=100', 'B=500', 'B=1000'), 
       col = c('red','black','blue','green'), 
       lty=c(1,1,1,1), bty='n')
dev.off()