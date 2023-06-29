data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> n_pieces;
  int delta[N];
  real offset[N];
  int interval_num[N];
  matrix[N, P] xmat;
}

parameters {
  vector[n_pieces] haz_eps;
  real eta;
  real<lower=0> sigma_haz;
  real<lower=0, upper=1> rho_eps;
  vector[P] beta;
}

transformed parameters{
  
  vector[n_pieces] log_haz;
  real<lower=-1, upper=1> rho;
  vector[N] lin_comb;
  
  rho = 2*(rho_eps - .5);
  
  // specify smoothing prior over increments 
  
  log_haz[1] = eta + sigma_haz*haz_eps[1];
  for(i in 2:n_pieces){
    log_haz[i] = eta*(1-rho) + rho*log_haz[i-1] + sigma_haz*haz_eps[i];
  }

  // independent priors 
  // log_haz[1] = eta + sigma*haz_eps[1];
  // for(i in 2:n_pieces){
  //   log_haz[i] = eta + sigma*haz_eps[i];
  // }

  
  // specify unsmoothed prior (i.e. independent hazards across increments )
  
  // log_haz[1] = -1*((sigma^2)/2) + sigma*haz_eps[1];
  // for(i in 2:n_pieces){
  //   log_haz[i] =  -1*((sigma^2)/2) + sigma*haz_eps[i];
  // }

  
  lin_comb = xmat * beta;
}

model {

  eta ~ normal(0, 1);
  rho_eps ~ beta(2,2);
  haz_eps ~ normal(0, 1);
  beta ~ normal(0,3);
  
  // likelihood contributions
  for(i in 1:N){
    delta[i] ~ poisson( exp( log_haz[ interval_num[i] ] + lin_comb[i] + log(offset[i]) ) );  
  }
  
  
}

