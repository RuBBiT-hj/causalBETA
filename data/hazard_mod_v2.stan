data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> n_pieces;
  array[N] int delta;
  array[N] real off_set;
  array[N] int interval_num;
  matrix[N, P] xmat;
  real<lower=0> sigma_beta;
}

parameters {
  vector[n_pieces] haz_eps;
  real eta;
  real<lower=0> sigma_haz [n_pieces];
  real<lower=0, upper=1> rho_eps;
  vector[P] beta;
}

transformed parameters{
  
  vector[n_pieces] log_haz;
  real<lower=-1, upper=1> rho;
  vector[N] lin_comb;
  
  rho = 2*(rho_eps - .5);
  
  // specify smoothing prior over increments 
  log_haz[1] = eta + sigma_haz[1]*haz_eps[1];
  for(i in 2:n_pieces){
    log_haz[i] = eta*(1-rho) + rho*log_haz[i-1] + sigma_haz[i]*haz_eps[i];
  }

  lin_comb = xmat * beta;
}

model {

  eta ~ normal(0, 1);
  rho_eps ~ beta(2,2);
  haz_eps ~ normal(0, 1);
  // user-defined odds ratio prior, B should be less than 3 but greater than 0
  beta ~ normal(0, sigma_beta);
  sigma_haz ~ gamma(1,1);
  
  // likelihood contributions
  for(i in 1:N){
    delta[i] ~ poisson( exp( log_haz[ interval_num[i] ] + lin_comb[i] + log(off_set[i]) ) );  
  }
  
  
}

