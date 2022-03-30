functions {
#include /stan_files/gen_d.stan
}
data {
#include /stan_files/data.stan
  
  real sigma_sd;
  int par_ind;
  real threshold;
}
transformed data {
  real<lower=0> sigma_sim;
  vector[N] y_sim;

#include /stan_files/tdata_sim.stan  

  sigma_sim = fabs(normal_rng(0,sigma_sd));
  for(i in 1:N){
    y_sim[i] = normal_rng(X[i,]*beta_sim + Z[i,]*g_sim,sigma_sim);
  } 
}
parameters {
#include /stan_files/params.stan
  
  real<lower=0> sigma;
}
model{

#include /stan_files/priors.stan
  
  sigma ~ normal(0,sigma_sd);
  y_sim ~ normal(X*beta + Z*g,sigma);
} 
generated quantities {
  int<lower=0,upper=1> betaIn = beta[par_ind] < beta_sim[par_ind];
  int<lower=0,upper=1> betaThresh = beta[par_ind] > threshold;
}

