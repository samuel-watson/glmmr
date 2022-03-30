functions {
#include /stan_files/gen_d.stan
}
data {
#include /stan_files/data.stan
  
#include /stan_files/data_m.stan
  
  real sigma_sd;
  real sigma_sd_m;
  int par_ind;
  real threshold;
}
transformed data {
  real<lower=0> sigma_sim;
  vector[N] y_sim;
  
#include /stan_files/tdata_sim_m.stan  
  
  sigma_sim = fabs(normal_rng(0,sigma_sd_m));
  for(i in 1:N){
    y_sim[i] = normal_rng(X_m[i,]*beta_sim + Z_m[i,]*g_sim,sigma_sim);
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

