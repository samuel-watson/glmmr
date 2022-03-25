functions {
#include /stan_files/gen_d.stan
}
data {
#include /stan_files/data.stan
  
#include /stan_files/data_m.stan
  
  int type;
  int par_ind;
  real threshold;
}
transformed data {
  vector[N] y_sim;
  
#include /stan_files/tdata_sim_m.stan  
  
  for(i in 1:N){
    if(type==1) y_sim[i] = poisson_log_rng(X*beta_sim + Z*g_sim);
  } 
}
parameters {
#include /stan_files/params.stan
}
model{
  
#include /stan_files/priors.stan
  
  if(type==1) y ~ poisson_log(X*beta + Z*g);
} 
generated quantities {
  int<lower=0,upper=1> betaIn = beta[par_ind] < beta_sim[par_ind];
  int<lower=0,upper=1> betaThresh = beta[par_ind] > threshold;
}