functions {
#include /stan_files/gen_d.stan
}
data {
#include /stan_files/data.stan

  vector[N] y;
  real sigma_sd;
}
parameters {
#include /stan_files/params.stan

  real<lower=0> sigma;
}
model{
#include /stan_files/priors.stan
  
  sigma ~ normal(0,sigma_sd);
  y ~ normal(X*beta + Z*g,sigma);
} 
