functions {
#include /stan_files/gen_d.stan
}
data {
#include /stan_files/data.stan
  
  int y[N];
  int type;
}
parameters {
#include /stan_files/params.stan
  
  real<lower=0> sigma;
}
model{
#include /stan_files/priors.stan
  
  if(type==1) y ~ poisson_log(X*beta + Z*g);
} 

