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
  
  if(type==1) y ~ bernoulli_logit(X*beta + Z*g);
  if(type==2) y~bernoulli(exp(X*beta + Z*g));
  if(type==3) y~bernoulli(X*beta + Z*g);
  if(type==4) y~bernoulli(Phi_approx(X*beta + Z*g));
} 

