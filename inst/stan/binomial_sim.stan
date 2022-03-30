functions {
#include /stan_files/gen_d.stan
}
data {
#include /stan_files/data.stan
  
  int type;
  int par_ind;
  real threshold;
}
transformed data {
  int y_sim[N];
  
#include /stan_files/tdata_sim.stan  
  
  for(i in 1:N){
    if(type==1) y_sim[i] = bernoulli_logit_rng(X[i,]*beta_sim + Z[i,]*g_sim);
    if(type==2) y_sim[i]=bernoulli_rng(exp(X[i,]*beta_sim + Z[i,]*g_sim));
    if(type==3) y_sim[i]=bernoulli_rng(X[i,]*beta_sim + Z[i,]*g_sim);
    if(type==4) y_sim[i]=bernoulli_rng(Phi_approx(X[i,]*beta_sim + Z[i,]*g_sim));
  } 
}
parameters {
#include /stan_files/params.stan
}
model{
  
#include /stan_files/priors.stan
  
  if(type==1) y_sim ~ bernoulli_logit(X*beta + Z*g);
  if(type==2) y_sim~bernoulli(exp(X*beta + Z*g));
  if(type==3) y_sim~bernoulli(X*beta + Z*g);
  if(type==4) y_sim~bernoulli(Phi_approx(X*beta + Z*g));
} 
generated quantities {
  int<lower=0,upper=1> betaIn = beta[par_ind] < beta_sim[par_ind];
  int<lower=0,upper=1> betaThresh = beta[par_ind] > threshold;
}

