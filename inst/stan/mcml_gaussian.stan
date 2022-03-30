data {
#include /stan_files/mcml_data.stan

  vector[N] y;
  real sigma;
  int type;
}
parameters {
  vector[Q] gamma;
}
model {
  vector[Q] zeroes = rep_vector(0,Q);
  gamma ~ multi_normal_cholesky(zeroes,L);
  if(type==1)y ~ normal(Xb + Z*gamma,sigma);
}

