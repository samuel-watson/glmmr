data {
#include /stan_files/mcml_data.stan

  int y[N];
  int type; // 1 = log
}
parameters {
  vector[Q] gamma;
}
model {
  vector[Q] zeroes = rep_vector(0,Q);
  gamma ~ multi_normal_cholesky(zeroes,L);
  if(type==1) y~poisson_log(Xb + Z*gamma);
}

