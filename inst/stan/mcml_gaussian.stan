data {
  int N; // sample size
  int P; // columns of X
  int Q; // columns of Z, size of RE terms
  vector[N] Xb;
  matrix[Q,Q] L;
  matrix[N,Q] Z;
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
