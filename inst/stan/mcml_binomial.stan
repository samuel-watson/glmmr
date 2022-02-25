data {
  int N; // sample size
  int P; // columns of X
  int Q; // columns of Z, size of RE terms
  vector[N] Xb;
  matrix[Q,Q] L;
  matrix[N,Q] Z;
  int y[N];
  int type; // 1 = logit, 2= log, 3=identity, 4= probit
}
parameters {
  vector[Q] gamma;
}
model {
  vector[Q] zeroes = rep_vector(0,Q);
  gamma ~ multi_normal_cholesky(zeroes,L);
  if(type==1)  y ~ bernoulli_logit(Xb + Z*gamma);
  if(type==2) y~bernoulli(exp(Xb + Z*gamma));
  if(type==3) y~bernoulli(Xb + Z*gamma);
  if(type==4) y~bernoulli(Phi_approx(Xb + Z*gamma));
}
