vector[P] beta_sim;
vector[Q] g_sim;
vector[Q] eta_sim;
vector<lower=0>[sum_N_par] gamma_sim;
matrix[sum(N_dim),sum(N_dim)] L_D_sim;

for(p in 1:P){
  beta_sim[p] = normal_rng(prior_b_mean[p],prior_b_sd[p]);
}
for(q in 1:sum_N_par){
  gamma_sim[q] = fabs(normal_rng(0,prior_g_sd[q]));
}
for(q in 1:Q){
  eta_sim[q] = normal_rng(0,1);
}

L_D_sim = gen_d(gamma_sim,
                B,
                N_dim, 
                N_func, 
                func_def,
                N_var_func, 
                col_id,
                N_par,
                sum_N_par,
                cov_data);
g_sim = L_D_sim * eta_sim;

