vector[P_m] beta_sim;
vector[Q_m] g_sim;
vector[Q_m] eta_sim;
vector<lower=0>[sum_N_par_m] gamma_sim;
matrix[sum(N_dim_m),sum(N_dim_m)] L_D_sim;

for(p in 1:P_m){
  beta_sim[p] = normal_rng(prior_b_mean_m[p],prior_b_sd_m[p]);
}
for(q in 1:sum_N_par_m){
  gamma_sim[q] = fabs(normal_rng(0,prior_g_sd_m[q]));
}
for(q in 1:Q_m){
  eta_sim[q] = normal_rng(0,1);
}

L_D_sim = gen_d(gamma_sim,
                B_m,
                N_dim_m, 
                N_func_m, 
                func_def_m,
                N_var_func_m, 
                col_id_m,
                N_par_m,
                sum_N_par_m,
                cov_data_m);
g_sim = L_D_sim * eta_sim;

