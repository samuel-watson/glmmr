vector[Q] g;
{
  matrix[sum(N_dim),sum(N_dim)] L_D = gen_d(gamma,
                                        B,
                                        N_dim, 
                                        N_func, 
                                        func_def,
                                        N_var_func, 
                                        col_id,
                                        N_par,
                                        sum_N_par,
                                        cov_data);
  g = L_D * eta;
}
for(p in 1:P){
  beta[p] ~ normal(prior_b_mean[p],prior_b_sd[p]);
}
for(q in 1:sum_N_par){
  gamma[q] ~ normal(0,prior_g_sd[q]);
}
eta ~ std_normal();

