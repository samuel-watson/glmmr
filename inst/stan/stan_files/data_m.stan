// data to define the covariance function
int B_m; // number of blocks
int max_N_func_m;
int max_N_var_func_m;
int max_N_dim_m;
int max_N_var_m;
int N_dim_m[B_m+1]; //dimension of each block
int N_func_m[B_m+1]; //number of functions for each block
int func_def_m[B_m,max_N_func_m]; //type codes for each function
int N_var_func_m[B_m,max_N_func_m]; //number of variables in the dataset for each block
int col_id_m[max_N_func_m,max_N_var_func_m,B_m]; // column IDs of the data
int N_par_m[B_m,max_N_func_m]; // number of parameters for each function
int sum_N_par_m; // total number of covariance parameters
real cov_data_m[max_N_dim_m,max_N_var_m,B_m]; // data defining position

//data to define the model
int P_m; // columns of X
int Q_m; // columns of Z, size of RE terms
matrix[N,P_m] X_m;
matrix[N,Q_m] Z_m;

//prior parameters
vector[P_m] prior_b_mean_m;
vector[P_m] prior_b_sd_m;
vector[sum_N_par_m] prior_g_sd_m;

