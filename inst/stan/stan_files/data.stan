// data to define the covariance function
int B; // number of blocks
int max_N_func;
int max_N_var_func;
int max_N_dim;
int max_N_var;
int N_dim[B+1]; //dimension of each block - need to pad these to prevent problems for B=1
int N_func[B+1]; //number of functions for each block
int func_def[B,max_N_func]; //type codes for each function
int N_var_func[B,max_N_func]; //number of variables in the dataset for each block
int col_id[max_N_func,max_N_var_func,B]; // column IDs of the data
int N_par[B,max_N_func]; // number of parameters for each function
int sum_N_par; // total number of covariance parameters
real cov_data[max_N_dim,max_N_var,B]; // data defining position

//data to define the model
int N; // sample size
int P; // columns of X
int Q; // columns of Z, size of RE terms
matrix[N,P] X;
matrix[N,Q] Z;

//prior parameters
vector[P+1] prior_b_mean;
vector[P+1] prior_b_sd;
vector[sum_N_par+1] prior_g_sd; //padding to read data from R of length 1

