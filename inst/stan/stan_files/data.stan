// data to define the covariance function
int B; // number of blocks
int N_dim[B]; //dimension of each block
int max_N_dim; //maximum of N_dim
int N_func[B]; //number of functions for each block
int max_N_func; //max of Fu
int func_def[B,max_N_func]; //type codes for each function
int N_var_func[B,max_N_func]; //number of variables in the dataset for each block
int max_N_var; //maximum number of columns of data for each block
int max_N_var_func; //maximum number of variables for a function
int col_id[B,max_N_func,max_N_var_func]; // column IDs of the data
int N_par[B,max_N_func]; // number of parameters for each function
int sum_N_par; // total number of covariance parameters
matrix[max_N_dim,max_N_var] cov_data[B]; // data defining position

//data to define the model
int N; // sample size
int P; // columns of X
int Q; // columns of Z, size of RE terms
matrix[N,P] X;
matrix[N,Q] Z;

//prior parameters
vector[P] prior_b_mean;
vector[P] prior_b_sd;
vector[sum_N_par] prior_g_sd;