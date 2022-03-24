matrix gen_d(vector gamma,
             int B,
             int[] N_dim, 
             int max_N_dim, 
             int[] N_func, 
             int max_N_func,
             int[,] func_def,
             int[,] N_var_func, 
             int max_N_var,
             int max_N_var_func,
             int[,,] col_id,
             int[,] N_par,
             int sum_N_par,
             matrix[] cov_data){
               
  int sum_N_dim = sum(N_dim);             
  matrix[sum_N_dim,sum_N_dim] D = rep_matrix(0,sum_N_dim,sum_N_dim);
  matrix[sum_N_dim,sum_N_dim] L_D;
  int idx = 1;
  int g1_idx = 1;
  int g1_idx_ii = 1;
  
//loop over blocks
for(b in 1:B){
  //loop over the elements of the submatrix
  for(i in idx:(idx+N_dim[b]-2)){
    for(j in (i+1):(idx+N_dim[b]-1)){
      real val = 1;
      //loop over all the functions
      int gamma_idx = g1_idx;
      for(k in 1:N_func[b]){
        // generate the distance
        real dist = 0;
        for(p in 1:N_var_func[b,k]){
          dist += pow(cov_data[b,i+1-idx,col_id[b,k,p]] - 
            cov_data[b,j+1-idx,col_id[b,k,p]],2);
        }
        dist = sqrt(dist);
         // now to generate the right function
      if(func_def[b,k] == 1){
        // group member ship
        if(dist == 0){
          val = val*pow(gamma[gamma_idx],2);
        } else {
          val = 0;
        }
      } else if(func_def[b,k] == 2){
        // exponential 1
        val = val * gamma[gamma_idx]*exp(-1*dist*gamma[gamma_idx+1]);
      } else if(func_def[b,k] == 3){
        // exponential 2 power
        val = val * pow(gamma[gamma_idx],dist);
      }
     D[i,j] = val;
     D[j,i] = val;
     gamma_idx += N_par[b,k]; 
      }
      
    }
  }
  
  for(i in idx:(idx+N_dim[b]-1)){
      real val = 1;
      //loop over all the functions
      int gamma_idx_ii = g1_idx_ii;
      for(k in 1:N_func[b]){
        // generate the distance
        real dist = 0;
         // now to generate the right function
      if(func_def[b,k] == 1){
        // group member ship
        if(dist == 0){
          val = val*pow(gamma[gamma_idx_ii],2);
        } else {
          val = 0;
        }
      } else if(func_def[b,k] == 2){
        // exponential 1
        val = val * gamma[gamma_idx_ii]*exp(-1*dist*gamma[gamma_idx_ii+1]);
      } else if(func_def[b,k] == 3){
        // exponential 2 power
        val = val * pow(gamma[gamma_idx_ii],dist);
      }
     D[i,i] = val;
     gamma_idx_ii += N_par[b,k]; 
      }
    }
  g1_idx += sum(N_par[b,]);
  g1_idx_ii += sum(N_par[b,]);
  idx += N_dim[b];
}
return cholesky_decompose(D);
//g ~ multi_normal_cholesky(zeros,L_D);
} 
