matrix gen_d(vector gamma,
             int B,
             int[] N_dim, 
             int[] N_func, 
             int[,] func_def,
             int[,] N_var_func, 
             int[,,] col_id,
             int[,] N_par,
             int sum_N_par,
             real[,,] cov_data){
               
  int sum_N_dim = sum(N_dim);             
  matrix[sum_N_dim,sum_N_dim] D = rep_matrix(0,sum_N_dim,sum_N_dim);
  matrix[sum_N_dim,sum_N_dim] L_D;
  int idx = 1;
  int g1_idx = 1;
  int g1_idx_ii = 1;
  
//loop over blocks
for(b in 1:B){
  //loop over the elements of the submatrix
  int gr1 = 0;
  for(p in 1:N_func[b]) gr1 += func_def[b,p] == 1;
  if(gr1 == 0){
    for(i in idx:(idx+N_dim[b]-2)){
    for(j in (i+1):(idx+N_dim[b]-1)){
      real val = 1;
      //loop over all the functions
      int gamma_idx = g1_idx;
      for(k in 1:N_func[b]){
        // generate the distance
        real dist = 0;
        for(p in 1:N_var_func[b,k]){
          dist += pow(cov_data[i+1-idx,col_id[k,p,b],b] - 
            cov_data[j+1-idx,col_id[k,p,b],b],2);
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
        val = val * exp(-1*dist*gamma[gamma_idx]);
      } else if(func_def[b,k] == 3){
        // ar1
        val = val * pow(gamma[gamma_idx],dist);
      } else if(func_def[b,k] == 4){
        // squared exponential
        val = val * gamma[gamma_idx]*exp(-1*pow(dist,2)*
          gamma[gamma_idx+1]^2);
      } else if(func_def[b,k] == 5){
        // matern
        real xr = pow(2*gamma[gamma_idx],0.5)*dist/gamma[gamma_idx+1];
        val = val * (pow(2, -1*(gamma[gamma_idx]-1))/tgamma(gamma[gamma_idx]))*
          pow(xr, gamma[gamma_idx])*modified_bessel_second_kind(1,xr);
      } else if(func_def[b,k] == 6){
        // bessel
        real xr = dist/gamma[gamma_idx];
        val = val * xr * modified_bessel_second_kind(1,xr);
      } 
     
     gamma_idx += N_par[b,k]; 
      }
      D[i,j] = val;
     D[j,i] = val;
    }
  }
  }
  
  
  for(i in idx:(idx+N_dim[b]-1)){
      real val = 1;
      //loop over all the functions
      int gamma_idx_ii = g1_idx_ii;
      for(k in 1:N_func[b]){
         // now to generate the right function
      if(func_def[b,k] == 1){
        // group member ship
        val = val*pow(gamma[gamma_idx_ii],2);
      } 
      // else if(func_def[b,k] == 4){
      //   // exponential 1func_def[b,k] == 2 || 
      //   val = val * gamma[gamma_idx_ii];
      // } 
     gamma_idx_ii += N_par[b,k]; 
      }
      D[i,i] = val;
    }
  g1_idx += sum(N_par[b,]);
  g1_idx_ii += sum(N_par[b,]);
  idx += N_dim[b];
}
return cholesky_decompose(D);
} 

