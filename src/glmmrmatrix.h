#ifndef GLMMRMATRIX_H
#define GLMMRMATRIX_H

#include <RcppArmadillo.h>
#include "glmmrmath.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Combines a field of matrices into a block diagonal matrix
//' 
//' Combines a field of matrices into a block diagonal matrix. Used on
//' the output of `genD`
//' @param matfield A field of matrices
//' @return A block diagonal matrix
// [[Rcpp::export]]
arma::mat blockMat(arma::field<arma::mat> matfield){
  arma::uword nmat = matfield.n_rows;
  if(nmat==1){
    return matfield(0);
  } else {
    arma::mat mat1 = matfield(0);
    arma::mat mat2;
    if(nmat==2){
      mat2 = matfield(1);
    } else {
      mat2 = blockMat(matfield.rows(1,nmat-1));
    }
    arma::uword n1 = mat1.n_rows;
    arma::uword n2 = mat2.n_rows;
    arma::mat dmat(n1+n2,n1+n2);
    dmat.fill(0);
    dmat.submat(0,0,n1-1,n1-1) = mat1;
    dmat.submat(n1,n1,n1+n2-1,n1+n2-1) = mat2;
    return dmat;
  }
}

//' Generates a block of the random effects covariance matrix
//' 
//' Generates a block of the random effects covariance matrix
//' @details 
//' Using the sparse representation of the random effects covariance matrix, constructs
//' one of the blocks. The function definitions are: 1 indicator, 2 exponential,
//' 3 AR-1, 4 squared exponential, 5 matern, 6 Bessel.
//' @param N_dim Integer specifying the dimension of the matrix
//' @param N_func Integer specifying the number of functions in the covariance function 
//' for this block.
//' @param func_def Vector of integers of same length as `func_def` specifying the function definition for each function. 
//' @param N_var_func Vector of integers of same length as `func_def` specying the number 
//' of variables in the argument to the function
//' @param col_id Matrix of integers of dimension length(func_def) x max(N_var_func) that indicates
//' the respective column indexes of `cov_data` 
//' @param N_par Vector of integers of same length as `func_def` specifying the number
//' of parameters in the function
//' @param cov_data Matrix holding the data for the covariance matrix
//' @param gamma Vector of covariance parameters specified in order they appear in the functions 
//' specified by `func_def`
//' @return A symmetric positive definite matrix
// [[Rcpp::export]]
arma::mat genBlockD(size_t N_dim,
                    size_t N_func,
                    const arma::uvec &func_def,
                    const arma::uvec &N_var_func,
                    const arma::umat &col_id,
                    const arma::uvec &N_par,
                    const arma::mat &cov_data,
                    const arma::vec &gamma){
  arma::mat D(N_dim,N_dim,fill::zeros);
  if(!all(func_def == 1)){
    for(arma::uword i=0;i<(N_dim-1);i++){
      for(arma::uword j=i+1;j<N_dim;j++){
        double val = 1;
        size_t gamma_idx = 0;
        for(arma::uword k=0;k<N_func;k++){
          double dist = 0;
          for(arma::uword p=0; p<N_var_func(k); p++){
            dist += pow(cov_data(i,col_id(k,p)-1) - cov_data(j,col_id(k,p)-1),2);
          }
          dist= pow(dist,0.5);
          
          if(func_def(k)==1){
            if(dist==0){
              val = val*pow(gamma(gamma_idx),2);
            } else {
              val = 0;
            }
          } else if(func_def(k)==2){
            val = val*fexp(dist,gamma(gamma_idx),gamma(gamma_idx+1));
          }else if(func_def(k)==3){
            val = val*pow(gamma(gamma_idx),dist);
          } else if(func_def(k)==4){
            val = val*sqexp(dist,gamma(gamma_idx),gamma(gamma_idx+1));
          } else if(func_def(k)==5){
            val = val*matern(dist,gamma(gamma_idx),gamma(gamma_idx+1));
          } else if(func_def(k)==6){
            val = val*bessel1(dist,gamma(gamma_idx));
          }
          gamma_idx += N_par(k);      
        }
        
        D(i,j) = val;
        D(j,i) = val;
      }
    }
  }
  
  
  for(arma::uword i=0;i<N_dim;i++){
    double val = 1;
    size_t gamma_idx_ii = 0;
    for(arma::uword k=0;k<N_func;k++){
      if(func_def(k)==1){
        val = val*pow(gamma(gamma_idx_ii),2);
      } else if(func_def(k)==2 || func_def(k)==4){
        val = val*gamma(gamma_idx_ii);
      }
    }
    D(i,i) = val;
  }
  
  return D;
}

//' Generates the covariance matrix of the random effects
//' 
//' Generates the covariance matrix of the random effects from a sparse representation
//' @param B Integer specifying the number of blocks in the matrix
//' @param N_dim Vector of integers, which each value specifying the dimension of each block
//' @param N_func Vector of integers specifying the number of functions in the covariance function 
//' for each block.
//' @param func_def Matrix of integers where each column specifies the function definition for each function in each block. 
//' @param N_var_func Matrix of integers of same size as `func_def` with each column specying the number 
//' of variables in the argument to each function in each block
//' @param col_id 3D array (cube) of integers of dimension length(func_def) x max(N_var_func) x B 
//' where each slice the respective column indexes of `cov_data` for each function in the block
//' @param N_par Matrix of integers of same size as `func_def` with each column specifying the number
//' of parameters in the function in each block
//' @param cov_data 3D array (cube) holding the data for the covariance matrix where each of the B slices
//' is the data required for each block
//' @param gamma Vector of covariance parameters specified in order they appear column wise in the functions 
//' specified by `func_def`
//' @return A symmetric positive definite covariance matrix
// [[Rcpp::export]]
arma::field<arma::mat> genD(const arma::uword &B,
                            const arma::uvec &N_dim,
                            const arma::uvec &N_func,
                            const arma::umat &func_def,
                            const arma::umat &N_var_func,
                            const arma::ucube &col_id,
                            const arma::umat &N_par,
                            const arma::uword &sum_N_par,
                            const arma::cube &cov_data,
                            const arma::vec &gamma){
  arma::field<arma::mat> DBlocks(B);
  arma::uword g_idx = 0;
  arma::uword sumpar;
  for(arma::uword b=0;b<B;b++){
    sumpar = sum(N_par.row(b));
    DBlocks[b] = genBlockD(N_dim(b),
                           N_func(b),
                           func_def.row(b).t(),
                           N_var_func.row(b).t(),
                           col_id.slice(b),
                           N_par.row(b).t(),
                           cov_data.slice(b),
                           gamma.subvec(g_idx,g_idx+sumpar-1));
    g_idx += sumpar;
  }
  return(DBlocks);
}

//' Efficiently calculate the inverse of a sub-matrix
//' 
//' Efficiently calculate the inverse of a sub-matrix
//' @details
//' For a given matrix \eqn{A = B^-1}, this will calculate the inverse of a submatrix of B 
//' given only matrix A and requiring only vector multiplication. B typically represents a 
//' covariance matrix for observations 1,...,N and the function is used to calculate the 
//' inverse covariance matrix for a subset of those observations.
//' @param A A square inverse matrix
//' @param i Vector of integers specifying the rows/columns to remove from B, see details
//' @return The inverse of a submatrix of B
//' @examples
//' B <- matrix(runif(16),nrow=4,ncol=4)
//' diag(B) <- diag(B)+1
//' A <- solve(B)
//' remove_one_many_mat(A,1)
//' #equal to
//' solve(B[2:4,2:4])
// [[Rcpp::export]]
arma::mat remove_one_many_mat(const arma::mat &A, 
                              const arma::uvec &i) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  arma::uvec isort = arma::sort(i,"descend"); // sort descending to avoid row conflicts
  arma::mat A2 = A;
  
  for(arma::uword j=0; j<i.n_elem; j++){
    arma::vec idx_new = arma::linspace(0, uidx.n_elem - 1, uidx.n_elem);
    arma::uvec uidx_new = arma::conv_to<arma::uvec>::from(idx_new);
    uidx.shed_row(isort(j));
    uidx_new.shed_row(isort(j));
    double d = A2(isort(j), isort(j));
    arma::vec b = A2.submat(uidx_new, arma::uvec({isort(j)}));
    A2 = A2.submat(uidx_new, uidx_new) - b * b.t() / d;
  }
  
  return A2;
}

double remove_one_many(const arma::mat &A, 
                       const arma::uvec &i,
                       const arma::vec &u) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  arma::uvec isort = arma::sort(i,"descend"); // sort descending to avoid row conflicts
  arma::mat A2 = A;
  for(arma::uword j=0; j<i.n_elem; j++){
    arma::vec idx_new = arma::linspace(0, uidx.n_elem - 1, uidx.n_elem);
    arma::uvec uidx_new = arma::conv_to<arma::uvec>::from(idx_new);
    uidx.shed_row(isort(j));
    uidx_new.shed_row(isort(j));
    double d = A2(isort(j), isort(j));
    arma::vec b = A2.submat(uidx_new, arma::uvec({isort(j)}));
    A2 = A2.submat(uidx_new, uidx_new) - b * b.t() / d;
  }
  return obj_fun(A2, u(uidx));
}


double add_one(const arma::mat &A, 
               double sigma_jj, 
               const arma::vec &f,
               const arma::vec &u) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;
  
  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
    ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));
  
  // step 4: compute K3_inv
  A2 = A2 -
    ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));
  
  return obj_fun(A2, u);
}

//' Efficiently calculates the inverse of a super-matrix 
//' 
//' Efficiently calculates the inverse of a super-matrix by adding an observation
//' @details
//' Given matrix \eqn{A = B^-1} where B is a submatrix of a matrix C, this function efficiently calculates
//' the inverse the matrix B+, which is B adding another row/column from C. For example, if C
//' is the covariance matrix of observations 1,...,N, and B the covariance matrix of observations
//' 1,...,n where n<N, then this function will calculate the inverse of the covariance matrix of 
//' the observations 1,...,n+1.
//' @param A Inverse matrix of dimensions `n`
//' @param sigma_jj The element of C corresponding to the element [j,j] in matrix C where j is the index
//' of the row/column we want to add, see Details
//' @param f A vector of dimension n x 1, corresponding to the elements [b,j] in C where b is the indexes
//' that make up submatrix B
//' @return A matrix of size dim(A)+1
//' @examples
//' B <- matrix(runif(16),nrow=4,ncol=4)
//' diag(B) <- diag(B)+1
//' A <- solve(B[2:4,2:4])
//' add_one_mat(A,B[1,1],B[2:4,1])
//' #equal to
//' solve(B)
// [[Rcpp::export]]
arma::mat add_one_mat(const arma::mat &A, 
                      double sigma_jj, 
                      const arma::vec &f) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  // step 1: compute A*
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  for (arma::uword j = 0; j < A2.n_rows - 1; j++) {
    A2(j, A2.n_cols - 1) = 0;
    A2(A2.n_rows - 1, j) = 0;
  }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;
  
  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
    ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));
  
  // step 4: compute K3_inv
  A2 = A2 -
    ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));
  
  return A2;
}

arma::uvec match_uvec(arma::uvec x, arma::uword val){
  arma::uvec match_vec(x.n_elem, fill::value(-1));
  arma::uword count = 0;
  for(arma::uword j=0; j<x.n_elem; j++){
    if(x(j)==val){
      match_vec(count) = j;
      count += 1;
    }
  }
  return match_vec.rows(0,count-1);
}


arma::uvec std_setdiff(arma::uvec &x, arma::uvec &y) {
  std::vector<int> a = arma::conv_to<std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to<std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}


arma::uvec uvec_minus(const arma::uvec &v, arma::uword rm_idx) {
  arma::uword n = v.size();
  if (rm_idx == 0) return v.tail(n-1);
  if (rm_idx == n-1) return v.head(n-1);
  arma::uvec res(v.size()-1);
  res.head(rm_idx) = v.head(rm_idx);
  res.tail(n-1-rm_idx) = v.tail(n-1-rm_idx);
  return res;
}

#endif