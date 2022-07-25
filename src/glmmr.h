#ifndef GLMMR_H
#define GLMMR_H

#include <cmath> 
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// class dmat{
//   arma::uword B_;
//   arma::uvec N_dim_;
//   arma::uvec N_func_;
//   arma::umat func_def_;
//   arma::umat N_var_func_;
//   arma::ucube col_id_;
//   arma::umat N_par_;
//   arma::uword sum_N_par_;
//   arma::cube cov_data_;
//   arma::vec gamma_;
// public :
//   dmat (const arma::uword &B,
//         const arma::uvec &N_dim,
//         const arma::uvec &N_func,
//         const arma::umat &func_def,
//         const arma::umat &N_var_func,
//         const arma::ucube &col_id,
//         const arma::umat &N_par,
//         const arma::uword &sum_N_par,
//         const arma::cube &cov_data,
//         const arma::vec &gamma) :
//   B_(B), N_dim_(N_dim), 
//   N_func_(N_func), 
//   func_def_(func_def), N_var_func_(N_var_func),
//   col_id_(col_id), N_par_(N_par), sum_N_par_(sum_N_par),
//   cov_data_(cov_data), gamma_(gamma) {}
//   arma::field<arma::mat>
//   
// };

inline double gaussian_cdf(double x){
  return R::pnorm(x, 0, 1, true, false);
}

inline arma::vec gaussian_cdf_vec(const arma::vec& v){
  arma::vec res = arma::zeros<arma::vec>(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i)
    res[i] = gaussian_cdf(v[i]);
  return res;
}

inline double gaussian_pdf(double x){
  return R::dnorm(x, 0, 1, false);
}

inline arma::vec gaussian_pdf_vec(const arma::vec& v){
  arma::vec res = arma::zeros<arma::vec>(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i)
    res[i] = gaussian_pdf(v[i]);
  return res;
}

//' Log multivariate Gaussian probability density funciton
//' 
//' Log multivariate Gaussian probability density funciton
//' @param u Vector of realisations from the distribution
//' @param D Inverse covariance matrix
//' @param logdetD Log determinant of the covariance matrix
// [[Rcpp::export]]
inline double log_mv_gaussian_pdf(const arma::vec& u,
                                  const arma::mat& D,
                                  const double& logdetD){
  arma::uword Q = u.n_elem;
  return (-0.5*Q*log(2*arma::datum::pi)-
          0.5*logdetD - 0.5*arma::as_scalar(u.t()*D*u));
}

// //' Exponential covariance function
// //' 
// //' Exponential covariance function
// //' @details
// //' \deqn{\theta_1 exp(-x/\theta_2)}
// //' @param x Numeric value 
// //' @param par1 First parameter of the distribution
// // [[Rcpp::export]]
// inline double fexp(const double &x, 
//                    double par1) {
//   return exp(-1*x/par1);
// }

//' Squared exponential covariance function
//' 
//' Squared exponential covariance function
//' @details
//' \deqn{\theta_1 exp(-x^2/\theta_2^2)}
//' @param x Numeric value 
//' @param par1 First parameter of the distribution
//' @param par2 Second parameter of the distribution
// [[Rcpp::export]]
inline double sqexp(const double &x, 
                    double par1,
                    double par2) {
  return par1*exp(-1*pow(x,2)/pow(par2,2));
}

//' Matern covariance function
//' 
//' Matern covariance function
//' @details
//' TBC
//' @param x Numeric value 
//' @param rho First parameter of the distribution
//' @param nu Second parameter of the distribution
// [[Rcpp::export]]
inline double matern(const double &x,
                     double rho, 
                     double nu){
  double xr = pow(2*nu,0.5)*x/rho;
  double ans = 1;
  if(xr!=0){
    if(nu == 0.5){
      ans = exp(-xr);
    } else {
      double cte = pow(2,-1*(nu-1))/R::gammafn(nu);
      ans = cte*pow(xr, nu)*R::bessel_k(xr,nu,1);
    }
  }
  return ans;
}

//' Bessel covariance function
//' 
//' Bessel covariance function
//' @details
//' TBC
//' @param x Numeric value 
//' @param rho First parameter of the distribution
// [[Rcpp::export]]
inline double bessel1(const double &x,
                      double rho){
  double xr = x/rho;
  return xr* R::bessel_k(xr,1,1);
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
inline arma::mat genBlockD(size_t N_dim,
                    size_t N_func,
                    const arma::uvec &func_def,
                    const arma::uvec &N_var_func,
                    const arma::umat &col_id,
                    const arma::uvec &N_par,
                    const arma::mat &cov_data,
                    const arma::vec &gamma){
  arma::mat D(N_dim,N_dim,fill::zeros);
  if(!all(func_def == 1)){
#pragma omp parallel for
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
            val = val*exp(-1*dist/gamma(gamma_idx));//fexp(dist,gamma(gamma_idx));
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
  
#pragma omp parallel for
  for(arma::uword i=0;i<N_dim;i++){
    double val = 1;
    size_t gamma_idx_ii = 0;
    for(arma::uword k=0;k<N_func;k++){
      if(func_def(k)==1){
        val = val*pow(gamma(gamma_idx_ii),2);
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
inline arma::field<arma::mat> genD(const arma::uword &B,
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

#endif