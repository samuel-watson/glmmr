#ifndef GLMMRMATRIX_H
#define GLMMRMATRIX_H

#include <RcppArmadillo.h>
#include "glmmrmath.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double findFunc(arma::vec funcdetail,
                arma::rowvec data1,
                arma::rowvec data2,
                arma::vec pars){
  //generate Euclidean distance
  double val = 0;
  arma::uword f1 = funcdetail(0);
  
  for(arma::uword j=0;j<4;j++){
    arma::sword idx = funcdetail(j+1);
    if(idx >= 0){
      val += pow(data1(idx) - data2(idx),2);
    }
  }
  val = sqrt(val);
  double out;
  //group
  if(f1==1){
    if(val==0){
      out = pow(pars(funcdetail(5)),2);
    } else {
      out = 0;
    }
  }
  
  //exponential 1
  if(f1==2){
    out = fexp(val,pars(funcdetail(5)),pars(funcdetail(6)));
  }
  
  //ar1
  if(f1==3){
    out = pow(pars(funcdetail(5)),val);
  }
  //c("gr","fexp","ar1","sqexp","matern","bessel",)
  //squared exponential
  if(f1==4){
    out = sqexp(val,pars(funcdetail(5)),pars(funcdetail(6)));
  }
  
  if(f1==5){
    out = matern(val,pars(funcdetail(5)),pars(funcdetail(6)));
  }
  
  if(f1==6){
    out = bessel1(val,pars(funcdetail(5)));
  }
  
  return(out);
}

// [[Rcpp::export]]
arma::mat blockMat(arma::mat mat1,
                   arma::mat mat2){
  arma::uword n1 = mat1.n_rows;
  arma::uword n2 = mat2.n_rows;
  arma::mat dmat(n1+n2,n1+n2);
  dmat.fill(0);
  dmat.submat(0,0,n1-1,n1-1) = mat1;
  dmat.submat(n1,n1,n1+n2-1,n1+n2-1) = mat2;
  return(dmat);
}

// [[Rcpp::export]]
arma::mat genSubD(arma::mat &func,
                  arma::mat &data,
                  arma::vec pars){
  arma::uword nF = func.n_cols;
  arma::uword N = data.n_rows;
  
  //initialise the matrix
  arma::mat D(N,N);
  D.fill(1);
  double val;
  //loop through the matrix and produce all the elements
  for(arma::uword j=0; j<N; j++){
    for(arma::uword k=j+1; k < N; k++){
      for(arma::uword l=0; l<nF; l++){
        val = findFunc(func.col(l),
                       data.row(j),
                       data.row(k),
                       pars);
        D(j,k) = D(j,k)*val;
        D(k,j) = D(k,j)*val;
      }
      
    }
    for(arma::uword l=0; l<nF; l++){
      val = findFunc(func.col(l),
                     data.row(j),
                     data.row(j),
                     pars);
      D(j,j) = D(j,j)*val;
    }
  }
  return(D);
}

// [[Rcpp::export]]
arma::mat genD(Rcpp::List &func,
               Rcpp::List &data,
               arma::vec pars){
  //arma::uword n = func.n_slices;
  int n = func.length();
  arma::mat f1 = func[0];
  arma::mat d1 = data[0];
  arma::mat D1 = genSubD(f1,d1,pars);
  if(n > 1){
    for(int j=1;j<n;j++){
      arma::mat f1 = func[j];
      arma::mat d1 = data[j];
      arma::mat D2 = genSubD(f1,d1,pars);
      D1 = blockMat(D1,D2);
    }
  }
  return(D1);
}

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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
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


// [[Rcpp::export]]
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