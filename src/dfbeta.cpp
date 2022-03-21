#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

  // [[Rcpp::export]]
arma::rowvec dfbeta_stat(const arma::mat &sigma,
                       const arma::mat &X,
                       const arma::vec &y,
                       arma::uword par){
  //initialise the matrices
  arma::mat invS = arma::inv_sympd(sigma);
  arma::mat invSx = invS * X;
  arma::mat M = X.t() * invSx;
  arma::mat B = M.i() * invSx.t();
  invS += (-1)*invSx*B;
  arma::mat dQ = diagmat(invS.diag());
  dQ = inv_sympd(dQ);
  arma::mat e = diagmat(dQ * invS.t() * y);
  B = B*e;
  
  return B.row(par-1);
}
