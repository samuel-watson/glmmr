#ifndef GLMMRMATH_H
#define GLMMRMATH_H

#include <cmath>  // std::pow
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double log_factorial_approx(int n){
  double ans;
  if(n==0){
    ans = 0;
  } else {
    ans = n*log(n) - n + log(n*(1+4*n*(1+2*n)))/6 + log(arma::datum::pi)/2;
  }
  return ans;
}

// [[Rcpp::export]]
double gaussian_cdf(const double& x){
  double k = 1.0/(1.0 + 0.2316419*x);
  double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
  
  if (x >= 0.0) {
    return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
  } else {
    return 1.0 - gaussian_cdf(-x);
  }
}

// [[Rcpp::export]]
double gaussian_pdf(const double& x){
  return (1.0/sqrt(2.0 * arma::datum::pi)) * exp(-0.5*x*x);
}

// [[Rcpp::export]]
double log_likelihood(arma::vec y,
                      arma::vec mu,
                      double var_par,
                      std::string family,
                      std::string link) {
  // generate the log-likelihood function
  // for a range of models
  
  double logl = 0;
  arma::uword n = y.n_elem;
  
  if(family=="gaussian"){
    if(link=="identity"){
      for(arma::uword i=0; i<n; i++){
        logl += -0.5*log(var_par) -0.5*log(2*arma::datum::pi) -
          0.5*pow((y(i) - mu(i))/var_par,2);
      }
    }
  }
  if(family=="binomial"){
    if(link=="logit"){
      for(arma::uword i=0; i<n; i++){
        if(y(i)==1){
          logl += log(1/(1+exp(-mu[i])));
        } else if(y(i)==0){
          logl += log(1 - 1/(1+exp(-mu[i])));
        }
      }
    }
    if(link=="log"){
      for(arma::uword i=0; i<n; i++){
        if(y(i)==1){
          logl += mu(i);
        } else if(y(i)==0){
          logl += log(1 - exp(mu(i)));
        }
      }
    }
  }
  if(family=="poisson"){
    if(link=="log"){
      for(arma::uword i=0;i<n; i++){
        double lf1 = log_factorial_approx(y[i]);
        logl += y(i)*mu(i) - exp(mu(i))-lf1;
      }
    }
  }
  
  return logl;
}


// [[Rcpp::export]]
arma::vec mod_inv_func(arma::vec mu,
                       std::string link){
  arma::uword n = mu.n_elem;
  if(link=="logit"){
    for(arma::uword j = 0; j<n; j++){
      mu(j) = exp(mu(j))/(1+exp(mu(j)));
    }
  }
  if(link=="log"){
    for(arma::uword j = 0; j<n; j++){
      mu(j) = exp(mu(j));
    }
  }
  if(link=="probit"){
    for(arma::uword j = 0; j<n; j++){
      mu(j) = gaussian_cdf(mu(j));
    }
  }
    
    return mu;
  
}

// [[Rcpp::export]]
arma::vec gen_dhdmu(arma::vec xb,
                    std::string family,
                    std::string link){
  
  arma::vec wdiag(xb.n_elem, fill::value(1));
  arma::vec p(xb.n_elem, fill::zeros);
  
  if(family=="poisson"){
    if(link=="log"){
      for(arma::uword j; j<xb.n_elem;j++){
        wdiag(j) = 1/exp(xb(j));
      }
    } else if(link =="identity"){
      wdiag = exp(xb);
    }
  } else if(family=="binomial"){
    p = mod_inv_func(xb,"logit");
    if(link=="logit"){
      for(arma::uword j; j<xb.n_elem;j++){
        wdiag(j) = 1/(p(j)*(1-p(j)));
      }
    } else if(link=="log"){
      for(arma::uword j; j<xb.n_elem;j++){
        wdiag(j) = (1-p(j))/p(j);
      }
    } else if(link=="identity"){
      for(arma::uword j; j<xb.n_elem;j++){
        wdiag(j) = (p(j)*(1-p(j)));
      }
    } else if(link=="probit"){
      p = mod_inv_func(xb,"probit");
      arma::vec p2(xb.n_elem,fill::zeros);
      for(arma::uword j; j<xb.n_elem;j++){
        wdiag(j) = (p(j)*(1-p(j)))/gaussian_pdf(xb(j));
      }
    }
  } else if(link=="gaussian"){
    // if identity do nothin
    if(link=="log"){
      for(arma::uword j; j<xb.n_elem;j++){
        wdiag(j) = 1/exp(xb(j));
      }
    }
  } // for gamma- inverse do nothing
  return wdiag;
}

// [[Rcpp::export]]
double obj_fun(const arma::mat &A, const arma::vec &U2) {
  // this is the directional derivative
  return arma::as_scalar(U2.t() * A * U2);
}

// [[Rcpp::export]]
double c_obj_fun(arma::mat M, arma::vec C) {
  // this is the objective function c-optimal
  arma::mat M_inv = arma::inv_sympd(M);
  return arma::as_scalar(C.t() * M_inv * C);
}

// [[Rcpp::export]]
arma::mat gen_m(const arma::mat &X, const arma::mat &A) {
  //generate information matrix
  return X.t() * A * X;
}

#endif