#ifndef GLMMRMATH_H
#define GLMMRMATH_H

#include <cmath>  // std::pow

#include <Rmath.h>
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
double gaussian_cdf(double x){
  return R::pnorm(x, 0, 1, true, false);
}

// [[Rcpp::export]]
arma::vec gaussian_cdf_vec(const arma::vec& v){
  arma::vec res = arma::zeros<arma::vec>(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i)
    res[i] = gaussian_cdf(v[i]);
  return res;
}

// [[Rcpp::export]]
double gaussian_pdf(double x){
  return R::dnorm(x, 0, 1, false);
}

// [[Rcpp::export]]
arma::vec gaussian_pdf_vec(const arma::vec& v){
  arma::vec res = arma::zeros<arma::vec>(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i)
    res[i] = gaussian_pdf(v[i]);
  return res;
}

// [[Rcpp::export]]
double log_mv_gaussian_pdf(const arma::rowvec& u,
                           const arma::mat& D,
                           const double& logdetD){
  arma::uword Q = u.n_elem;
  return (-0.5*Q*log(2*arma::datum::pi)-
    0.5*logdetD - 0.5*arma::as_scalar(u*D*u.t()));
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
  //arma::uword n = mu.n_elem;
  if(link=="logit"){
    mu = exp(mu) / (1+exp(mu));
  }
  if(link=="log"){
    mu = exp(mu);
  }
  if(link=="probit"){
    mu = gaussian_cdf_vec(mu);
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
      wdiag = 1/exp(xb);
    } else if(link =="identity"){
      wdiag = exp(xb);
    }
  } else if(family=="binomial"){
    p = mod_inv_func(xb,"logit");
    if(link=="logit"){
      wdiag = 1/(p % (1-p));
    } else if(link=="log"){
      wdiag = (1-p)/p;
    } else if(link=="identity"){
      wdiag = p % (1-p);
    } else if(link=="probit"){
      p = mod_inv_func(xb,"probit");
      arma::vec p2(xb.n_elem,fill::zeros);
      wdiag = (p % (1-p))/gaussian_pdf_vec(xb);
    }
  } else if(link=="gaussian"){
    // if identity do nothin
    if(link=="log"){
      wdiag = 1/exp(xb);
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
