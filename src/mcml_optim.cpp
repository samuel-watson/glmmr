#include <cmath>  // std::pow
#include <RcppArmadillo.h>
#include "rbobyqa.h"
#include "glmmrmath.h"
#include "glmmrmatrix.h"
using namespace rminqa;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

class D_likelihood : public ObjFun {
  arma::uword B_;
  arma::uvec N_dim_;
  arma::uvec N_func_;
  arma::umat func_def_;
  arma::umat N_var_func_;
  arma::ucube col_id_;
  arma::umat N_par_;
  arma::uword sum_N_par_;
  arma::cube cov_data_;
  arma::mat u_;
  
public:
  D_likelihood(const arma::uword &B,
               const arma::uvec &N_dim,
               const arma::uvec &N_func,
               const arma::umat &func_def,
               const arma::umat &N_var_func,
               const arma::ucube &col_id,
               const arma::umat &N_par,
               const arma::uword &sum_N_par,
               const arma::cube &cov_data, 
               const arma::mat &u) : 
    B_(B), N_dim_(N_dim), 
    N_func_(N_func), 
    func_def_(func_def), N_var_func_(N_var_func),
    col_id_(col_id), N_par_(N_par), sum_N_par_(sum_N_par),
    cov_data_(cov_data), u_(u) {}
  double operator()(const vec &par) override{
    arma::uword nrow = u_.n_cols;
    arma::field<arma::mat> Dfield = genD(B_,N_dim_,
                                         N_func_,
                                         func_def_,N_var_func_,
                                         col_id_,N_par_,sum_N_par_,
                                         cov_data_,par);
    double dmv = 0;
    double logdetD;
    for(arma::uword b=0;b<B_;b++){
      arma::uword ndim_idx = 0;
      logdetD = arma::log_det_sympd(Dfield[b]);
      for(arma::uword j=0;j<nrow;j++){
        dmv += log_mv_gaussian_pdf(u_.col(j).subvec(ndim_idx,ndim_idx+N_dim_(b)-1),
                                   inv_sympd(Dfield[b]),logdetD);
      }
      ndim_idx += N_dim_(b);
    }
    
    return -1 * dmv/nrow;
  }
};


// [[Rcpp::export]]
arma::vec d_lik_optim(const arma::uword &B,
                      const arma::uvec &N_dim,
                      const arma::uvec &N_func,
                      const arma::umat &func_def,
                      const arma::umat &N_var_func,
                      const arma::ucube &col_id,
                      const arma::umat &N_par,
                      const arma::uword &sum_N_par,
                      const arma::cube &cov_data, 
                      const arma::mat &u,
                      arma::vec start,
                      const arma::vec &lower,
                      const arma::vec &upper,
                      int trace = 0){
  D_likelihood dl(B,N_dim,
                  N_func,
                  func_def,N_var_func,
                  col_id,N_par,sum_N_par,
                  cov_data, u);
  
  Rbobyqa<D_likelihood> opt;
  opt.set_upper(upper);
  opt.set_lower(lower);
  opt.control.iprint = trace;
  opt.minimize(dl, start);
  return opt.par();
}




class L_likelihood : public ObjFun {
  arma::mat Z_;
  arma::mat X_;
  arma::vec y_;
  arma::mat u_;
  std::string family_;
  std::string link_;
  
public:
  L_likelihood(const arma::mat &Z, 
               const arma::mat &X,
               const arma::vec &y, 
               const arma::mat &u, 
               std::string family, 
               std::string link) : 
  Z_(Z), X_(X), y_(y), u_(u),family_(family) , link_(link) {}
  double operator()(const vec &par) override{
    arma::uword niter = u_.n_cols;
    arma::uword n = y_.n_elem;
    arma::vec zd(n);
    arma::uword P = par.n_elem;
    arma::vec xb(n);
    double var_par;
    
    if(family_=="gaussian"){
      var_par = par(P-1);
      xb = X_*par.subvec(0,P-2);
    } else {
      var_par = 0;
      xb = X_*par;
    }
    
    double lfa;
    
    for(arma::uword j=0; j<niter ; j++){
      zd = Z_ * u_.col(j);
      double ll = log_likelihood(y_,
                                 xb + zd,
                                 var_par,
                                 family_,
                                 link_);
      lfa += ll;
    }
    
    return -1 * lfa/niter;
  }
};

// [[Rcpp::export]]
arma::vec l_lik_optim(const arma::mat &Z, 
                      const arma::mat &X,
                      const arma::vec &y, 
                      const arma::mat &u, 
                      std::string family, 
                      std::string link,
                      arma::vec start,
                      const arma::vec &lower,
                      const arma::vec &upper,
                      int trace){
  L_likelihood dl(Z,X,y,u,family,link);
  
  Rbobyqa<L_likelihood> opt;
  opt.set_upper(upper);
  opt.set_lower(lower);
  opt.control.iprint = trace;
  
  // Roptim<L_likelihood> opt("Nelder-Mead");
  // opt.control.trace = 0;
  //opt.set_hessian(true);
  opt.minimize(dl, start);
  
  return opt.par();
}


class F_likelihood : public ObjFun {
  arma::uword B_;
  arma::uvec N_dim_;
  arma::uvec N_func_;
  arma::umat func_def_;
  arma::umat N_var_func_;
  arma::ucube col_id_;
  arma::umat N_par_;
  arma::uword sum_N_par_;
  arma::cube cov_data_;
  arma::mat Z_;
  arma::mat X_;
  arma::vec y_;
  arma::mat u_;
  arma::vec cov_par_fix_;
  std::string family_;
  std::string link_;
  bool importance_;
  
public:
  F_likelihood(const arma::uword &B,
               const arma::uvec &N_dim,
               const arma::uvec &N_func,
               const arma::umat &func_def,
               const arma::umat &N_var_func,
               const arma::ucube &col_id,
               const arma::umat &N_par,
               const arma::uword &sum_N_par,
               const arma::cube &cov_data,
               arma::mat Z, 
               arma::mat X,
               arma::vec y, 
               arma::mat u,
               arma::vec cov_par_fix,
               std::string family, 
               std::string link,
               bool importance) : 
  B_(B), N_dim_(N_dim), 
  N_func_(N_func), 
  func_def_(func_def), N_var_func_(N_var_func),
  col_id_(col_id), N_par_(N_par), sum_N_par_(sum_N_par),
  cov_data_(cov_data),Z_(Z), X_(X), y_(y), 
  u_(u),cov_par_fix_(cov_par_fix), family_(family) , link_(link),
  importance_(importance){}
  double operator()(const vec &par) override{
    arma::uword niter = u_.n_cols;
    arma::uword n = y_.n_elem;
    arma::uword P = X_.n_cols;
    arma::uword Q = cov_par_fix_.n_elem;
    double du;
    
    //now loop through all the iterations
    // log likelihood of D
    arma::field<arma::mat> Dfield = genD(B_,N_dim_,
                                         N_func_,
                                         func_def_,N_var_func_,
                                         col_id_,N_par_,sum_N_par_,
                                         cov_data_,par.subvec(P,P+Q-1));
    arma::vec numerD(niter,fill::zeros);
    double logdetD;
    for(arma::uword b=0;b<B_;b++){
      arma::uword ndim_idx = 0;
      logdetD = arma::log_det_sympd(Dfield[b]);
      for(arma::uword j=0;j<niter;j++){
        numerD(j) += log_mv_gaussian_pdf(u_.col(j).subvec(ndim_idx,ndim_idx+N_dim_(b)-1),
               inv_sympd(Dfield[b]),logdetD);
      }
      ndim_idx += N_dim_(b);
    }
    
    // log likelihood for observations
    arma::vec zd(n);
    arma::vec xb(n);
    double var_par;
    
    if(family_=="gaussian"){
      var_par = par(P+Q);
    } else {
      var_par = 0;
    }
    
    xb = X_*par.subvec(0,P-1);
    arma::vec lfa(niter,fill::zeros);
    
    for(arma::uword j=0; j<niter ; j++){
      zd = Z_ * u_.col(j);
      lfa(j) += log_likelihood(y_,
          xb + zd,
          var_par,
          family_,
          link_);
    }
    
    if(importance_){
      // denominator density for importance sampling
      Dfield = genD(B_,N_dim_,
                    N_func_,
                    func_def_,N_var_func_,
                    col_id_,N_par_,sum_N_par_,
                    cov_data_,cov_par_fix_);
      arma::vec denomD(niter,fill::zeros);
      for(arma::uword b=0;b<B_;b++){
        arma::uword ndim_idx = 0;
        logdetD = arma::log_det_sympd(Dfield[b]);
        for(arma::uword j=0;j<niter;j++){
          denomD(j) += log_mv_gaussian_pdf(u_.col(j).subvec(ndim_idx,ndim_idx+N_dim_(b)-1),
                 inv_sympd(Dfield[b]),logdetD);
        }
        ndim_idx += N_dim_(b);
      }
      du = 0;
      for(arma::uword j=0;j<niter;j++){
        du  += exp(lfa(j)+numerD(j))/exp(denomD(j));
      }
      
      du = -1 * log(du/niter);
    } else {
      du = -1* (mean(numerD) + mean(lfa));
    }
    
    return du;
  }
};

// [[Rcpp::export]]
arma::vec f_lik_grad(const arma::uword &B,
                     const arma::uvec &N_dim,
                     const arma::uvec &N_func,
                     const arma::umat &func_def,
                     const arma::umat &N_var_func,
                     const arma::ucube &col_id,
                     const arma::umat &N_par,
                     const arma::uword &sum_N_par,
                     const arma::cube &cov_data,
                     const arma::mat &Z, 
                     const arma::mat &X,
                     const arma::vec &y, 
                     const arma::mat &u,
                     const arma::vec &cov_par_fix,
                     std::string family, 
                     std::string link,
                     arma::vec start,
                     const arma::vec &lower,
                     const arma::vec &upper,
                     double tol){
  
  F_likelihood dl(B,N_dim,
                  N_func,
                  func_def,N_var_func,
                  col_id,N_par,sum_N_par,
                  cov_data,Z,X,y,u,
                  cov_par_fix,family,
                  link,false);
  
  dl.os.usebounds_ = 1;
  if(!lower.is_empty()){
    dl.os.lower_ = lower;
  }
  if(!upper.is_empty()){
    dl.os.upper_ = upper;
  }
  dl.os.ndeps_ = arma::ones<arma::vec>(start.size()) * tol;
  arma::vec gradient(start.n_elem,fill::zeros);
  dl.Gradient(start,gradient);
  return gradient;
}

// [[Rcpp::export]]
arma::mat f_lik_hess(const arma::uword &B,
                     const arma::uvec &N_dim,
                     const arma::uvec &N_func,
                     const arma::umat &func_def,
                     const arma::umat &N_var_func,
                     const arma::ucube &col_id,
                     const arma::umat &N_par,
                     const arma::uword &sum_N_par,
                     const arma::cube &cov_data,
                     const arma::mat &Z,
                     const arma::mat &X,
                     const arma::vec &y,
                     const arma::mat &u,
                     const arma::vec &cov_par_fix,
                     std::string family,
                     std::string link,
                     arma::vec start,
                     const arma::vec &lower,
                     const arma::vec &upper,
                     double tol = 1e-4){
  F_likelihood dl(B,N_dim,
                  N_func,
                  func_def,N_var_func,
                  col_id,N_par,sum_N_par,
                  cov_data,Z,X,y,u,
                  cov_par_fix,family,
                  link,false);
  dl.os.usebounds_ = 1;
  if(!lower.is_empty()){
    dl.os.lower_ = lower;
  }
  if(!upper.is_empty()){
    dl.os.upper_ = upper;
  }
  dl.os.ndeps_ = arma::ones<arma::vec>(start.size()) * tol;
  arma::mat hessian(start.n_elem,start.n_elem,fill::zeros);
  dl.Hessian(start,hessian);
  return hessian;
}

// [[Rcpp::export]]
arma::mat f_lik_optim(const arma::uword &B,
                      const arma::uvec &N_dim,
                      const arma::uvec &N_func,
                      const arma::umat &func_def,
                      const arma::umat &N_var_func,
                      const arma::ucube &col_id,
                      const arma::umat &N_par,
                      const arma::uword &sum_N_par,
                      const arma::cube &cov_data,
                      const arma::mat &Z, 
                      const arma::mat &X,
                      const arma::vec &y, 
                      const arma::mat &u,
                      const arma::vec &cov_par_fix,
                      std::string family, 
                      std::string link,
                      arma::vec start,
                      const arma::vec &lower,
                      const arma::vec &upper,
                      int trace){
  
  F_likelihood dl(B,N_dim,
                  N_func,
                  func_def,N_var_func,
                  col_id,N_par,sum_N_par,
                  cov_data,Z,X,y,u,
                  cov_par_fix,family,
                  link,true);
  
  Rbobyqa<F_likelihood> opt;
  opt.set_lower(lower);
  opt.control.iprint = trace;
  opt.minimize(dl, start);
  
  return opt.par();
}



// [[Rcpp::export]]
Rcpp::List mcnr_step(const arma::vec &y, const arma::mat &X, const arma::mat &Z,
                     const arma::vec &beta, const arma::mat &u,
                     const std::string &family, const std::string &link){
  const arma::uword n = y.n_elem;
  const arma::uword P = X.n_cols;
  const arma::uword niter = u.n_cols;
  
  //generate residuals
  arma::vec xb = X*beta;
  arma::vec sigmas(niter);
  arma::mat XtWX(P,P,fill::zeros);
  arma::vec Wu(n,fill::zeros);
  arma::vec zd(n,fill::zeros);
  arma::vec resid(n,fill::zeros);
  arma::vec wdiag(n,fill::zeros);
  arma::vec wdiag2(n,fill::zeros);
  arma::mat W(n,n,fill::zeros);
  arma::mat W2(n,n,fill::zeros);
  for(arma::uword i = 0; i < niter; ++i){
    zd = Z * u.col(i);
    zd = mod_inv_func(xb + zd, link);
    resid = y - zd;
    sigmas(i) = arma::stddev(resid);
    wdiag = gen_dhdmu(xb + zd,family,link);
    wdiag2 = 1/arma::pow(wdiag, 2);
    if(family=="gaussian" || family=="gamma") wdiag2 *= sigmas(i);
    for(arma::uword j = 0; j<n; j++){
      XtWX += wdiag2(j)*(X.row(j).t()*X.row(j));
      Wu(j) += wdiag(j)*wdiag2(j)*resid(j); 
    }
  }
  //Rcpp::Rcout<< XtWX;
  XtWX = arma::inv_sympd(XtWX/niter);
  Wu = Wu/niter;
  
  arma::vec bincr = XtWX*X.t()*Wu;
  Rcpp::List L = List::create(_["beta_step"] = bincr , _["sigmahat"] = arma::mean(sigmas));
  return L;
}

// [[Rcpp::export]]
double aic_mcml(const arma::mat &Z, 
                const arma::mat &X,
                const arma::vec &y, 
                const arma::mat &u, 
                std::string family, 
                std::string link,
                const arma::uword &B,
                const arma::uvec &N_dim,
                const arma::uvec &N_func,
                const arma::umat &func_def,
                const arma::umat &N_var_func,
                const arma::ucube &col_id,
                const arma::umat &N_par,
                const arma::uword &sum_N_par,
                const arma::cube &cov_data,
                const arma::vec& beta_par,
                const arma::vec& cov_par){
  arma::uword niter = u.n_cols;
  arma::uword n = y.n_elem;
  arma::vec zd(n);
  arma::uword P = beta_par.n_elem;
  arma::vec xb(n);
  double var_par;
  arma::uword dof = beta_par.n_elem + cov_par.n_elem;
  
  if(family=="gaussian"){
    var_par = beta_par(P-1);
    xb = X*beta_par.subvec(0,P-2);
  } else {
    var_par = 0;
    xb = X*beta_par;
  }
  
  arma::field<arma::mat> Dfield = genD(B,N_dim,
                                       N_func,
                                       func_def,N_var_func,
                                       col_id,N_par,sum_N_par,
                                       cov_data,cov_par);
  double dmv = 0;
  double logdetD;
  for(arma::uword b=0;b<B;b++){
    arma::uword ndim_idx = 0;
    logdetD = arma::log_det_sympd(Dfield[b]);
    for(arma::uword j=0;j<niter;j++){
      dmv += log_mv_gaussian_pdf(u.col(j).subvec(ndim_idx,ndim_idx+N_dim(b)-1),
                                 inv_sympd(Dfield[b]),logdetD);
    }
    ndim_idx += N_dim(b);
  }
  
  double lfa = 0;
  
  for(arma::uword j=0; j<niter ; j++){
    zd = Z * u.col(j);
    lfa += log_likelihood(y,
                          xb + zd,
                          var_par,
                          family,
                          link);
  }
  
  return (-2*( lfa/niter + dmv/niter ) + 2*arma::as_scalar(dof)); 
  
}

