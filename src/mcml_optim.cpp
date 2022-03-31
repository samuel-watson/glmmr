#include <cmath>  
#include <RcppArmadillo.h>
#include "rbobyqa.h"
#include "glmmrmath.h"
#include "glmmrmatrix.h"
using namespace rminqa;
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

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
    arma::vec dmvvec(nrow,fill::zeros);
    double logdetD;
    arma::uword ndim_idx = 0;
    for(arma::uword b=0;b<B_;b++){
      if(all(func_def_.row(b)==1)){
#pragma omp parallel for collapse(2)
        for(arma::uword j=0;j<nrow;j++){
          for(arma::uword k=0; k<Dfield[b].n_rows; k++){
            dmvvec(j) += -0.5*log(Dfield[b](k,k)) -0.5*log(2*arma::datum::pi) -
              0.5*pow(u_(ndim_idx+k,j),2)/Dfield[b](k,k);
          }
        }
        
      } else {
        arma::mat invD = inv_sympd(Dfield[b]);
        logdetD = arma::log_det_sympd(Dfield[b]);
#pragma omp parallel for
        for(arma::uword j=0;j<nrow;j++){
          dmvvec(j) += log_mv_gaussian_pdf(u_.col(j).subvec(ndim_idx,ndim_idx+N_dim_(b)-1),
                 invD,logdetD);
        }
      }
      ndim_idx += N_dim_(b);
    }
    return -1 * mean(dmvvec);
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

// [[Rcpp::export]]
arma::mat d_lik_hess(const arma::uword &B,
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
                      double tol = 1e-4){
  D_likelihood dl(B,N_dim,
                  N_func,
                  func_def,N_var_func,
                  col_id,N_par,sum_N_par,
                  cov_data, u);
  
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
arma::vec d_lik_grad(const arma::uword &B,
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
                     double tol = 1e-4){
  D_likelihood dl(B,N_dim,
                  N_func,
                  func_def,N_var_func,
                  col_id,N_par,sum_N_par,
                  cov_data, u);
  
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
    arma::uword P = par.n_elem;
    arma::vec xb(n);
    double var_par;
    
    if(family_=="gaussian"){
      xb = X_*par.subvec(0,P-2);
      var_par = par(P-1);
    } else {
      xb = X_*par;
      var_par = 0;
    }
    
    // double lfa;
    arma::vec ll(niter,fill::zeros);
    arma::mat zd = Z_ * u_;
#pragma omp parallel for
    for(arma::uword j=0; j<niter ; j++){
      ll(j) += log_likelihood(y_,
                             xb + zd.col(j),
                             var_par,
                             family_,
                             link_);
      
    }
    
    return -1 * mean(ll);
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
  
  opt.minimize(dl, start);
  
  return opt.par();
}

// [[Rcpp::export]]
arma::mat l_lik_hess(const arma::mat &Z, 
                      const arma::mat &X,
                      const arma::vec &y, 
                      const arma::mat &u, 
                      std::string family, 
                      std::string link,
                      arma::vec start,
                      const arma::vec &lower,
                      const arma::vec &upper,
                      double tol = 1e-4){
  L_likelihood dl(Z,X,y,u,family,link);
  
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
arma::vec l_lik_grad(const arma::mat &Z, 
                     const arma::mat &X,
                     const arma::vec &y, 
                     const arma::mat &u, 
                     std::string family, 
                     std::string link,
                     arma::vec start,
                     const arma::vec &lower,
                     const arma::vec &upper,
                     double tol = 1e-4){
  L_likelihood dl(Z,X,y,u,family,link);
  
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
  importance_(importance) {}
  double operator()(const vec &par) override{
    
    arma::uword niter = u_.n_cols;
    arma::uword n = y_.n_elem;
    arma::uword P = X_.n_cols;
    arma::uword Q = cov_par_fix_.n_elem;
    double du;
    
    arma::field<arma::mat> Dfield = genD(B_,N_dim_,
                                         N_func_,
                                         func_def_,N_var_func_,
                                         col_id_,N_par_,sum_N_par_,
                                         cov_data_,par.subvec(P,P+Q-1));
    arma::vec numerD(niter,fill::zeros);
    double logdetD;
    arma::uword ndim_idx = 0;
    for(arma::uword b=0;b<B_;b++){
      if(all(func_def_.row(b)==1)){
#pragma omp parallel for collapse(2)
        for(arma::uword j=0;j<niter;j++){
          for(arma::uword k=0; k<Dfield[b].n_rows; k++){
            numerD(j) += -0.5*log(Dfield[b](k,k)) -0.5*log(2*arma::datum::pi) -
              0.5*pow(u_(ndim_idx+k,j),2)/Dfield[b](k,k);
          }
        }
        
      } else {
        arma::mat invD = inv_sympd(Dfield[b]);
        logdetD = arma::log_det_sympd(Dfield[b]);
#pragma omp parallel for
        for(arma::uword j=0;j<niter;j++){
          numerD(j) += log_mv_gaussian_pdf(u_.col(j).subvec(ndim_idx,ndim_idx+N_dim_(b)-1),
                 invD,logdetD);
        }
      }
      ndim_idx += N_dim_(b);
    }
    
    // log likelihood for observations
    //arma::vec zd(n);
    arma::vec xb(n);
    double var_par;
    
    if(family_=="gaussian"){
      var_par = par(P+Q);
    } else {
      var_par = 0;
    }
    
    xb = X_*par.subvec(0,P-1);
    arma::vec lfa(niter,fill::zeros);
    arma::mat zd = Z_ * u_;
#pragma omp parallel for
    for(arma::uword j=0; j<niter ; j++){
      lfa(j) += log_likelihood(y_,
          xb + zd.col(j),
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
      double logdetD;
      arma::uword ndim_idx = 0;
      for(arma::uword b=0;b<B_;b++){
        if(all(func_def_.row(b)==1)){
#pragma omp parallel for collapse(2)
          for(arma::uword j=0;j<niter;j++){
            for(arma::uword k=0; k<Dfield[b].n_rows; k++){
              denomD(j) += -0.5*log(Dfield[b](k,k)) -0.5*log(2*arma::datum::pi) -
                0.5*pow(u_(ndim_idx+k,j),2)/Dfield[b](k,k);
            }
          }
          
        } else {
          arma::mat invD = inv_sympd(Dfield[b]);
          logdetD = arma::log_det_sympd(Dfield[b]);
#pragma omp parallel for
          for(arma::uword j=0;j<niter;j++){
            denomD(j) += log_mv_gaussian_pdf(u_.col(j).subvec(ndim_idx,ndim_idx+N_dim_(b)-1),
                   invD,logdetD);
          }
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
    //sigmas(i) = pow(sigmas(i),2);
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
  //arma::vec zd(n);
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
  arma::vec dmvvec(niter,fill::zeros);
  double logdetD;
  arma::uword ndim_idx = 0;
  for(arma::uword b=0;b<B;b++){
    if(all(func_def.row(b)==1)){
#pragma omp parallel for collapse(2)
      for(arma::uword j=0;j<niter;j++){
        for(arma::uword k=0; k<Dfield[b].n_rows; k++){
          dmvvec(j) += -0.5*log(Dfield[b](k,k)) -0.5*log(2*arma::datum::pi) -
            0.5*pow(u(ndim_idx+k,j),2)/Dfield[b](k,k);
        }
      }
      
    } else {
      arma::mat invD = inv_sympd(Dfield[b]);
      logdetD = arma::log_det_sympd(Dfield[b]);
#pragma omp parallel for
      for(arma::uword j=0;j<niter;j++){
        dmvvec(j) += log_mv_gaussian_pdf(u.col(j).subvec(ndim_idx,ndim_idx+N_dim(b)-1),
               invD,logdetD);
      }
    }
    ndim_idx += N_dim(b);
  }
  
  arma::vec ll(niter,fill::zeros);
  arma::mat zd = Z * u;
#pragma omp parallel for
  for(arma::uword j=0; j<niter ; j++){
    ll(j) += log_likelihood(y,
       xb + zd.col(j),
       var_par,
       family,
       link);
  }
  
  return (-2*( mean(ll) + mean(dmvvec) ) + 2*arma::as_scalar(dof)); 
  
}

