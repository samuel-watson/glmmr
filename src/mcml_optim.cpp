#include <cmath>  // std::pow
#include <RcppArmadillo.h>
// #include <roptim.h>
#include "rbobyqa.h"
#include "glmmrmath.h"
#include "glmmrmatrix.h"
using namespace rminqa;
using namespace Rcpp;
using namespace arma;
// using namespace roptim;

// [[Rcpp::depends(RcppArmadillo)]]
// not include // [[Rcpp::depends(roptim)]]

class D_likelihood : public ObjFun {
  Rcpp::List func_;
  Rcpp::List data_;
  arma::mat u_;
  
public:
  D_likelihood(Rcpp::List func, Rcpp::List data, arma::mat u) : 
    func_(func), data_(data), u_(u) {}
  double operator()(const vec &par) override{
    arma::uword nrow = u_.n_rows;
    // arma::uword Q = u_.n_cols;
    
    arma::mat D = genD(func_,
                       data_,
                       par);
    double logdetD = arma::log_det_sympd(D);
    D = arma::inv_sympd(D);
    double dmv = 0;
    
    for(arma::uword j=0;j<nrow;j++){
      dmv += log_mv_gaussian_pdf(u_.row(j),D,logdetD);
    }
    // Rcpp::Rcout << "\n" << dmv;
    //double check that the matrix u is orientated
    // as draws * random effects (cols)
    return -1 * dmv/nrow;
  }
};


// [[Rcpp::export]]
arma::vec d_lik_optim(const Rcpp::List &func,
                      const Rcpp::List &data,
                      const arma::mat &u,
                      arma::vec start,
                      const arma::vec &lower,
                      const arma::vec &upper,
                      int trace = 0){
  D_likelihood dl(func, data, u);
  
  //Roptim<D_likelihood> opt("Nelder-Mead");
  Rbobyqa<D_likelihood> opt;
  opt.set_upper(upper);
  opt.set_lower(lower);
  opt.control.iprint = trace;
  //opt.set_hessian(true);
  opt.minimize(dl, start);
  //opt.par().print("par = ");
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
  L_likelihood(arma::mat Z, arma::mat X,
               arma::vec y, arma::mat u, 
               std::string family, std::string link) : 
  Z_(Z), X_(X), y_(y), u_(u),family_(family) , link_(link) {}
  double operator()(const vec &par) override{
    arma::uword niter = u_.n_rows;
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
      zd = Z_ * u_.row(j).t();
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
  Rcpp::List func_;
  Rcpp::List data_;
  arma::mat Z_;
  arma::mat X_;
  arma::vec y_;
  arma::mat u_;
  arma::vec cov_par_fix_;
  std::string family_;
  std::string link_;
  bool importance_;
  
public:
  F_likelihood(Rcpp::List func,
               Rcpp::List data,
               arma::mat Z, 
               arma::mat X,
               arma::vec y, 
               arma::mat u,
               arma::vec cov_par_fix,
               std::string family, 
               std::string link,
               bool importance) : 
  func_(func), data_(data), Z_(Z), X_(X), y_(y), 
  u_(u),cov_par_fix_(cov_par_fix), family_(family) , link_(link),
  importance_(importance){}
  double operator()(const vec &par) override{
    arma::uword niter = u_.n_rows;
    arma::uword n = y_.n_elem;
    arma::uword P = X_.n_cols;
    arma::uword Q = cov_par_fix_.n_elem;
    double du;
    
    //now loop through all the iterations
    // log likelihood of D
    arma::mat Dnew = genD(func_,
                          data_,
                          par.subvec(P,P+Q-1));
    double logdetDnew = arma::log_det_sympd(Dnew);
    Dnew = arma::inv_sympd(Dnew);
    arma::vec numerD(niter);
    //Rcpp::Rcout << logdetDnew  << std::endl;
    for(arma::uword j=0;j<niter;j++){
      double lldn = -0.5*(Q*log(2*arma::datum::pi) +
                          logdetDnew + arma::as_scalar(u_.row(j)*Dnew*u_.row(j).t()));
      if(importance_){
        numerD(j) = exp(lldn);
      } else {
        numerD(j) = lldn;
      }
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
    arma::vec lfa(niter);
    
    for(arma::uword j=0; j<niter ; j++){
      zd = Z_ * u_.row(j).t();
      double ll = log_likelihood(y_,
                                 xb + zd,
                                 var_par,
                                 family_,
                                 link_);
      if(importance_){
        lfa(j) = exp(ll);
      } else {
        lfa(j) = ll;
      }
      
    }
    
    if(importance_){
      // denominator density for importance sampling
      arma::mat D = genD(func_,
                         data_,
                         cov_par_fix_);
      
      double logdetD = arma::log_det_sympd(D);
      D = arma::inv_sympd(D);
      arma::vec denom(niter);
      du = 0;
      for(arma::uword j=0;j<niter;j++){
        double lld = -0.5*Q*log(2*arma::datum::pi)-
          0.5*logdetD - 0.5*arma::as_scalar(u_.row(j)*D*u_.row(j).t());
        du  += lfa(j)*numerD(j)/exp(lld);
      }
      
      // for(arma::uword j=0; j<niter ; j++){
      //   du += lfa(j)*numerD(j)/denom(j);
      // }
      du = -1 * log(du/niter);
    } else {
      du = -1* (mean(numerD) + mean(lfa));
    }
    //Rcpp::Rcout << du << numerD << lfa << std::endl;
    return du;
  }
};

// [[Rcpp::export]]
arma::vec f_lik_grad(Rcpp::List func,
                     Rcpp::List data,
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
  
  F_likelihood dl(func,data,Z,X,y,u,
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
arma::mat f_lik_hess(Rcpp::List func,
                     Rcpp::List data,
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
  F_likelihood dl(func,data,Z,X,y,u,
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
arma::mat f_lik_optim(Rcpp::List func,
                      Rcpp::List data,
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
  
  F_likelihood dl(func,data,Z,X,y,u,
                  cov_par_fix,family,
                  link,true);
  
  Rbobyqa<F_likelihood> opt;
  opt.set_lower(lower);
  opt.control.iprint = trace;
  opt.minimize(dl, start);
  
  return opt.par();
  
  // Roptim<F_likelihood> opt;
  // opt.control.trace = 0;
  // opt.set_lower(lower);
  // if(!importance){
  //   opt.set_hessian(true);
  // }
  // opt.minimize(dl, start);
  // 
  // if(importance){
  //   return opt.par();
  // } else {
  //   return opt.hessian();
  // }
  
}

// // [[Rcpp::export]]
// arma::mat f_lik_optim2(Rcpp::List func,
//                       Rcpp::List data,
//                       const arma::mat &Z,
//                       const arma::mat &X,
//                       const arma::vec &y,
//                       const arma::mat &u,
//                       const arma::vec &cov_par_fix,
//                       std::string family,
//                       std::string link,
//                       arma::vec start,
//                       const arma::vec &lower,
//                       int trace){
// 
//   F_likelihood dl(func,data,Z,X,y,u,
//                   cov_par_fix,family,
//                   link,false);
// 
// 
//   Roptim<F_likelihood> opt("L-BFGS-B");
//   opt.control.trace = 0;
//   opt.set_lower(lower);
//   opt.set_hessian(true);
//   opt.minimize(dl, start);
// 
//   return opt.hessian();
// 
// }

// [[Rcpp::export]]
Rcpp::List mcnr_step(const arma::vec &y, const arma::mat &X, const arma::mat &Z,
                     const arma::vec &beta, const arma::mat &u,
                     const std::string &family, const std::string &link){
  const arma::uword n = y.n_elem;
  const arma::uword P = X.n_cols;
  const arma::uword niter = u.n_rows;
  
  //generate residuals
  arma::vec xb = X*beta;
  arma::vec sigmas(niter);
  arma::mat XtWX(P,P,fill::zeros);
  arma::vec Wu(n,fill::zeros);
  for(arma::uword i = 0; i < niter; ++i){
    arma::vec zd = Z * u.row(i).t();
    arma::vec mu = mod_inv_func(xb + zd, link);
    arma::vec resid = y - mu;
    sigmas(i) = arma::stddev(resid);
    arma::vec wdiag = gen_dhdmu(xb + zd,family,link);
    arma::vec wdiag2 = 1/arma::pow(wdiag, 2);
    if(family=="gaussian" || family=="gamma") wdiag2 *= sigmas(i);
    arma::mat W = diagmat(wdiag2);
    arma::mat W2 = diagmat(wdiag);
    XtWX += X.t()*W*X;
    Wu += W*W2*resid;
  }
  //Rcpp::Rcout<< XtWX;
  XtWX = arma::inv_sympd(XtWX/niter);
  Wu = Wu/niter;
  
  arma::vec bincr = XtWX*X.t()*Wu;
  Rcpp::List L = List::create(_["beta_step"] = bincr , _["sigmahat"] = arma::mean(sigmas));
  return L;
}

// // [[Rcpp::export]]
// Rcpp::List mcnr_step(arma::vec &y,
//                      arma::mat &X,
//                      arma::mat &Z,
//                      arma::vec &beta,
//                      arma::mat &u,
//                      std::string family,
//                      std::string link){
//   arma::uword n = y.n_elem;
//   arma::uword P = X.n_cols;
//   arma::uword niter = u.n_rows;
//   arma::mat XtWX(P,P,fill::zeros);
//   arma::vec Wu(n,fill::zeros);
//   
//   //generate residuals
//   arma::vec xb = X*beta;
//   arma::vec resid(n);
//   arma::vec mu(n);
//   arma::vec zd(n);
//   arma::vec sigmas(niter);
//   arma::vec wdiag(xb.n_elem,fill::zeros);
//   arma::vec wdiag2(xb.n_elem,fill::zeros);
//   //arma::sp_mat W(n,n);
//   
//   for(arma::uword i=0;i<niter;i++){
//     zd = Z * u.row(i).t();
//     mu = mod_inv_func(xb + zd,
//                       link);
//     resid = y - mu;
//     sigmas(i) = arma::stddev(resid);
//     wdiag = gen_dhdmu(xb + zd,family,link);
//     for(arma::uword j=0;j<wdiag2.n_elem;j++){
//       wdiag2(j) = 1/(pow(wdiag(j),2));
//       if(family=="gaussian" || family=="gamma"){
//         wdiag2(j) = sigmas(i)*wdiag2(j);
//       }
//     }
//     arma::mat W = diagmat(wdiag2);
//     arma::mat W2 = diagmat(wdiag);
//     XtWX += X.t()*W*X;
//     Wu += W*W2*resid;
//   }
//   //Rcpp::Rcout<< XtWX;
//   XtWX = arma::inv_sympd(XtWX/niter);
//   Wu = Wu/niter;
//   
//   arma::vec bincr = XtWX*X.t()*Wu;
//   Rcpp::List L = List::create(_["beta_step"] = bincr , _["sigmahat"] = arma::mean(sigmas));
//   return L;
// }

// [[Rcpp::export]]
double aic_mcml(const arma::mat &Z, 
                const arma::mat &X,
                const arma::vec &y, 
                const arma::mat &u, 
                std::string family, 
                std::string link,
                Rcpp::List func,
                Rcpp::List data,
                const arma::vec& beta_par,
                const arma::vec& cov_par){
  arma::uword niter = u.n_rows;
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
  
  arma::mat D = genD(func,
                     data,
                     cov_par);
  double logdetD = arma::log_det_sympd(D);
  D = arma::inv_sympd(D);
  double dmv = 0;
  
  double lfa;
  
  for(arma::uword j=0; j<niter ; j++){
    zd = Z * u.row(j).t();
    double ll = log_likelihood(y,
                               xb + zd,
                               var_par,
                               family,
                               link);
    
    dmv += log_mv_gaussian_pdf(u.row(j),D,logdetD);
    lfa += ll;
  }
  
  return (-2*( lfa/niter + dmv/niter ) + 2*arma::as_scalar(dof)); 
  
}

