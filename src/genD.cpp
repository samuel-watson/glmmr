#include <cmath>  // std::pow
#include <RcppArmadillo.h>
#include <roptim.h>
#include "rbobyqa.h"
using namespace rminqa;
using namespace Rcpp;
using namespace arma;
using namespace roptim;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(roptim)]]

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
    out = pars(funcdetail(5))*exp(-1*val*pars(funcdetail(6)));
  }
  
  //power exponential
  if(f1==3){
    out = pow(pars(funcdetail(5)),val);
  }
  
  return(out);
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
double log_factorial_approx(int n){
  double ans;
  if(n==0){
    ans = 0;
  } else {
    ans = n*log(n) - n + log(n*(1+4*n*(1+2*n)))/6 + log(arma::datum::pi)/2;
  }
  return ans;
}

class D_likelihood : public ObjFun {
  Rcpp::List func_;
  Rcpp::List data_;
  arma::mat u_;
  
public:
  D_likelihood(Rcpp::List func, Rcpp::List data, arma::mat u) : func_(func), data_(data), u_(u) {}
  double operator()(const vec &par) override{
    arma::uword nrow = u_.n_rows;
    arma::uword Q = u_.n_cols;
    
    arma::mat D = genD(func_,
                       data_,
                       par);
    double logdetD = arma::log_det_sympd(D);
    D = arma::inv_sympd(D);
    double dmv = 0;
    
    for(arma::uword j=0;j<nrow;j++){
      dmv += -0.5*Q*log(2*arma::datum::pi)-
        0.5*logdetD - 0.5*arma::as_scalar(u_.row(j)*D*u_.row(j).t());
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


class F_likelihood : public Functor {
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
    //Rcpp::Rcout << par  << std::endl;
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
arma::mat f_lik_optim(Rcpp::List func,
                      Rcpp::List data,
                      arma::mat Z, 
                      arma::mat X,
                      arma::vec y, 
                      arma::mat u,
                      arma::vec cov_par_fix,
                      std::string family, 
                      std::string link,
                      arma::vec start,
                      arma::vec lower,
                      bool importance){
  
  F_likelihood dl(func,data,Z,X,y,u,
                  cov_par_fix,family,
                  link,importance);
  
  Roptim<F_likelihood> opt("L-BFGS-B");
  opt.control.trace = 0;
  opt.set_lower(lower);
  if(!importance){
    opt.set_hessian(true);
  }
  opt.minimize(dl, start);
  
  if(importance){
    return opt.par();
  } else {
    return opt.hessian();
  }
  
}

// [[Rcpp::export]]
arma::vec inv_func(arma::vec mu,
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
  
  return mu;
  
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
  
  return mu;
  
}

// [[Rcpp::export]]
Rcpp::List mcnr_step(arma::vec &y,
                    arma::mat &X,
                    arma::mat &Z,
                    arma::vec &beta,
                    arma::mat &u,
                    std::string link){
  arma::uword n = y.n_elem;
  arma::uword P = X.n_cols;
  arma::uword niter = u.n_rows;
  arma::mat XtWX(P,P,fill::zeros);
  arma::vec Wu(n,fill::zeros);
  
  //generate residuals
  arma::vec xb = X*beta;
  arma::vec resid(n);
  arma::vec mu(n);
  arma::vec zd(n);
  arma::vec sigmas(niter);
  //arma::sp_mat W(n,n);
  
  for(arma::uword i=0;i<niter;i++){
    zd = Z * u.row(i).t();
    mu = mod_inv_func(xb + zd,
                  link);
    resid = y - mu;
    sigmas(i) = arma::stddev(resid);
    arma::mat W = diagmat(arma::vec(n,fill::value(arma::as_scalar(sigmas(i)))));
    XtWX += X.t()*W*X;
    Wu += W*resid;
  }
  XtWX = arma::inv_sympd(XtWX/niter);
  Wu = Wu/niter;
  
  arma::vec bincr = XtWX*X.t()*Wu;
  Rcpp::List L = List::create(_["beta_step"] = bincr , _["sigmahat"] = arma::mean(sigmas));
  return L;
}

