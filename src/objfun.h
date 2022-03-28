
#ifndef OBJFUN_H_
#define OBJFUN_H_

#include <RcppArmadillo.h>

namespace rminqa {

class ObjFun {
public:
  struct OptStruct {
    bool has_grad_ = false;
    bool has_hess_ = false;
    arma::vec ndeps_;       // tolerances for numerical derivatives
    double fnscale_ = 1.0;  // scaling for objective
    arma::vec parscale_;    // scaling for parameters
    int usebounds_ = 0;
    arma::vec lower_, upper_;
    bool sann_use_custom_function_ = false;
  } os;
public:
  ObjFun() {}
  virtual ~ObjFun() {}
  virtual double operator()(const arma::vec &par) = 0;
  virtual void Gradient(const arma::vec &par, arma::vec &grad) {
    ApproximateGradient(par, grad);
  }
  virtual void Hessian(const arma::vec &par, arma::mat &hess) {
    ApproximateHessian(par, hess);
  }
  
  // Returns forward-difference approximation of Gradient
  void ApproximateGradient(const arma::vec &par, arma::vec &grad);
  
  // Returns forward-difference approximation of Hessian
  void ApproximateHessian(const arma::vec &par, arma::mat &hess);
};

inline void ObjFun::ApproximateGradient(const arma::vec &par,
                                        arma::vec &grad) {
  if (os.parscale_.is_empty()) os.parscale_ = arma::ones<arma::vec>(par.size());
  if (os.ndeps_.is_empty())
    os.ndeps_ = arma::ones<arma::vec>(par.size()) * 1e-6;
  
  grad = arma::zeros<arma::vec>(par.size());
  arma::vec x = par % os.parscale_;
  
  if (os.usebounds_ == 0) {
    for (std::size_t i = 0; i != par.size(); ++i) {
      double eps = os.ndeps_(i);
      
      x(i) = (par(i) + eps) * os.parscale_(i);
      double val1 = operator()(x) / os.fnscale_;
      
      x(i) = (par(i) - eps) * os.parscale_(i);
      double val2 = operator()(x) / os.fnscale_;
      
      grad(i) = (val1 - val2) / (2 * eps);
      
      x(i) = par(i) * os.parscale_(i);
    }
  } else {  // use bounds
    for (std::size_t i = 0; i != par.size(); ++i) {
      double epsused = os.ndeps_(i);
      double eps = os.ndeps_(i);
      
      double tmp = par(i) + eps;
      if (tmp > os.upper_(i)) {
        tmp = os.upper_(i);
        epsused = tmp - par(i);
      }
      
      x(i) = tmp * os.parscale_(i);
      double val1 = operator()(x) / os.fnscale_;
      
      tmp = par(i) - eps;
      if (tmp < os.lower_(i)) {
        tmp = os.lower_(i);
        eps = par(i) - tmp;
      }
      
      x(i) = tmp * os.parscale_(i);
      double val2 = operator()(x) / os.fnscale_;
      
      grad(i) = (val1 - val2) / (epsused + eps);
      
      x(i) = par(i) * os.parscale_(i);
    }
  }
}

inline void ObjFun::ApproximateHessian(const arma::vec &par, arma::mat &hess) {
  if (os.parscale_.is_empty()) os.parscale_ = arma::ones<arma::vec>(par.size());
  if (os.ndeps_.is_empty())
    os.ndeps_ = arma::ones<arma::vec>(par.size()) * 1e-3;
  
  hess = arma::zeros<arma::mat>(par.size(), par.size());
  arma::vec dpar = par / os.parscale_;
  arma::vec df1 = arma::zeros<arma::vec>(par.size());
  arma::vec df2 = arma::zeros<arma::vec>(par.size());
  
  for (std::size_t i = 0; i != par.size(); ++i) {
    double eps = os.ndeps_(i) / os.parscale_(i);
    dpar(i) += eps;
    Gradient(dpar, df1);
    dpar(i) -= 2 * eps;
    Gradient(dpar, df2);
    for (std::size_t j = 0; j != par.size(); ++j)
      hess(i, j) = os.fnscale_ * (df1(j) - df2(j)) /
        (2 * eps * os.parscale_(i) * os.parscale_(j));
    dpar(i) = dpar(i) + eps;
  }
  
  // now symmetrize
  for (std::size_t i = 0; i != par.size(); ++i) {
    for (std::size_t j = 0; j != par.size(); ++j) {
      double tmp = 0.5 * (hess(i, j) + hess(j, i));
      
      hess(i, j) = tmp;
      hess(j, i) = tmp;
    }
  }
}

inline double minqa_objfun(long n, const double *x, void *data) {
  arma::vec par(x, n);
  return static_cast<ObjFun *>(data)->operator()(par);
}

} // namespace rminqa

#endif
