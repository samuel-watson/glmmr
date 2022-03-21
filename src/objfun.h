
#ifndef OBJFUN_H_
#define OBJFUN_H_

#include <RcppArmadillo.h>

namespace rminqa {
class ObjFun {
public:
  ObjFun() {}
  virtual ~ObjFun() {}
  virtual double operator()(const arma::vec &par) = 0;
};

inline double minqa_objfun(long n, const double *x, void *data) {
  arma::vec par(x, n);
  return static_cast<ObjFun *>(data)->operator()(par);
}

}

#endif
