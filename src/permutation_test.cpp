#include <string>

#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "glm.h"
#include "glmmr.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;

typedef MatrixXd::Index Index;

// [[Rcpp::export]]
List fast_glm_impl(Rcpp::NumericMatrix Xs,
             Rcpp::NumericVector ys,
             Rcpp::NumericVector weightss,
             Rcpp::NumericVector offsets,
             Rcpp::NumericVector starts,
             Rcpp::NumericVector mus,
             Rcpp::NumericVector etas,
             Function var,
             Function mu_eta,
             Function linkinv,
             Function dev_resids,
             Function valideta,
             Function validmu,
             int type,
             double tol,
             int maxit)
{
  const Map<MatrixXd>  X(as<Map<MatrixXd> >(Xs));
  const Map<VectorXd>  y(as<Map<VectorXd> >(ys));
  const Map<VectorXd>  weights(as<Map<VectorXd> >(weightss));
  const Map<VectorXd>  offset(as<Map<VectorXd> >(offsets));
  const Map<VectorXd>  beta_init(as<Map<VectorXd> >(starts));
  const Map<VectorXd>  mu_init(as<Map<VectorXd> >(mus));
  const Map<VectorXd>  eta_init(as<Map<VectorXd> >(etas));
  Index                n = X.rows();
  if ((Index)y.size() != n) throw invalid_argument("size mismatch");

  // instantiate fitting class
  GlmBase<Eigen::VectorXd, Eigen::MatrixXd> *glm_solver = NULL;

  bool is_big_matrix = false;

  glm_solver = new glm(X, y, weights, offset,
                       var, mu_eta, linkinv, dev_resids,
                       valideta, validmu, tol, maxit, type,
                       is_big_matrix);

  // initialize parameters
  glm_solver->init_parms(beta_init, mu_init, eta_init);


  // maximize likelihood
  int iters = glm_solver->solve(maxit);

  VectorXd beta      = glm_solver->get_beta();
  VectorXd se        = glm_solver->get_se();
  VectorXd mu        = glm_solver->get_mu();
  VectorXd eta       = glm_solver->get_eta();
  VectorXd wts       = glm_solver->get_w();
  VectorXd pweights  = glm_solver->get_weights();

  double dev         = glm_solver->get_dev();
  int rank           = glm_solver->get_rank();
  bool converged     = glm_solver->get_converged();

  int df = X.rows() - rank;

  delete glm_solver;

  return List::create(_["coefficients"]      = beta,
                      _["se"]                = se,
                      _["fitted.values"]     = mu,
                      _["linear.predictors"] = eta,
                      _["deviance"]          = dev,
                      _["weights"]           = wts,
                      _["prior.weights"]     = pweights,
                      _["rank"]              = rank,
                      _["df.residual"]       = df,
                      _["iter"]              = iters,
                      _["converged"]         = converged);
}

//' Simplified version of fastglm's `fastglm` function
//' 
//' Fast generalized model fitting with a simplified version of fastglm's `fastglm` function
//' @details
//' This is a simplified wrapper to the `fastglm` function in the fastglm package, please access that package for more details.
//' @param Xs input model matrix. Must be a matrix object
//' @param y numeric response vector
//' @param weightss an optional vector of 'prior weights' to be used in the fitting process. Should be a numeric vector.
//' @param offsets this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be a numeric vector of length equal to the number of cases
//' @param A \link{stats}[family] object. See `fastglm` for details
//' @return A list with the elements:
//' * `coefficients` a vector of coefficients
//' * `se` a vector of the standard errors of the coefficient estimates
//' * `rank` a scalar denoting the computed rank of the model matrix
//' * `df.residual` a scalar denoting the degrees of freedom in the model
//' * `residuals` a vector of residuals
//' * `s` a numeric scalar - the root mean square for residuals
//' * `fitted.values` the fitted values
// [[Rcpp::export]]
Rcpp::List myglm(Rcpp::NumericMatrix Xs,
         Rcpp::NumericVector ys,
         Rcpp::NumericVector weightss,
         Rcpp::NumericVector offsets,
	 Rcpp::List family) {
    Rcpp::NumericVector mustart = ys;
    Rcpp::NumericVector startval(Xs.ncol());
    // Rcpp::Rcout << "startval = " << startval << std::endl;

    Function linkfun = family["linkfun"];
    Function linkinv = family["linkinv"];
    Rcpp::NumericVector eta = linkfun(mustart);
    Rcpp::NumericVector mu = linkinv(eta);

    return fast_glm_impl(Xs, ys, weightss, offsets,
                                   startval, mu, eta,
                                   family["variance"], family["mu.eta"],
                                   family["linkinv"], family["dev.resids"],
                                   family["valideta"], family["validmu"],
                                   0, 1e-7, 100);
}

arma::vec get_G(const arma::vec &x, 
                       const Rcpp::String &family2){
  arma::vec dx;
  if(family2 == "identity"){
    dx = arma::ones<arma::vec>(x.n_elem);
  } else if(family2 == "log"){
    dx = exp(x);
  } else if(family2 == "logit"){
    dx = (exp(x)/(1+exp(x))) % (1-exp(x)/(1+exp(x)));
  } else if(family2 == "probit"){
    dx = -1/gaussian_pdf_vec(x);
  }
  return dx;
}

//' The quasi-score statistic for a generalised linear mixed model
//' 
//' Generates the quasi-score statistic for a generalised linear mixed model
//' 
//' @param resids A numeric vector of generalised residuals
//' @param tr A numeric vector of 1s (treatment group) and -1s (control group)
//' @param xb A numeric vector of fitted linear predictors
//' @param invS A matrix. If using the weighted statistic then it should be the inverse covariance matrix of the observations
//' @param family2 A string naming the link function
//' @param weight Logical value indicating whether to use the weighted statistic (TRUE) or the unweighted statistic (FALSE)
//' @return A scalar value with the value of the statistic
// [[Rcpp::export]]
double qscore_impl(const arma::vec &resids,
                   arma::vec tr,
                   const arma::vec &xb,
                   const arma::mat &invS,
                   const std::string &family2,
                   bool weight=true) {
  arma::mat g;
  double q;
  if (weight){
    //Rcpp::String str = family2(0);
    g = arma::trans(get_G(xb, family2));
    q = arma::as_scalar((g * invS) * (tr % resids));
  } else {
    q = arma::as_scalar(arma::dot(tr, resids));
  }
  return std::abs(q);
}

//' Generates realisations of the permutational test statistic distribution 
//' 
//' Generates realisations of the permutational test statistic distribution from a given matrix of permutations
//' 
//' @param resids A numeric vector of generalised residuals
//' @param tr_mat A matrix. Each column is a new random treatment allocation with 1s (treatment group) and 0s (control group)
//' @param xb A numeric vector of fitted linear predictors
//' @param invS A matrix. If using the weighted statistic then it should be the inverse covariance matrix of the observations
//' @param family2 A string naming the link function
//' @param weight Logical value indicating whether to use the weighted statistic (TRUE) or the unweighted statistic (FALSE)
//' @param verbose Logical indicating whether to report detailed output
//' @return A numeric vector of quasi-score test statistics for each of the permutations
// [[Rcpp::export]]
arma::vec permutation_test_impl(const arma::vec &resids,
                                const arma::mat &tr_mat,
                                const arma::vec &xb,
                                const arma::mat &invS,
                                const std::string &family2,
                                bool weight,
                                int iter = 1000,
                                bool verbose = true) {
  if (verbose) Rcpp::Rcout << "Starting permutations\n" << std::endl;

  arma::vec qtest = arma::zeros<arma::vec>(iter);
#pragma omp parallel for
  for (int i = 0; i < iter; ++i) {
      arma::vec tr = tr_mat.col(i);
      tr.replace(0, -1);
      qtest[i] = qscore_impl(resids, tr, xb, invS, family2, weight);
  }
  return qtest;
}

//' Confidence interval search procedure
//' 
//' Search for the bound of a confidence interval using permutation test statistics
//' 
//' @param start Numeric value indicating the starting value for the search procedure
//' @param b Numeric value indicating the parameter estimate
//' @param Xnull_ Numeric matrix. The covariate design matrix with the treatment variable removed
//' @param y_ Numeric vector of response variables
//' @param tr_ Numeric vector. The original random allocation (0s and 1s)
//' @param new_tr_mat A matrix. Each column is a new random treatment allocation with 1s (treatment group) and 0s (control group)
//' @param xb A numeric vector of fitted linear predictors
//' @param invS A matrix. If using the weighted statistic then it should be the inverse covariance matrix of the observations
//' @param family A \link{stats}[family] object
//' @param family2 A string naming the link function
//' @param nsteps Integer specifying the number of steps of the search procedure
//' @param weight Logical indicating whether to use the weighted (TRUE) or unweighted (FALSE) test statistic
//' @param alpha The function generates (1-alpha)*100% confidence intervals. Default is 0.05.
//' @param verbose Logical indicating whether to provide detailed output.
//' @return The estimated confidence interval bound
// [[Rcpp::export]]
double confint_search(double start,
                      double b,
                      Rcpp::NumericMatrix Xnull_,
                      Rcpp::NumericVector y_,
                      Rcpp::NumericVector tr_,
                      const arma::mat &new_tr_mat,
                      const arma::vec &xb,
                      const arma::mat &invS,
                      Rcpp::List family,
                      const std::string &family2,
                      int nsteps = 1000,
                      bool weight = true,
                      double alpha = 0.05,
                      bool verbose = true) {
  arma::mat Xull = Rcpp::as<arma::mat>(Xnull_);
  arma::vec y = Rcpp::as<arma::vec>(y_);
  arma::vec tr = Rcpp::as<arma::vec>(tr_);
  tr.replace(0,-1);

  Rcpp::NumericVector weights_(y.n_elem);
  weights_.fill(1);
  Function linkinv = family["linkinv"];

  double bound = start;
  double k_tmp = gaussian_cdf(1-alpha);
  double k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));

  double qstat, qtest;
  for (int i = 1; i <= nsteps; ++i) {
    Rcpp::List null_fit = myglm(Xnull_, y_, weights_, bound*tr_, family);
    Rcpp::NumericVector xb = null_fit["linear.predictors"];
    arma::vec ypred = Rcpp::as<arma::vec>(linkinv(xb));
    arma::vec resids = y - ypred;
    qstat = qscore_impl(resids, tr, xb, invS, family2, weight);

    arma::vec new_tr = new_tr_mat.col(i-1);
    new_tr.replace(0, -1);
    qtest = qscore_impl(resids, new_tr, xb, invS, family2, weight);
    
    bool rjct = qstat > qtest;
    double step = k * (b - bound);
    if (rjct) {
      bound += step*alpha/i;
    } else {
      bound -= step*(1-alpha)/i;
    }
    
    if(verbose){
      Rcpp::Rcout << "\r Step = " << i << " bound: " << bound << std::endl;
    }
  }
  return bound;
}
