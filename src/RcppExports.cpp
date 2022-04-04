// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dfbeta_stat
arma::rowvec dfbeta_stat(const arma::mat& sigma, const arma::mat& X, const arma::vec& y, arma::uword par);
RcppExport SEXP _glmmr_dfbeta_stat(SEXP sigmaSEXP, SEXP XSEXP, SEXP ySEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uword >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(dfbeta_stat(sigma, X, y, par));
    return rcpp_result_gen;
END_RCPP
}
// GradRobustStep
Rcpp::List GradRobustStep(arma::uword N, arma::uvec idx_in, Rcpp::List C_list, Rcpp::List X_list, Rcpp::List sig_list, arma::vec weights, arma::uvec exp_cond, arma::uvec nfix, arma::uword any_fix, arma::uword rd_mode, bool trace);
RcppExport SEXP _glmmr_GradRobustStep(SEXP NSEXP, SEXP idx_inSEXP, SEXP C_listSEXP, SEXP X_listSEXP, SEXP sig_listSEXP, SEXP weightsSEXP, SEXP exp_condSEXP, SEXP nfixSEXP, SEXP any_fixSEXP, SEXP rd_modeSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type idx_in(idx_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type C_list(C_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sig_list(sig_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type exp_cond(exp_condSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type nfix(nfixSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type any_fix(any_fixSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type rd_mode(rd_modeSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(GradRobustStep(N, idx_in, C_list, X_list, sig_list, weights, exp_cond, nfix, any_fix, rd_mode, trace));
    return rcpp_result_gen;
END_RCPP
}
// log_factorial_approx
double log_factorial_approx(int n);
RcppExport SEXP _glmmr_log_factorial_approx(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(log_factorial_approx(n));
    return rcpp_result_gen;
END_RCPP
}
// log_mv_gaussian_pdf
double log_mv_gaussian_pdf(const arma::vec& u, const arma::mat& D, const double& logdetD);
RcppExport SEXP _glmmr_log_mv_gaussian_pdf(SEXP uSEXP, SEXP DSEXP, SEXP logdetDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const double& >::type logdetD(logdetDSEXP);
    rcpp_result_gen = Rcpp::wrap(log_mv_gaussian_pdf(u, D, logdetD));
    return rcpp_result_gen;
END_RCPP
}
// gen_dhdmu
arma::vec gen_dhdmu(arma::vec xb, std::string family, std::string link);
RcppExport SEXP _glmmr_gen_dhdmu(SEXP xbSEXP, SEXP familySEXP, SEXP linkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_dhdmu(xb, family, link));
    return rcpp_result_gen;
END_RCPP
}
// fexp
double fexp(const double& x, double par1, double par2);
RcppExport SEXP _glmmr_fexp(SEXP xSEXP, SEXP par1SEXP, SEXP par2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type par1(par1SEXP);
    Rcpp::traits::input_parameter< double >::type par2(par2SEXP);
    rcpp_result_gen = Rcpp::wrap(fexp(x, par1, par2));
    return rcpp_result_gen;
END_RCPP
}
// sqexp
double sqexp(const double& x, double par1, double par2);
RcppExport SEXP _glmmr_sqexp(SEXP xSEXP, SEXP par1SEXP, SEXP par2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type par1(par1SEXP);
    Rcpp::traits::input_parameter< double >::type par2(par2SEXP);
    rcpp_result_gen = Rcpp::wrap(sqexp(x, par1, par2));
    return rcpp_result_gen;
END_RCPP
}
// matern
double matern(const double& x, double rho, double nu);
RcppExport SEXP _glmmr_matern(SEXP xSEXP, SEXP rhoSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(matern(x, rho, nu));
    return rcpp_result_gen;
END_RCPP
}
// bessel1
double bessel1(const double& x, double rho);
RcppExport SEXP _glmmr_bessel1(SEXP xSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(bessel1(x, rho));
    return rcpp_result_gen;
END_RCPP
}
// blockMat
arma::mat blockMat(arma::field<arma::mat> matfield);
RcppExport SEXP _glmmr_blockMat(SEXP matfieldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type matfield(matfieldSEXP);
    rcpp_result_gen = Rcpp::wrap(blockMat(matfield));
    return rcpp_result_gen;
END_RCPP
}
// genBlockD
arma::mat genBlockD(size_t N_dim, size_t N_func, const arma::uvec& func_def, const arma::uvec& N_var_func, const arma::umat& col_id, const arma::uvec& N_par, const arma::mat& cov_data, const arma::vec& gamma);
RcppExport SEXP _glmmr_genBlockD(SEXP N_dimSEXP, SEXP N_funcSEXP, SEXP func_defSEXP, SEXP N_var_funcSEXP, SEXP col_idSEXP, SEXP N_parSEXP, SEXP cov_dataSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< size_t >::type N_dim(N_dimSEXP);
    Rcpp::traits::input_parameter< size_t >::type N_func(N_funcSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type func_def(func_defSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_var_func(N_var_funcSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type col_id(col_idSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_par(N_parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cov_data(cov_dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(genBlockD(N_dim, N_func, func_def, N_var_func, col_id, N_par, cov_data, gamma));
    return rcpp_result_gen;
END_RCPP
}
// genD
arma::field<arma::mat> genD(const arma::uword& B, const arma::uvec& N_dim, const arma::uvec& N_func, const arma::umat& func_def, const arma::umat& N_var_func, const arma::ucube& col_id, const arma::umat& N_par, const arma::uword& sum_N_par, const arma::cube& cov_data, const arma::vec& gamma);
RcppExport SEXP _glmmr_genD(SEXP BSEXP, SEXP N_dimSEXP, SEXP N_funcSEXP, SEXP func_defSEXP, SEXP N_var_funcSEXP, SEXP col_idSEXP, SEXP N_parSEXP, SEXP sum_N_parSEXP, SEXP cov_dataSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_dim(N_dimSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_func(N_funcSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type func_def(func_defSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_var_func(N_var_funcSEXP);
    Rcpp::traits::input_parameter< const arma::ucube& >::type col_id(col_idSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_par(N_parSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type sum_N_par(sum_N_parSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cov_data(cov_dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(genD(B, N_dim, N_func, func_def, N_var_func, col_id, N_par, sum_N_par, cov_data, gamma));
    return rcpp_result_gen;
END_RCPP
}
// remove_one_many_mat
arma::mat remove_one_many_mat(const arma::mat& A, const arma::uvec& i);
RcppExport SEXP _glmmr_remove_one_many_mat(SEXP ASEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(remove_one_many_mat(A, i));
    return rcpp_result_gen;
END_RCPP
}
// add_one_mat
arma::mat add_one_mat(const arma::mat& A, double sigma_jj, const arma::vec& f);
RcppExport SEXP _glmmr_add_one_mat(SEXP ASEXP, SEXP sigma_jjSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type sigma_jj(sigma_jjSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(add_one_mat(A, sigma_jj, f));
    return rcpp_result_gen;
END_RCPP
}
// d_lik_optim
arma::vec d_lik_optim(const arma::uword& B, const arma::uvec& N_dim, const arma::uvec& N_func, const arma::umat& func_def, const arma::umat& N_var_func, const arma::ucube& col_id, const arma::umat& N_par, const arma::uword& sum_N_par, const arma::cube& cov_data, const arma::mat& u, arma::vec start, const arma::vec& lower, const arma::vec& upper, int trace);
RcppExport SEXP _glmmr_d_lik_optim(SEXP BSEXP, SEXP N_dimSEXP, SEXP N_funcSEXP, SEXP func_defSEXP, SEXP N_var_funcSEXP, SEXP col_idSEXP, SEXP N_parSEXP, SEXP sum_N_parSEXP, SEXP cov_dataSEXP, SEXP uSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_dim(N_dimSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_func(N_funcSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type func_def(func_defSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_var_func(N_var_funcSEXP);
    Rcpp::traits::input_parameter< const arma::ucube& >::type col_id(col_idSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_par(N_parSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type sum_N_par(sum_N_parSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cov_data(cov_dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< int >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(d_lik_optim(B, N_dim, N_func, func_def, N_var_func, col_id, N_par, sum_N_par, cov_data, u, start, lower, upper, trace));
    return rcpp_result_gen;
END_RCPP
}
// l_lik_optim
arma::vec l_lik_optim(const arma::mat& Z, const arma::mat& X, const arma::vec& y, const arma::mat& u, std::string family, std::string link, arma::vec start, const arma::vec& lower, const arma::vec& upper, int trace);
RcppExport SEXP _glmmr_l_lik_optim(SEXP ZSEXP, SEXP XSEXP, SEXP ySEXP, SEXP uSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< int >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lik_optim(Z, X, y, u, family, link, start, lower, upper, trace));
    return rcpp_result_gen;
END_RCPP
}
// f_lik_grad
arma::vec f_lik_grad(const arma::uword& B, const arma::uvec& N_dim, const arma::uvec& N_func, const arma::umat& func_def, const arma::umat& N_var_func, const arma::ucube& col_id, const arma::umat& N_par, const arma::uword& sum_N_par, const arma::cube& cov_data, const arma::mat& Z, const arma::mat& X, const arma::vec& y, const arma::mat& u, const arma::vec& cov_par_fix, std::string family, std::string link, arma::vec start, const arma::vec& lower, const arma::vec& upper, double tol);
RcppExport SEXP _glmmr_f_lik_grad(SEXP BSEXP, SEXP N_dimSEXP, SEXP N_funcSEXP, SEXP func_defSEXP, SEXP N_var_funcSEXP, SEXP col_idSEXP, SEXP N_parSEXP, SEXP sum_N_parSEXP, SEXP cov_dataSEXP, SEXP ZSEXP, SEXP XSEXP, SEXP ySEXP, SEXP uSEXP, SEXP cov_par_fixSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_dim(N_dimSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_func(N_funcSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type func_def(func_defSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_var_func(N_var_funcSEXP);
    Rcpp::traits::input_parameter< const arma::ucube& >::type col_id(col_idSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_par(N_parSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type sum_N_par(sum_N_parSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cov_data(cov_dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cov_par_fix(cov_par_fixSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(f_lik_grad(B, N_dim, N_func, func_def, N_var_func, col_id, N_par, sum_N_par, cov_data, Z, X, y, u, cov_par_fix, family, link, start, lower, upper, tol));
    return rcpp_result_gen;
END_RCPP
}
// f_lik_hess
arma::mat f_lik_hess(const arma::uword& B, const arma::uvec& N_dim, const arma::uvec& N_func, const arma::umat& func_def, const arma::umat& N_var_func, const arma::ucube& col_id, const arma::umat& N_par, const arma::uword& sum_N_par, const arma::cube& cov_data, const arma::mat& Z, const arma::mat& X, const arma::vec& y, const arma::mat& u, const arma::vec& cov_par_fix, std::string family, std::string link, arma::vec start, const arma::vec& lower, const arma::vec& upper, double tol);
RcppExport SEXP _glmmr_f_lik_hess(SEXP BSEXP, SEXP N_dimSEXP, SEXP N_funcSEXP, SEXP func_defSEXP, SEXP N_var_funcSEXP, SEXP col_idSEXP, SEXP N_parSEXP, SEXP sum_N_parSEXP, SEXP cov_dataSEXP, SEXP ZSEXP, SEXP XSEXP, SEXP ySEXP, SEXP uSEXP, SEXP cov_par_fixSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_dim(N_dimSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_func(N_funcSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type func_def(func_defSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_var_func(N_var_funcSEXP);
    Rcpp::traits::input_parameter< const arma::ucube& >::type col_id(col_idSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_par(N_parSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type sum_N_par(sum_N_parSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cov_data(cov_dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cov_par_fix(cov_par_fixSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(f_lik_hess(B, N_dim, N_func, func_def, N_var_func, col_id, N_par, sum_N_par, cov_data, Z, X, y, u, cov_par_fix, family, link, start, lower, upper, tol));
    return rcpp_result_gen;
END_RCPP
}
// f_lik_optim
arma::mat f_lik_optim(const arma::uword& B, const arma::uvec& N_dim, const arma::uvec& N_func, const arma::umat& func_def, const arma::umat& N_var_func, const arma::ucube& col_id, const arma::umat& N_par, const arma::uword& sum_N_par, const arma::cube& cov_data, const arma::mat& Z, const arma::mat& X, const arma::vec& y, const arma::mat& u, const arma::vec& cov_par_fix, std::string family, std::string link, arma::vec start, const arma::vec& lower, const arma::vec& upper, int trace);
RcppExport SEXP _glmmr_f_lik_optim(SEXP BSEXP, SEXP N_dimSEXP, SEXP N_funcSEXP, SEXP func_defSEXP, SEXP N_var_funcSEXP, SEXP col_idSEXP, SEXP N_parSEXP, SEXP sum_N_parSEXP, SEXP cov_dataSEXP, SEXP ZSEXP, SEXP XSEXP, SEXP ySEXP, SEXP uSEXP, SEXP cov_par_fixSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_dim(N_dimSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_func(N_funcSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type func_def(func_defSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_var_func(N_var_funcSEXP);
    Rcpp::traits::input_parameter< const arma::ucube& >::type col_id(col_idSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_par(N_parSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type sum_N_par(sum_N_parSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cov_data(cov_dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cov_par_fix(cov_par_fixSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< int >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(f_lik_optim(B, N_dim, N_func, func_def, N_var_func, col_id, N_par, sum_N_par, cov_data, Z, X, y, u, cov_par_fix, family, link, start, lower, upper, trace));
    return rcpp_result_gen;
END_RCPP
}
// mcnr_step
Rcpp::List mcnr_step(const arma::vec& y, const arma::mat& X, const arma::mat& Z, const arma::vec& beta, const arma::mat& u, const std::string& family, const std::string& link);
RcppExport SEXP _glmmr_mcnr_step(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP uSEXP, SEXP familySEXP, SEXP linkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const std::string& >::type link(linkSEXP);
    rcpp_result_gen = Rcpp::wrap(mcnr_step(y, X, Z, beta, u, family, link));
    return rcpp_result_gen;
END_RCPP
}
// aic_mcml
double aic_mcml(const arma::mat& Z, const arma::mat& X, const arma::vec& y, const arma::mat& u, std::string family, std::string link, const arma::uword& B, const arma::uvec& N_dim, const arma::uvec& N_func, const arma::umat& func_def, const arma::umat& N_var_func, const arma::ucube& col_id, const arma::umat& N_par, const arma::uword& sum_N_par, const arma::cube& cov_data, const arma::vec& beta_par, const arma::vec& cov_par);
RcppExport SEXP _glmmr_aic_mcml(SEXP ZSEXP, SEXP XSEXP, SEXP ySEXP, SEXP uSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP BSEXP, SEXP N_dimSEXP, SEXP N_funcSEXP, SEXP func_defSEXP, SEXP N_var_funcSEXP, SEXP col_idSEXP, SEXP N_parSEXP, SEXP sum_N_parSEXP, SEXP cov_dataSEXP, SEXP beta_parSEXP, SEXP cov_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_dim(N_dimSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type N_func(N_funcSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type func_def(func_defSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_var_func(N_var_funcSEXP);
    Rcpp::traits::input_parameter< const arma::ucube& >::type col_id(col_idSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type N_par(N_parSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type sum_N_par(sum_N_parSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cov_data(cov_dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_par(beta_parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cov_par(cov_parSEXP);
    rcpp_result_gen = Rcpp::wrap(aic_mcml(Z, X, y, u, family, link, B, N_dim, N_func, func_def, N_var_func, col_id, N_par, sum_N_par, cov_data, beta_par, cov_par));
    return rcpp_result_gen;
END_RCPP
}
// fast_glm_impl
List fast_glm_impl(Rcpp::NumericMatrix Xs, Rcpp::NumericVector ys, Rcpp::NumericVector weightss, Rcpp::NumericVector offsets, Rcpp::NumericVector starts, Rcpp::NumericVector mus, Rcpp::NumericVector etas, Function var, Function mu_eta, Function linkinv, Function dev_resids, Function valideta, Function validmu, int type, double tol, int maxit);
RcppExport SEXP _glmmr_fast_glm_impl(SEXP XsSEXP, SEXP ysSEXP, SEXP weightssSEXP, SEXP offsetsSEXP, SEXP startsSEXP, SEXP musSEXP, SEXP etasSEXP, SEXP varSEXP, SEXP mu_etaSEXP, SEXP linkinvSEXP, SEXP dev_residsSEXP, SEXP validetaSEXP, SEXP validmuSEXP, SEXP typeSEXP, SEXP tolSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weightss(weightssSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mus(musSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type etas(etasSEXP);
    Rcpp::traits::input_parameter< Function >::type var(varSEXP);
    Rcpp::traits::input_parameter< Function >::type mu_eta(mu_etaSEXP);
    Rcpp::traits::input_parameter< Function >::type linkinv(linkinvSEXP);
    Rcpp::traits::input_parameter< Function >::type dev_resids(dev_residsSEXP);
    Rcpp::traits::input_parameter< Function >::type valideta(validetaSEXP);
    Rcpp::traits::input_parameter< Function >::type validmu(validmuSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_glm_impl(Xs, ys, weightss, offsets, starts, mus, etas, var, mu_eta, linkinv, dev_resids, valideta, validmu, type, tol, maxit));
    return rcpp_result_gen;
END_RCPP
}
// myglm
Rcpp::List myglm(Rcpp::NumericMatrix Xs, Rcpp::NumericVector ys, Rcpp::NumericVector weightss, Rcpp::NumericVector offsets, Rcpp::List family);
RcppExport SEXP _glmmr_myglm(SEXP XsSEXP, SEXP ysSEXP, SEXP weightssSEXP, SEXP offsetsSEXP, SEXP familySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weightss(weightssSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type family(familySEXP);
    rcpp_result_gen = Rcpp::wrap(myglm(Xs, ys, weightss, offsets, family));
    return rcpp_result_gen;
END_RCPP
}
// qscore_impl
double qscore_impl(const arma::vec& resids, arma::vec tr, const arma::vec& xb, const arma::mat& invS, const std::string& family2, bool weight);
RcppExport SEXP _glmmr_qscore_impl(SEXP residsSEXP, SEXP trSEXP, SEXP xbSEXP, SEXP invSSEXP, SEXP family2SEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tr(trSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invS(invSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family2(family2SEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(qscore_impl(resids, tr, xb, invS, family2, weight));
    return rcpp_result_gen;
END_RCPP
}
// permutation_test_impl
arma::vec permutation_test_impl(const arma::vec& resids, const arma::mat& tr_mat, const arma::vec& xb, const arma::mat& invS, const std::string& family2, bool weight, int iter, bool verbose);
RcppExport SEXP _glmmr_permutation_test_impl(SEXP residsSEXP, SEXP tr_matSEXP, SEXP xbSEXP, SEXP invSSEXP, SEXP family2SEXP, SEXP weightSEXP, SEXP iterSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tr_mat(tr_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invS(invSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family2(family2SEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(permutation_test_impl(resids, tr_mat, xb, invS, family2, weight, iter, verbose));
    return rcpp_result_gen;
END_RCPP
}
// confint_search
double confint_search(double start, double b, Rcpp::NumericMatrix Xnull_, Rcpp::NumericVector y_, Rcpp::NumericVector tr_, const arma::mat& new_tr_mat, const arma::vec& xb, const arma::mat& invS, Rcpp::List family, const std::string& family2, int nsteps, bool weight, double alpha, bool verbose);
RcppExport SEXP _glmmr_confint_search(SEXP startSEXP, SEXP bSEXP, SEXP Xnull_SEXP, SEXP y_SEXP, SEXP tr_SEXP, SEXP new_tr_matSEXP, SEXP xbSEXP, SEXP invSSEXP, SEXP familySEXP, SEXP family2SEXP, SEXP nstepsSEXP, SEXP weightSEXP, SEXP alphaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type start(startSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xnull_(Xnull_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tr_(tr_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type new_tr_mat(new_tr_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invS(invSSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type family(familySEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family2(family2SEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(confint_search(start, b, Xnull_, y_, tr_, new_tr_mat, xb, invS, family, family2, nsteps, weight, alpha, verbose));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4binomial_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4binomial_sim_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4binomial_sim_misspec_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gaussian_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gaussian_sim_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gaussian_sim_misspec_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4mcml_binomial_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4mcml_gaussian_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4mcml_poisson_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4poisson_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4poisson_sim_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4poisson_sim_misspec_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_glmmr_dfbeta_stat", (DL_FUNC) &_glmmr_dfbeta_stat, 4},
    {"_glmmr_GradRobustStep", (DL_FUNC) &_glmmr_GradRobustStep, 11},
    {"_glmmr_log_factorial_approx", (DL_FUNC) &_glmmr_log_factorial_approx, 1},
    {"_glmmr_log_mv_gaussian_pdf", (DL_FUNC) &_glmmr_log_mv_gaussian_pdf, 3},
    {"_glmmr_gen_dhdmu", (DL_FUNC) &_glmmr_gen_dhdmu, 3},
    {"_glmmr_fexp", (DL_FUNC) &_glmmr_fexp, 3},
    {"_glmmr_sqexp", (DL_FUNC) &_glmmr_sqexp, 3},
    {"_glmmr_matern", (DL_FUNC) &_glmmr_matern, 3},
    {"_glmmr_bessel1", (DL_FUNC) &_glmmr_bessel1, 2},
    {"_glmmr_blockMat", (DL_FUNC) &_glmmr_blockMat, 1},
    {"_glmmr_genBlockD", (DL_FUNC) &_glmmr_genBlockD, 8},
    {"_glmmr_genD", (DL_FUNC) &_glmmr_genD, 10},
    {"_glmmr_remove_one_many_mat", (DL_FUNC) &_glmmr_remove_one_many_mat, 2},
    {"_glmmr_add_one_mat", (DL_FUNC) &_glmmr_add_one_mat, 3},
    {"_glmmr_d_lik_optim", (DL_FUNC) &_glmmr_d_lik_optim, 14},
    {"_glmmr_l_lik_optim", (DL_FUNC) &_glmmr_l_lik_optim, 10},
    {"_glmmr_f_lik_grad", (DL_FUNC) &_glmmr_f_lik_grad, 20},
    {"_glmmr_f_lik_hess", (DL_FUNC) &_glmmr_f_lik_hess, 20},
    {"_glmmr_f_lik_optim", (DL_FUNC) &_glmmr_f_lik_optim, 20},
    {"_glmmr_mcnr_step", (DL_FUNC) &_glmmr_mcnr_step, 7},
    {"_glmmr_aic_mcml", (DL_FUNC) &_glmmr_aic_mcml, 17},
    {"_glmmr_fast_glm_impl", (DL_FUNC) &_glmmr_fast_glm_impl, 16},
    {"_glmmr_myglm", (DL_FUNC) &_glmmr_myglm, 5},
    {"_glmmr_qscore_impl", (DL_FUNC) &_glmmr_qscore_impl, 6},
    {"_glmmr_permutation_test_impl", (DL_FUNC) &_glmmr_permutation_test_impl, 8},
    {"_glmmr_confint_search", (DL_FUNC) &_glmmr_confint_search, 14},
    {"_rcpp_module_boot_stan_fit4binomial_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4binomial_mod, 0},
    {"_rcpp_module_boot_stan_fit4binomial_sim_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4binomial_sim_mod, 0},
    {"_rcpp_module_boot_stan_fit4binomial_sim_misspec_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4binomial_sim_misspec_mod, 0},
    {"_rcpp_module_boot_stan_fit4gaussian_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gaussian_mod, 0},
    {"_rcpp_module_boot_stan_fit4gaussian_sim_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gaussian_sim_mod, 0},
    {"_rcpp_module_boot_stan_fit4gaussian_sim_misspec_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gaussian_sim_misspec_mod, 0},
    {"_rcpp_module_boot_stan_fit4mcml_binomial_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mcml_binomial_mod, 0},
    {"_rcpp_module_boot_stan_fit4mcml_gaussian_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mcml_gaussian_mod, 0},
    {"_rcpp_module_boot_stan_fit4mcml_poisson_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mcml_poisson_mod, 0},
    {"_rcpp_module_boot_stan_fit4poisson_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4poisson_mod, 0},
    {"_rcpp_module_boot_stan_fit4poisson_sim_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4poisson_sim_mod, 0},
    {"_rcpp_module_boot_stan_fit4poisson_sim_misspec_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4poisson_sim_misspec_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_glmmr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
