// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gradient_descent_BB_lsq
Rcpp::List gradient_descent_BB_lsq(const arma::vec& y, const arma::mat& A, const arma::vec& x0, double lambda, double tol, int max_iter, bool printing);
RcppExport SEXP _code815_gradient_descent_BB_lsq(SEXP ySEXP, SEXP ASEXP, SEXP x0SEXP, SEXP lambdaSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP printingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type printing(printingSEXP);
    rcpp_result_gen = Rcpp::wrap(gradient_descent_BB_lsq(y, A, x0, lambda, tol, max_iter, printing));
    return rcpp_result_gen;
END_RCPP
}
// irls_poisson
arma::vec irls_poisson(const arma::vec& y, const arma::mat& X, const arma::vec& beta0, double tol, int max_iter, bool printing);
RcppExport SEXP _code815_irls_poisson(SEXP ySEXP, SEXP XSEXP, SEXP beta0SEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP printingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type printing(printingSEXP);
    rcpp_result_gen = Rcpp::wrap(irls_poisson(y, X, beta0, tol, max_iter, printing));
    return rcpp_result_gen;
END_RCPP
}
// loss_ridge
double loss_ridge(const arma::vec& y, const arma::mat& A, const arma::vec& x, double lambda);
RcppExport SEXP _code815_loss_ridge(SEXP ySEXP, SEXP ASEXP, SEXP xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_ridge(y, A, x, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_code815_gradient_descent_BB_lsq", (DL_FUNC) &_code815_gradient_descent_BB_lsq, 7},
    {"_code815_irls_poisson", (DL_FUNC) &_code815_irls_poisson, 6},
    {"_code815_loss_ridge", (DL_FUNC) &_code815_loss_ridge, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_code815(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
