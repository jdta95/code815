#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat VI_normal(
    const arma::vec& X,
    double mumu0,
    double musigma0,
    double A0,
    double B0,
    double tol = 1e-6
) {
  // Initialize parameters
  double n = X.n_elem;
  double inf = std::numeric_limits<double>::infinity();
  double Xbar = arma::mean(X);
  double mumu = mumu0;
  double musigma = musigma0;
  double A = A0 + n / 2;
  double B = A0 / B0;
  double logp_cur = inf;
  
  // create output matrix
  arma::mat output(5, 1);
  output(0, 0) = mumu;
  output(1, 0) = musigma;
  output(2, 0) = A;
  output(3, 0) = B;
  output(4, 0) = inf;
  
  // loop as long as div > tol
  double div = inf;
  
  while (div > tol) {
    // update musigma
    musigma = 1 / (n * A / B + 1 / musigma0);
    
    // update mumu
    mumu = (n * Xbar * A / B + mumu0 / musigma0) * musigma;
    
    // update B
    B = B0 + 0.5 * arma::as_scalar((X - mumu * arma::ones(n)).t() * (X - mumu * arma::ones(n))) +
      n * musigma;
    
    double logp_new = 0.5 * std::log(musigma / musigma0) -
      (pow(mumu - mumu0, 2) - musigma) / (2 * musigma0) -
      A * std::log(B);
    
    div = std::abs(logp_new - logp_cur);
    
    logp_cur = logp_new;
    
    // add new column to output with updated parameters
    output.insert_cols(output.n_cols, 1);
    output(0, output.n_cols - 1) = mumu;
    output(1, output.n_cols - 1) = musigma;
    output(2, output.n_cols - 1) = A;
    output(3, output.n_cols - 1) = B;
    output(4, output.n_cols - 1) = div;
  }
  
  return output;
}

