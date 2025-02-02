#include <RcppArmadillo.h>
// #include <losses.h>

// IRLS function for Poisson regression
// [[Rcpp::export]]
arma::vec irls_poisson(
  const arma::vec& y,
  const arma::mat& X,
  const arma::vec& beta0,
  double tol = 0.0001,
  int max_iter = 10000,
  bool printing = false
) {
  arma::uword n = y.n_elem;
  arma::uword p = X.n_cols;
  
  if (n != X.n_rows) {
    Rcpp::stop("y and X must have the same number of rows.");
  }
  
  if (n < p) {
    Rcpp::stop("y must have at least as many rows as columns in X.");
  }
  
  arma::vec beta = beta0;
  
  for (int iter = 0; iter < max_iter; iter++) {
    arma::vec mu = exp(X * beta);
    arma::vec z = X * beta + (y - mu) / mu;
    arma::mat W = diagmat(mu);
    arma::vec beta_new = solve(W * X, W * z);
    
    if (norm(beta_new - beta, "inf") < tol) {
      if (printing) {
        Rcpp::Rcout << "Converged in " << iter << " iterations." << std::endl;
      }
      return beta_new;
    }
    
    beta = beta_new;
  }
  
  if (printing) {
    Rcpp::Rcout << "Did not converge within the maximum number of iterations." << std::endl;
  }
  
  return beta;
}

// compile package and load in library for testing
/*** R
# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()
*/

// test the irls_poisson function
/*** R
# simulate Poisson distributed data for a Poisson regression
set.seed(123)
n = 1000
p = 4
X = cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
beta = rnorm(p)
y = rpois(n, lambda = exp(X %*% beta))
beta0 = rep(0, p)

# test irls_poisson function
beta_hat = code815::irls_poisson(y, X, beta0 = rep(0, p))
*/
