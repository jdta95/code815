#include <RcppArmadillo.h>
// #include <losses.h>

// IRLS function for Poisson regression
arma::vec irls_poisson(
  const arma::vec& y,
  const arma::mat& X,
  const arma::vec& beta0,
  double tol = 0.0001,
  int max_iter = 10000,
  bool printing = false
) {
  int n = y.n_elem;
  int p = X.n_cols;
  
  // stopifnot(n == nrow(X))
  // stopifnot(n >= p)
  
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


// Rcpp::List gradient_descent_BB_lsq(
//   const arma::vec& y,
//   const arma::mat& A,
//   const arma::vec& x0,
//   double lambda,
//   double tol = 0.0001,
//   int max.iter = 10000,
//   bool printing = FALSE
// ) {
//   int n = y.n_elem;
//   int p = A.n_cols;
//   // stopifnot(n == nrow(A))
//   // stopifnot(n >= p)
//   
//   arma::mat AA = A.t() * A;
//   arma::mat Ay = A.t() * y;
//   arma::vec grad = AA * x0 - Ay;
//   
//   double loss = loss_ridge(y, A, x0, lambda);
//   arma::vec = grad + 2 * lambda * x0;
//   
//   arma::vec x = x0 - grad;
//   arma::vec prevx = x0;
//   arma::vec prevgrad = grad;
//   arma::vec prevloss = loss;
//   int iter = 0;
//   double diff = inf;
//   
//   
//   
//   
//   
//   
//   
//   
//   
//   
//   
//   
//   return 
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
