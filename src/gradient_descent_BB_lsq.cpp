#include <RcppArmadillo.h>
#include "losses.h"

// gradient descent for least squares with Barzilai-Borwein step size
// [[Rcpp::export]]
Rcpp::List gradient_descent_BB_lsq(
  const arma::vec& y,
  const arma::mat& A,
  const arma::vec& x0,
  double lambda,
  double tol = 0.0001,
  int max_iter = 10000,
  bool printing = false
) {
  arma::uword n = y.n_elem;
  arma::uword p = A.n_cols;
  
  if (n != A.n_rows) {
    Rcpp::stop("y and A must have the same number of rows.");
  }
  
  if (n < p) {
    Rcpp::stop("y must have at least as many rows as columns in A.");
  }

  arma::mat AA = A.t() * A; // p x p
  arma::mat Ay = A.t() * y; // p x 1
  arma::vec grad = AA * x0 - Ay; // p x 1

  double loss = loss_ridge(y, A, x0, lambda);
  grad += 2 * lambda * x0; // p x 1

  arma::vec x = x0 - grad; // p x 1
  arma::vec prevx = x0; // p x 1
  arma::vec prevgrad = grad; // p x 1
  double prevloss = loss;
  int iter = 0;
  double diff = std::numeric_limits<double>::infinity();
  vector<double> diff_rec;
  double gamma;
  
  while (diff > tol && iter < max_iter) {
    grad = AA * x - Ay;
    grad += 2 * lambda * x;
    
    double denom = dot(x - prevx, grad - prevgrad);
    if (denom == 0) {
      gamma = 1.0; // default or some other heuristic
    } else {
      gamma = dot(x - prevx, x - prevx) / denom;
    }
    
    // double gamma = dot(x - prevx, x - prevx) / dot(x - prevx, grad - prevgrad);
    prevx = x;
    prevgrad = grad;
    prevloss = loss;
    x = prevx - gamma * prevgrad;
    loss = loss_ridge(y, A, x, lambda);
    double diff_val = (prevloss - loss) / std::abs(prevloss);
    diff_rec.push_back(diff_val);
    // diff_rec[iter] = (prevloss - loss) / std::abs(prevloss);
    diff = std::abs(diff_val);
    
    // diff
    
    prevloss = loss;
    iter++;
  }
  
  if (printing) {
    Rcpp::Rcout << "Converged in " << iter << " iterations." << std::endl;
  }
  if (iter == max_iter) {
    Rcpp::Rcout << "Did not converge within the maximum number of iterations." << std::endl;
  };
  return( Rcpp::List::create(
    Rcpp::Named("x") = x,
    Rcpp::Named("loss") = loss,
    Rcpp::Named("diff") = diff_rec
    )
  );
}

/*** R
# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()
*/

/*** R
set.seed(123)

n = 1000
x = c(0.3, 0.1, 0.03, 0, 0)
p = length(x)
A = matrix(rnorm(n * p), n, p)
y = A %*% x + rnorm(n)
lambda = 10   #-- lambda is the tuning parameter, not step size
x0 = rnorm(p)


# closed form solution to ridge
solve(t(A) %*% A + 2*lambda*diag(p)) %*% t(A) %*% y

out = gradient_descent_BB_lsq(y, A, x0, lambda, tol = 0.0001, max_iter = 10000, printing = TRUE)
*/
