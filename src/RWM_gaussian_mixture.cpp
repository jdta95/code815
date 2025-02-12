#include <RcppArmadillo.h>
#include "set_seed.h"

// // set seed
// // [[Rcpp::export]]
// void set_seed(double seed) {
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(std::floor(std::fabs(seed)));
// }

// use random walk metropolis with Gaussian proposal to sample from mixture of Gaussians
// [[Rcpp::export]]
arma::vec RWM_gaussian_mixture(
  int n,
  double x0,
  double mu1,
  double mu2,
  double sigma1,
  double sigma2,
  double w1,
  double w2,
  double qsigma,
  double seed
) {
  // set seed
  set_seed(seed);
  // initialize vector to hold samples
  arma::vec x_vec(n);
  // set initial value
  double x = x0;
  // record initial value
  x_vec(0) = x0;
  // loop over samples
  for (int i = 1; i < n; i++) {
    // sample from normal distribution with variance qsigma^2
    double x_new = x + qsigma * arma::randn();
    // compute acceptance ratio
    double alpha = (w1 * arma::normpdf(x_new, mu1, sigma1) + w2 * arma::normpdf(x_new, mu2, sigma2)) /
      (w1 * arma::normpdf(x, mu1, sigma1) + w2 * arma::normpdf(x, mu2, sigma2));
    // accept or reject new value
    double u = arma::randu();
    if (u < alpha) {
      x = x_new;
    }
    // record x
    x_vec(i) = x;
  }
  return x_vec;
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
x = RWM_gaussian_mixture(
  n = 10000,
  x0 = 0,
  w1 = 0.3,
  mu1 = -3,
  sigma1 = 1,
  w2 = 0.7,
  mu2 = 3,
  sigma2 = 1,
  qsigma = 1,
  seed = 123
)

hist(x)
plot(x)
*/

