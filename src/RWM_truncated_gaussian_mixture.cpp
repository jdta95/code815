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
arma::vec RWM_truncated_gaussian_mixture(
    int n,
    double x0,
    double mu1,
    double mu2,
    double sigma1,
    double sigma2,
    double w1,
    double w2,
    double qsigma,
    double radius,
    double seed
) {
  // set seed
  set_seed(seed);
  // initialize vector to hold samples
  arma::vec x_vec(n);
  // transform from [-radius, radius] support to [0, 1] support
  double x = (x0 + radius) / (2 * radius);
  // record initial value
  x_vec(0) = x;
  // loop over samples
  for (int i = 1; i < n; i++) {
    // transform to logit space
    double z = log(x) - log(1 - x);
    // sample new z from N(z, qsigma^2)
    double z_new = z + qsigma * arma::randn();
    // transform back to new z to new x
    double x_new = 1 / (1 + exp(-1 * z_new));
    // compute acceptance probability alpha
    double q1 = x_new / x;
    double q2 = (1 - x_new) / (1 - x);
    double q_ratio = q1 * q2;
    double xp = x * 2 * radius - radius;
    double xp_new = x_new * 2 * radius - radius;
    double p_ratio = (w1 * arma::normpdf(xp_new, mu1, sigma1) + w2 * arma::normpdf(xp_new, mu2, sigma2)) /
      (w1 * arma::normpdf(xp, mu1, sigma1) + w2 * arma::normpdf(xp, mu2, sigma2));
    double alpha = p_ratio * q_ratio;
    // accept or reject new value
    double u = arma::randu();
    if (u < alpha) {
      x = x_new;
    }
    // record x
    x_vec(i) = x;
  }
  // transform from [0, 1] to [-radius, radius]
  x_vec = x_vec * 2 * radius - radius;
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
x = RWM_truncated_gaussian_mixture(
  n = 10000,
  x0 = 3,
  w1 = 0.3,
  mu1 = -3,
  sigma1 = 1,
  w2 = 0.7,
  mu2 = 3,
  sigma2 = 1,
  qsigma = 1,
  radius = 4,
  seed = 123
)

hist(x)
plot(x)
*/
