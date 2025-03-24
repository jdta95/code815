#include <RcppArmadillo.h>
#include "ramadapt.h"

using namespace std;

// Function to compute the log-density of the Beta distribution
double beta_log_density(double x, double alpha, double beta) {
  if (x <= 0.0 || x >= 1.0 || alpha <= 0.0 || beta <= 0.0) {
    return -std::numeric_limits<double>::infinity(); // Log-density is undefined for these values
  }

  // Compute the log B(alpha, beta) using lgamma
  double log_beta_func = std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);

  // Compute the log-density
  double log_density = (alpha - 1.0) * std::log(x) + (beta - 1.0) * std::log(1.0 - x) - log_beta_func;

  return log_density;
}

//[[Rcpp::export]]
arma::vec rwm(double xstart, double alpha, double beta, int mcmc=1000, double lower=0, double upper=1){

  arma::vec samples = arma::zeros(mcmc);
  arma::mat theta_bounds = arma::zeros(1, 2);
  theta_bounds(0, 0) = lower;
  theta_bounds(0, 1) = upper;

  arma::vec theta_cur = xstart * arma::ones(1);

  arma::mat theta_metrop_sd = 0.05 * arma::eye(1, 1);
  RAMAdapt theta_adapt(1, theta_metrop_sd, 0.24);


  for(int i=0; i<mcmc; i++){


    theta_adapt.count_proposal();

    Rcpp::RNGScope scope;
    arma::vec U_update = arma::randn(1);

    arma::vec theta_alt = par_huvtransf_back(par_huvtransf_fwd(
      theta_cur, theta_bounds) + theta_adapt.paramsd * U_update, theta_bounds);

    double curr_logdens = beta_log_density(theta_cur(0), alpha, beta);

    double prop_logdens = beta_log_density(theta_alt(0), alpha, beta);

    // make move
    double jacobian  = calc_jacobian(theta_alt, theta_cur, theta_bounds);
    double logaccept = prop_logdens - curr_logdens + jacobian;

    bool accepted = do_I_accept(logaccept);

    if(accepted){
      theta_cur = theta_alt;
    }

    theta_adapt.update_ratios();

    if(true){
      theta_adapt.adapt(U_update, exp(logaccept), i);
    }

    samples(i) = theta_cur(0);
  }

  return samples;

}