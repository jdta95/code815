#include <RcppArmadillo.h>
#include "ramadapt.h"

using namespace std;

// // Function to compute the log-density of the Beta distribution
// double beta_log_density(double x, double alpha, double beta) {
//   if (x <= 0.0 || x >= 1.0 || alpha <= 0.0 || beta <= 0.0) {
//     return -std::numeric_limits<double>::infinity(); // Log-density is undefined for these values
//   }
//   
//   // Compute the log B(alpha, beta) using lgamma
//   double log_beta_func = std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);
//   
//   // Compute the log-density
//   double log_density = (alpha - 1.0) * std::log(x) + (beta - 1.0) * std::log(1.0 - x) - log_beta_func;
//   
//   return log_density;
// }

// double invgamma_log_density(double x, double alpha, double beta) {
//   if (x <= 0.0 || alpha <= 0.0 || beta <= 0.0) {
//     return -std::numeric_limits<double>::infinity(); // Log-density is undefined for these values
//   }
//   
//   // Compute the log-density
//   double log_density = (-1 * alpha - 1.0) * std::log(x) - beta / x - std::lgamma(alpha) + alpha * std::log(beta);
//   
//   return log_density;
// }

double GP_log_density(
  arma::vec Y,
  arma::mat X,
  arma::vec theta,
  int p
) {
  // Extract parameters from theta
  double sigmasq = theta(0);
  double tausq = theta(1);
  arma::vec phi = theta.subvec(2, p + 1);
  arma::mat C = arma::zeros(X.n_rows, X.n_rows);
  
  // Compute the covariance matrix
  for(int i = 0; i < X.n_rows; i++) {
    for (int j = 0; j < X.n_rows; j++) {
      C(i, j) = arma::as_scalar(exp(-0.5 * (X.row(i) - X.row(j)) * arma::diagmat(1 / phi) * (X.row(i) - X.row(j)).t()));
    }
  }
  
  // Add noise variance
  C += tausq * arma::eye(X.n_rows, X.n_rows);
  
  C += X * arma::eye(X.n_cols, X.n_cols) * 1000000 * X.t();
  // C += X * X.t();
  
  // Compute the log-likelihood
  double log_likelihood = -0.5 * arma::as_scalar(Y.t() * arma::solve(C, Y)) - 0.5 * std::log(arma::det(C)) - 0.5 * X.n_rows * std::log(2 * M_PI);
  
  return log_likelihood;
}

//[[Rcpp::export]]
arma::mat BayesGPRegARD(
    arma::vec Y,
    arma::mat X,
    double sigmasq_start,
    double tausq_start,
    arma::vec phi_start,
    double alpha_s,
    double beta_s,
    double alpha_t,
    double beta_t,
    arma::vec lower,
    arma::vec upper,
    int mcmc=1000
){
  
  int p = phi_start.n_elem;
  arma::mat samples = arma::zeros(mcmc, p + 2);
  arma::mat theta_bounds = arma::zeros(p + 2, 2);
  theta_bounds.col(0) = lower;
  theta_bounds.col(1) = upper;
  
  arma::vec sigmasq_cur = sigmasq_start * arma::eye(1, 1);
  arma::vec tausq_cur = tausq_start * arma::eye(1, 1);
  arma::vec phi_cur = phi_start;
  arma::vec theta_cur = arma::join_cols(sigmasq_cur, tausq_cur, phi_cur);
  
  arma::mat theta_metrop_sd = 0.05 * arma::eye(p + 2, p + 2);
  RAMAdapt theta_adapt(p + 2, theta_metrop_sd, 0.24);
  
  for(int i=0; i<mcmc; i++){
    
    
    theta_adapt.count_proposal();
    
    Rcpp::RNGScope scope;
    arma::vec U_update = arma::randn(p + 2);
    
    arma::vec theta_alt = par_huvtransf_back(par_huvtransf_fwd(
      theta_cur, theta_bounds) + theta_adapt.paramsd * U_update, theta_bounds);
    
    double curr_logdens_sigmasq = invgamma_logdens(theta_cur(0), alpha_s, beta_s);
    double curr_logdens_tausq = invgamma_logdens(theta_cur(1), alpha_t, beta_t);
    double curr_logdens_GP = GP_log_density(Y, X, theta_cur, p);
    
    double curr_logdens = curr_logdens_sigmasq + curr_logdens_tausq + curr_logdens_GP;
    
    double prop_logdens_sigmasq = invgamma_logdens(theta_alt(0), alpha_s, beta_s);
    double prop_logdens_tausq = invgamma_logdens(theta_alt(1), alpha_t, beta_t);
    double prop_logdens_GP = GP_log_density(Y, X, theta_alt, p);
    double prop_logdens = prop_logdens_sigmasq + prop_logdens_tausq + prop_logdens_GP;
    
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
    
    samples.row(i) = theta_cur.t();
  }
  
  return samples;
  
}