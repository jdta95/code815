#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List Bayes_Poisson(
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& beta0,
    double step,
    const int& mcmc
) {
  // Function to perform Bayesian Poisson regression
  int n = y.n_elem;
  int p = X.n_cols;
  double accept_count = 0;
  
  arma::mat beta_samples = arma::zeros(mcmc, p);
  arma::vec beta_cur = beta0;
  
  for(int i=0; i<mcmc; i++){
    
    arma::vec U_update = arma::randn(p);
    
    arma::vec gradient = arma::zeros(p);
    
    for(int j=0; j<p; j++){
      gradient(j) = arma::as_scalar(y.t() * X.col(j) - X.col(j).t() * exp(X.col(j) * beta_cur(j)));
    }
    
    arma::vec beta_alt = beta_cur + pow(step, 2) / 2 * gradient + step * U_update;
    
    arma::vec lambda_cur = arma::exp(X * beta_cur);
    arma::vec lambda_alt = arma::exp(X * beta_alt);
    
    double curr_logdens = arma::as_scalar(y.t() * X * beta_cur) - arma::accu(lambda_cur);
    double prop_logdens = arma::as_scalar(y.t() * X * beta_alt) - arma::accu(lambda_alt);
    
    double logaccept = prop_logdens - curr_logdens;
    bool accepted = exp(logaccept) > arma::randu();
    
    if(accepted){
      beta_cur = beta_alt;
      accept_count++;
    }

    beta_samples.row(i) = beta_cur.t();
  }
  
  Rcpp::List output;
  output["beta_samples"] = beta_samples;
  output["acceptance_rate"] = accept_count / mcmc;
  
  return output;
}

