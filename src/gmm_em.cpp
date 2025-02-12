#include <RcppArmadillo.h>

// look for log sum exp trick
// 
// // my function for EM algorithm for Gaussian mixture model
// // [[Rcpp::export]]
// Rcpp::List gmm_em(
//   arma::mat& X,
//   int k,
//   int max_iter = 10000,
//   double tol = 1e-6
// ) {
//   int n = X.n_rows;
//   int d = X.n_cols;
//   int iter = 0;
//   double diff = tol + 1;
// 
//   arma::vec q = arma::randu<arma::vec>(k);
//   arma::mat mu = arma::randu<arma::mat>(k, d);
//   arma::cube sigma = arma::randu<arma::cube>(d, d, k);
//   arma::mat w(n, k);
//   double pi = 3.14159265358979323846;
//   
//   while (diff > tol && iter < max_iter) {
//     // E-step
//     // Calculate w
//     for (int j = 0; j < k; ++j) {
//       for (int i = 0; i < n; ++i) {
//         w(i, j) = 
//           q(j) * pow(2 * pi, -1 * d / 2) *
//           exp(-0.5 * (X.row(i) - mu.row(j)).t() * inv(sigma.slice(j)) *
//           (X.row(i) - mu.row(j))) / 
//           sqrt(det(sigma.slice(j)));
//       }
//     }
//     
//     for (int i = 0; i < n; ++i) {
//       w.row(i) /= sum(w.row(i));
//     }
//     // M-step
//     for (int j = 0; j < k; ++j) {
//       mu.row(j) = sum(w.col(j).t() * X) / sum(w.col(j));
//       
//       for (int i = 0; i < n; ++i) {
//         sigma.slice(j) += w(i, j) * (X.row(i) - mu.row(j)).t() * (X.row(i) - mu.row(j));
//       }
//       sigma.slice(j) /= sum(w.col(j));
//       
//       q(j) = sum(w.col(j)) / n;
//     }
//     // calculate the norm for mu, sigma, and q combined
//     
//     iter++
//   }  
//   return Rcpp::List::create(Rcpp::Named("mu") = mu,
//                             Rcpp::Named("sigma") = sigma,
//                             Rcpp::Named("q") = q);
// }


// // function for EM algorithm for Gaussian mixture model
// // [[Rcpp::export]]
// Rcpp::List gmm_em_copilot(
//     arma::mat X,
//     int k,
//     int max_iter = 10000,
//     double tol = 1e-6
// ) {
//   int n = X.n_rows;
//   int d = X.n_cols;
// 
//   // initialize parameters
//   arma::mat mu = arma::randu<arma::mat>(k, d);
//   arma::mat sigma = arma::randu<arma::mat>(k, d);
//   arma::vec pi = arma::randu<arma::vec>(k);
//   pi /= sum(pi);
// 
//   // initialize responsibilities
//   arma::mat r(n, k);
// 
//   for (int iter = 0; iter < max_iter; ++iter) {
//     // E-step: compute responsibilities
//     for (int j = 0; j < k; ++j) {
//       r.col(j) = pi(j) * mvnpdf(X, mu.row(j), sigma.row(j));
//     }
//     r.each_row() /= sum(r, 1);
// 
//     // M-step: update parameters
//     for (int j = 0; j < k; ++j) {
//       double N_k = sum(r.col(j));
//       mu.row(j) = sum(r.col(j) % X) / N_k;
//       sigma.row(j) = sum(r.col(j) % (X - mu.row(j)).t() * (X - mu.row(j))) / N_k;
//       pi(j) = N_k / n;
//     }
// 
//     // check convergence
//     if (arma::norm(r - r, "fro") < tol) {
//       break;
//     }
//   }
// 
//   return Rcpp::List::create(Rcpp::Named("mu") = mu,
//                             Rcpp::Named("sigma") = sigma,
//                             Rcpp::Named("pi") = pi,
//                             Rcpp::Named("r") = r);
// }