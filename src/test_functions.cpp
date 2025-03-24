// #include <RcppArmadillo.h>
// 
// arma::mat inv_Chol2(
//     arma::mat A
// ) {
//   // try {
//   A = 0.5 * (A + A.t()); // Ensure symmetry
//   A += arma::eye(A.n_rows, A.n_cols) * 1e-6; // Add a small value to the diagonal for numerical stability
//   arma::mat R = arma::chol(A);
//   arma::mat Rinv = arma::inv(R);
//   arma::mat Ainv = Rinv * Rinv.t();
//   return Ainv;
//   // } catch (std::exception& e) {
//   //   A.save("matrix_debug.csv", arma::csv_ascii);
//   //   Rcpp::stop("Error in inv_Chol: %s", e.what());
//   // }
// }
// 
// double logdet2(
//     arma::mat A
// ) {
//   // try {
//   A = 0.5 * (A + A.t()); // Ensure symmetry
//   A += arma::eye(A.n_rows, A.n_cols) * 1e-6; // Add a small value to the diagonal for numerical stability
//   arma::mat R = arma::chol(A);
//   double logdet = 2 * arma::accu(arma::log(arma::diagvec(R)));
//   return logdet;
//   // } catch (std::exception& e) {
//   //   A.save("matrix_debug.csv", arma::csv_ascii);
//   //   Rcpp::stop("Error in logdet: %s", e.what());
//   // }
// }
// 
// //[[Rcpp::export]]
// double test_GP_log_density(
//     arma::mat X,
//     int n,
//     double sigmasq,
//     arma::vec w,
//     arma::vec phi
// ) {
//   arma::mat C_phi = arma::zeros(n, n);
//   
//   for(int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       C_phi(i, j) = arma::as_scalar(
//         exp(
//           -0.5 * (X.row(i) - X.row(j)) * arma::diagmat(1 / phi) *
//             (X.row(i) - X.row(j)).t()
//         )
//       );
//     }
//   }
//   
//   arma::mat Sigma = sigmasq * C_phi;
//   
//   double log_likelihood =  -0.5 * logdet2(Sigma) -
//     0.5 * arma::as_scalar(w.t() * inv_Chol2(Sigma) * w);
//   
//   return log_likelihood;
// }
// 
// double test_GP_log_density(
//     arma::vec Y,
//     arma::mat X,
//     int n,
//     int p,
//     double sigmasq,
//     double tausq,
//     arma::vec w,
//     arma::vec beta,
//     arma::mat V,
//     arma::vec phi
// ) {
//   // print checkpoint message
//   // Rcpp::Rcout << "Check 1 good. " << std::endl;
// 
//   arma::mat C_phi = arma::zeros(n, n);
// 
//   for(int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       C_phi(i, j) = arma::as_scalar(
//         exp(
//           -0.5 * (X.row(i) - X.row(j)) * arma::diagmat(1 / phi) *
//             (X.row(i) - X.row(j)).t()
//         )
//       );
//     }
//   }
// 
//   arma::mat Sigma = X * V * X.t() + sigmasq * C_phi + tausq * arma::eye(n, n);
//   arma::vec mu = X * beta + w;
// 
//   double log_likelihood =  -0.5 * logdet2(Sigma) -
//     0.5 * arma::as_scalar((Y - mu).t() * inv_Chol2(Sigma) * (Y - mu));
// 
//   // Rcpp::Rcout << "Check 2 good. " << std::endl;
// 
//   return log_likelihood;
// }