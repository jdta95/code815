# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

RWM_gaussian_mixture <- function(n, x0, mu1, mu2, sigma1, sigma2, w1, w2, qsigma, seed) {
    .Call(`_code815_RWM_gaussian_mixture`, n, x0, mu1, mu2, sigma1, sigma2, w1, w2, qsigma, seed)
}

RWM_truncated_gaussian_mixture <- function(n, x0, mu1, mu2, sigma1, sigma2, w1, w2, qsigma, radius, seed) {
    .Call(`_code815_RWM_truncated_gaussian_mixture`, n, x0, mu1, mu2, sigma1, sigma2, w1, w2, qsigma, radius, seed)
}

gradient_descent_BB_lsq <- function(y, A, x0, lambda, tol = 0.0001, max_iter = 10000L, printing = FALSE) {
    .Call(`_code815_gradient_descent_BB_lsq`, y, A, x0, lambda, tol, max_iter, printing)
}

irls_poisson <- function(y, X, beta0, tol = 0.0001, max_iter = 10000L, printing = FALSE) {
    .Call(`_code815_irls_poisson`, y, X, beta0, tol, max_iter, printing)
}

loss_ridge <- function(y, A, x, lambda) {
    .Call(`_code815_loss_ridge`, y, A, x, lambda)
}

