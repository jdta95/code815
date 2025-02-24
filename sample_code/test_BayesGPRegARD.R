# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

set.seed(123)

alpha_s = 1
beta_s = 1
alpha_t = 1
beta_t = 1
p = 5
l = 0.00001
u = 2

sigmasq = 1 / rgamma(1, shape = alpha_s, shape = beta_s)
tausq = 1 / rgamma(1, shape = alpha_t, shape = beta_t)
phi = runif(p, min = l, max = u)

# Generate data
n = 100
X = matrix(runif(n * p), nrow = n, ncol = p)

# Generate the covariance matrix with ARD kernel
K = matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    K[i, j] = t(X[i,] - X[j,]) %*% diag(1 / phi) %*% (X[i,] - X[j,])
  }
}

C = X %*% diag(1000000, p) %*% t(X) + sigmasq * exp(-0.5 * K) + tausq * diag(1, n)

Y = rnorm(n, rep(0, n), C)

post_samples = code815::BayesGPRegARD(
  Y,
  X,
  sigmasq_start = 1,
  tausq_start = 1,
  phi_start = rep(1, p),
  alpha_s = 1,
  beta_s = 1,
  alpha_t = 1,
  beta_t = 1,
  lower = rep(l, p + 2),
  upper = c(Inf, Inf, rep(u, p)),
  mcmc = 1000
)

plot(post_samples[, 1], type = "l", main = "Posterior samples of sigmasq", ylab = "sigmasq", xlab = "Iteration")
plot(post_samples[, 2], type = "l", main = "Posterior samples of tausq", ylab = "tausq", xlab = "Iteration")
plot(post_samples[, 3], type = "l", main = "Posterior samples of phi", ylab = "phi1", xlab = "Iteration")
plot(post_samples[, 4], type = "l", main = "Posterior samples of phi", ylab = "phi2", xlab = "Iteration")
plot(post_samples[, 5], type = "l", main = "Posterior samples of phi", ylab = "phi3", xlab = "Iteration")
plot(post_samples[, 6], type = "l", main = "Posterior samples of phi", ylab = "phi4", xlab = "Iteration")
plot(post_samples[, 7], type = "l", main = "Posterior samples of phi", ylab = "phi5", xlab = "Iteration")

estimates = apply(post_samples[200:1000,], 2, mean)
estimates
