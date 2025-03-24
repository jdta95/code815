# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

set.seed(123)

n = 200
p = 2

X = matrix(rnorm(n * p, 0, 10), nrow = n, ncol = p)

V = diag(1, p)
beta = MASS::mvrnorm(1, rep(0, p), V)

a_s = 2
b_s = 4
sigmasq = 4

l = 0.01
u = 10
phi = c(2, 8)

# Generate the covariance matrix with ARD kernel
C_phi = matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    C_phi[i, j] = exp(-0.5 * t(X[i,] - X[j,]) %*% diag(1 / phi) %*% (X[i,] - X[j,]))
  }
}

w = MASS::mvrnorm(1, rep(0, n), sigmasq * C_phi)

a_t = 2
b_t = 0.01
tausq = 0.01
epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))

Y = as.vector(X %*% beta + w + epsilon)

post_samples = code815::GP_Gibbs(
  Y,
  X,
  sigmasq_start = 1,
  tausq_start = 1,
  phi_start = rep(l + (u - l) / 2, p),
  w_start = rep(0, n),
  beta_start = rep(0, p),
  V,
  a_s,
  b_s,
  a_t,
  b_t,
  lower = rep(l, p),
  upper = rep(u, p),
  mcmc = 1000
)

start = 500
end = 1000

plot(post_samples$beta[start:end, 1], type = "l", main = "Posterior samples of beta 1", ylab = "beta 1", xlab = "Iteration")
plot(post_samples$beta[start:end, 2], type = "l", main = "Posterior samples of beta 2", ylab = "beta 2", xlab = "Iteration")
plot(rowMeans(post_samples$w)[start:end], type = "l", main = "Posterior samples of w means", ylab = "w", xlab = "Iteration")
plot(post_samples$sigmasq[start:end], type = "l", main = "Posterior samples of sigmasq", ylab = "sigmasq", xlab = "Iteration")
plot(post_samples$tausq[start:end], type = "l", main = "Posterior samples of tausq", ylab = "tausq", xlab = "Iteration")
plot(post_samples$phi[start:end, 1], type = "l", main = "Posterior samples of phi 1", ylab = "phi 1", xlab = "Iteration")
plot(post_samples$phi[start:end, 2], type = "l", main = "Posterior samples of phi 2", ylab = "phi 2", xlab = "Iteration")

# show true values and posterior estimates of beta, mean w, sigma^2, tau^2, phi
beta[1]
mean(post_samples$beta[start:end, 1])

beta[2]
mean(post_samples$beta[start:end, 2])

mean(w)
mean(rowMeans(post_samples$w[start:end, ]))

sigmasq
mean(post_samples$sigmasq[start:end])

tausq
mean(post_samples$tausq[start:end])

phi[1]
mean(post_samples$phi[start:end, 1])

phi[2]
mean(post_samples$phi[start:end, 2])
