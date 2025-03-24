# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

set.seed(123)

n = 200
p = 2

beta = runif(p + 1)
X = matrix(rnorm(n * p), n, p)
X = cbind(rep(1, n), X)
y = rpois(n, exp(X %*% beta))

# get acceptance rates and ESS with different step values
step_values = seq(0.01, 0.02, by = 0.001)
accept_ess = sapply(step_values, function(step) {
  out = code815::Bayes_Poisson(
    y = y,
    X = X,
    beta0 = rep(0, p + 1),
    step = step,
    mcmc = 10000
  )
  ess = apply(out$beta_samples, 2, function(x) {
    coda::effectiveSize(x)
  })
  return(c(out$acceptance_rate, ess))
})

# step = 0.013 yielded closest acceptance rate to 0.57
# optimal acceptance rate does NOT result in optimal ESS
cbind(step_values, t(accept_ess))

# test the function
out = code815::Bayes_Poisson(
  y = y,
  X = X,
  beta0 = rep(0, p + 1),
  step = 0.013,
  mcmc = 20000
)

# acceptance rate
out$acceptance_rate

# trace plot for beta
plot(out$beta_samples[, 1], type = "l", main = "Trace plot for beta[1]")
plot(out$beta_samples[, 2], type = "l", main = "Trace plot for beta[2]")
plot(out$beta_samples[, 3], type = "l", main = "Trace plot for beta[3]")

# beta estimates
beta_estimates = apply(out$beta_samples, 2, mean)

# ESS
ess = apply(out$beta_samples, 2, function(x) {
  coda::effectiveSize(x)
})
