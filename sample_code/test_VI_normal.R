# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

set.seed(123)

n = 10000
mu = rnorm(1)
sigma = rnorm(1, 1)

X = rnorm(n, mu, sqrt(sigma))

output = code815::VI_normal(
  X = X,
  mumu0 = 0,
  musigma0 = 1,
  A0 = 2,
  B0 = 1,
  tol = 1e-10
)

output

final_mumu = output[1, ncol(output)]
final_musigma = output[2, ncol(output)]
final_A = output[3, ncol(output)]
final_B = output[4, ncol(output)]
final_divergence = output[5, ncol(output)]

# estimated mu
final_mumu

# estimated mu variance
final_musigma

# estimated sigma
final_B / (final_A - 1)

# estimated sigma variance
final_B^2 / ((final_A - 1)^2 * (final_A - 2))

# percent error of mu
100 * (final_mumu - mu) / mu

# percent error of sigma
100 * (final_B / (final_A - 1) - sigma) / sigma
