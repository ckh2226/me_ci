# Function to simulate data- choose variance of errors, sample size, and coefficients for model (to simulate CI magnitude)
sim_data = function(sigmaU, n, alpha1, beta1) {
  # Error-free exposure
  X = rnorm(n = n, mean = 1.8, sd = 1)

  # Error-free fractional rank
  R = (rank(X) - 1) / n + 1 / (2 * n)

  # Health outcome | Fractional rank of true proximity
  eps = rnorm(n = n, mean = 0, sd = 1)
  Y = alpha1 + beta1 * R + eps

  # Errors and error-prone exposure
  U = rnorm(n = n, mean = 0, sd = sigmaU)
  Xstar = X + U

  # Error-prone fractional rank
  Rstar = (rank(Xstar) - 1) / n + 1 / (2 * n)
  w = Rstar - R

  # Return complete dataset, including mean of Y
  data.frame(Y, X, R, Z, U, Xstar, Rstar, w, Y_bar = mean(Y))
}
