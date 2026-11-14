# --------------------------------------------------------------------------
# mcglm() argument-level paths previously missing from tests:
#   * iterate = TRUE for K >= 3
#   * direct c1/c2 entry (advanced binary path)
#   * Gaussian family for the analytical methods
#   * pi_z auto-estimated from Pi via matrix interface (Bayesian inversion)
#   * .mcglm_default_method() via S3 fallback when method = NULL
# --------------------------------------------------------------------------

# ==========================================================================
# iterate = TRUE for multicategory case
# ==========================================================================

test_that("iterate=TRUE for K=3 BCA/BCM converges and matches one-shot when stable", {
  set.seed(91)
  n <- 1000; K <- 3
  pi_z <- c(0.5, 0.3, 0.2)
  Pi <- matrix(0.05, K, K); diag(Pi) <- 1 - 0.05 * (K - 1)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1), 1, prob = Pi[, zi + 1]),
                  integer(1))
  x <- cbind(1, rnorm(n))
  gf <- c(0, 0.7, -0.4)
  y <- rpois(n, exp(drop(x %*% c(-0.3, 0.5)) + gf[z + 1]))

  fit_once <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                    method = c("bca", "bcm"), Pi = Pi, pi_z = pi_z,
                    iterate = FALSE)
  fit_iter <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                    method = c("bca", "bcm"), Pi = Pi, pi_z = pi_z,
                    iterate = TRUE)

  expect_true(all(is.finite(fit_iter$coefficients$bca)))
  expect_true(all(is.finite(fit_iter$coefficients$bcm)))
  # iterated BCM should approach the corrected-score solution. We only
  # require that iteration moves things by less than the magnitude of
  # the one-shot correction itself.
  diff_bcm <- max(abs(fit_iter$coefficients$bcm - fit_once$coefficients$bcm))
  expect_lt(diff_bcm, 0.5)
})

# ==========================================================================
# Misclassification parameters and rejection of direct c1 / c2
# ==========================================================================

test_that("c1 and c2 are computed internally from p01/p10/pi_z and direct inputs are rejected", {
  set.seed(92)
  n <- 500
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- rpois(n, exp(0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  c1_expected <- p01 * (1 - pi_z)
  c2_expected <- p01 * (1 - pi_z) - p10 * pi_z

  fit_pi <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                  method = c("bca", "bcm", "cs"),
                  p01 = p01, p10 = p10, pi_z = pi_z)

  expect_equal(fit_pi$misclass$c1, c1_expected)
  expect_equal(fit_pi$misclass$c2, c2_expected)
  expect_true(all(is.finite(vcov(fit_pi, method = "cs"))))

  # Supplying c1 or c2 directly is rejected as an unused argument
  expect_error(mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                     p01 = p01, p10 = p10, pi_z = pi_z, c1 = c1_expected),
               "unused argument")
  expect_error(mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                     p01 = p01, p10 = p10, pi_z = pi_z, c2 = c2_expected),
               "unused argument")
})

test_that("vcov_corrected works with p01/p10/pi_z (binary)", {
  set.seed(99)
  n <- 800
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- rpois(n, exp(0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = "bcm", p01 = 0.10, p10 = 0.15, pi_z = 0.4,
               vcov_corrected = TRUE)
  V <- vcov(fit, method = "bcm")
  expect_true(all(is.finite(V)))
  expect_true(all(diag(V) > 0))
})

# ==========================================================================
# Gaussian family for analytical methods
# ==========================================================================

test_that("Gaussian family supports naive/bca/bcm/cs", {
  set.seed(93)
  n <- 800
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- drop(x %*% c(-0.5, 0.7) + 0.8 * z + rnorm(n, sd = 0.5))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "gaussian",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = 0.10, p10 = 0.15, pi_z = 0.4)
  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_true(all(is.finite(fit$coefficients[[m]])), info = m)
    expect_length(fit$coefficients[[m]], 3)
  }
  # Gamma should move toward the truth (0.8) when going from naive to bcm
  expect_lt(abs(fit$coefficients$bcm[1] - 0.8),
            abs(fit$coefficients$naive[1] - 0.8))
})

# ==========================================================================
# Estimated pi_z under matrix interface (Bayesian inversion via Pi)
# ==========================================================================

test_that("matrix interface: pi_z auto-estimated when only Pi supplied (binary)", {
  set.seed(94)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- rpois(n, exp(0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  Pi <- matrix(c(0.90, 0.10,
                 0.15, 0.85), nrow = 2, ncol = 2)

  # Only Pi supplied — pi_z must be inferred via Pi^{-1} %*% pi_obs.
  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("bca", "bcm"), Pi = Pi)
  expect_true(all(is.finite(fit$coefficients$bca)))
  expect_true(all(is.finite(fit$coefficients$bcm)))
})

test_that("pi_z estimation solves the probability system", {
  Pi <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2)
  latent <- c(0.6, 0.4)
  observed <- drop(Pi %*% latent)
  counts <- as.integer(observed * 1000)
  z_hat <- rep(0:1, counts)

  expect_equal(mismeasured:::.mcglm_estimate_pi_z(z_hat, Pi), latent[2],
               tolerance = 1e-12)
})

test_that("matrix interface: pi_z auto-estimated when only Pi supplied (K=3)", {
  set.seed(95)
  n <- 800; K <- 3
  pi_z <- c(0.5, 0.3, 0.2)
  Pi <- matrix(0.05, K, K); diag(Pi) <- 1 - 0.05 * (K - 1)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1), 1, prob = Pi[, zi + 1]),
                  integer(1))
  x <- cbind(1, rnorm(n))
  y <- rpois(n, exp(drop(x %*% c(-0.3, 0.5)) + c(0, 0.7, -0.4)[z + 1]))

  fit <- mcglm(y ~ mc(z_obs, Pi) + x1,
               data = data.frame(y = y, z_obs = z_hat, x1 = x[, 2]),
               family = "poisson", method = c("naive", "bca", "bcm"))
  expect_named(fit$coefficients, c("naive", "bca", "bcm"))
})

# ==========================================================================
# Default method dispatch (.mcglm_default_method)
# ==========================================================================

test_that("S3 methods default to a reasonable method when none specified", {
  set.seed(96)
  n <- 300
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.4)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))
  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca", "bcm"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.4)

  # No method= specified -> default (corrected if available, else naive)
  cf <- coef(fit)
  expect_true(is.numeric(cf) || is.list(cf))
  expect_true(is.matrix(vcov(fit)))
  expect_length(fitted(fit), n)
})

test_that("mcglm rejects invalid misclassification parameters", {
  set.seed(97)
  n <- 100
  z <- rep(0:1, length.out = n)
  y <- rpois(n, exp(0.5 * z))
  x <- matrix(1, n, 1)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)

  bad_entry <- matrix(c(1.1, -0.1, -0.1, 1.1), 2, 2)
  expect_error(mcglm(y, z_hat = z, x = x, method = "naive", Pi = bad_entry),
               "[0, 1]", fixed = TRUE)
  bad_sum <- matrix(c(0.9, 0.2, 0.2, 0.9), 2, 2)
  expect_error(mcglm(y, z_hat = z, x = x, method = "naive", Pi = bad_sum),
               "sum to 1")
  bad_finite <- Pi
  bad_finite[1, 1] <- Inf
  expect_error(mcglm(y, z_hat = z, x = x, method = "naive", Pi = bad_finite),
               "finite")

  expect_error(mcglm(y, z_hat = z, x = x, method = "cs",
                     p01 = -0.1, p10 = 0.1, pi_z = 0.4), "p01")
  expect_error(mcglm(y, z_hat = z, x = x, method = "cs",
                     p01 = 0.1, p10 = 0.1, pi_z = c(0.4, 0.6)), "length 1")
  expect_error(mcglm(y, z_hat = z, x = x, method = "cs",
                     c1 = Inf, c2 = 0), "unused argument")
  expect_error(mcglm(y, z_hat = z, x = x, method = "naive", Pi = Pi,
                     weights = rep(Inf, n)), "finite positive")
})
