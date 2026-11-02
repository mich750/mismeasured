# Tests for the AKN98 corrected-score method (cs_akn).

# ---- helpers ----

gen_binary <- function(n = 4000, gamma = 0.8, alpha = c(-0.5, 0.7),
                       p01 = 0.1, p10 = 0.1, prev = 0.4, family,
                       seed = 42) {
  set.seed(seed)
  Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)
  z  <- rbinom(n, 1, prev)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  x1  <- rnorm(n)
  eta <- gamma * z + alpha[1] + alpha[2] * x1
  y <- switch(family,
              poisson  = rpois(n, exp(eta)),
              binomial = rbinom(n, 1, plogis(eta)),
              gaussian = eta + rnorm(n, 0, 0.4))
  list(df = data.frame(y = y, z = factor(z_hat), x1 = x1),
       Pi = Pi, truth = c(gamma, alpha))
}

gen_multi <- function(n = 4000, K = 3,
                      gamma = c(1.0, -0.9), alpha = c(0.8, -0.7),
                      seed = 7) {
  set.seed(seed)
  Pi <- matrix(0.1, K, K); diag(Pi) <- 0.8
  pi_z <- c(0.5, 0.3, 0.2)
  Z  <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Zh <- vapply(Z, function(z) sample(0:(K - 1), 1, prob = Pi[, z + 1]), 0L)
  x1 <- rnorm(n)
  eta <- c(0, gamma)[Z + 1] + alpha[1] + alpha[2] * x1
  y  <- rpois(n, exp(eta))
  list(y = y, z_hat = Zh, x = cbind(1, x1),
       Pi = Pi, pi_z = pi_z, truth = c(gamma, alpha))
}

# ---- basic correctness ----

test_that("cs_akn is close to truth for binary Poisson", {
  d <- gen_binary(n = 4000, family = "poisson", seed = 42)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = c("naive", "cs", "cs_akn"))
  cf <- coef(fit, method = "cs_akn")
  expect_true(all(abs(cf - d$truth) < 0.07))
  # cs and cs_akn solve different forms of the same equation -> close.
  cf_drift <- coef(fit, method = "cs")
  expect_true(max(abs(cf - cf_drift)) < 0.05)
})

test_that("cs_akn is close to truth for binary binomial", {
  d <- gen_binary(n = 6000, family = "binomial", seed = 11,
                  gamma = 0.9, alpha = c(-0.4, 0.6), p01 = 0.15, p10 = 0.10)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "binomial",
               method = c("naive", "cs", "cs_akn"))
  cf <- coef(fit, method = "cs_akn")
  expect_true(all(abs(cf - d$truth) < 0.10))
})

test_that("cs_akn closed-form Gaussian is consistent with the drift CS", {
  d <- gen_binary(n = 4000, family = "gaussian", seed = 33)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "gaussian",
               method = c("cs", "cs_akn"))
  cf_akn   <- coef(fit, method = "cs_akn")
  cf_drift <- coef(fit, method = "cs")
  # Drift CS and AKN98 surrogate solve different unbiased equations
  # (different finite-sample iterates) but share a population limit.
  # On n = 4000 they agree to a few parts in 1e-4.
  expect_true(max(abs(cf_akn - cf_drift)) < 1e-3)
  expect_true(all(abs(cf_akn - d$truth) < 0.05))
})

test_that("cs_akn works for K=3 Poisson with matrix interface", {
  d <- gen_multi(n = 4000, K = 3, seed = 7)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = c("naive", "cs", "cs_akn"),
               Pi = d$Pi, pi_z = d$pi_z)
  cf <- coef(fit, method = "cs_akn")
  expect_named(cf, c("gamma1", "gamma2", "alpha0", "alpha1"))
  expect_true(all(abs(cf - d$truth) < 0.25))
  # cs_akn should agree closely with the drift cs at this n.
  expect_true(max(abs(cf - coef(fit, method = "cs"))) < 0.05)
})

# ---- variance ----

test_that("cs_akn vcov is positive semidefinite and SE is finite", {
  d <- gen_binary(n = 3000, family = "poisson", seed = 5)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = "cs_akn")
  V <- vcov(fit, method = "cs_akn")
  expect_true(all(is.finite(V)))
  expect_true(all(diag(V) > 0))
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(ev) > -1e-10)
  se <- sqrt(diag(V))
  expect_true(all(is.finite(se)))
  # formula fits carry model.matrix names for the alphas since 0.7.0
  expect_named(se, c("gamma", "(Intercept)", "x1"))
})

test_that("cs_akn matrix-interface accepts only Pi (no pi_z required)", {
  d <- gen_multi(n = 2000, K = 3, seed = 21)
  # cs_akn does not need pi_z; passing only Pi must work. The package
  # always returns the naive fit alongside requested methods.
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = "cs_akn", Pi = d$Pi)
  expect_s3_class(fit, "mcglm")
  expect_true("cs_akn" %in% names(fit$coefficients))
  expect_true(all(is.finite(coef(fit, method = "cs_akn"))))
})

# ---- input validation ----

test_that("cs_akn errors when neither Pi nor (p01,p10) is supplied (binary)", {
  d <- gen_binary(n = 1000, family = "poisson", seed = 1)
  expect_error(
    mcglm(d$df$y, z_hat = as.integer(as.character(d$df$z)),
          x = cbind(1, d$df$x1), family = "poisson", method = "cs_akn"),
    "supply Pi"
  )
})

test_that("cs_akn errors when Pi is missing (multicategory)", {
  d <- gen_multi(n = 1000, K = 3, seed = 2)
  expect_error(
    mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
          method = "cs_akn"),
    "supply Pi"
  )
})

test_that("cs_akn rejects multinomial responses", {
  set.seed(8)
  n <- 500
  Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  z <- rbinom(n, 1, 0.4); z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.1)
  x1 <- rnorm(n)
  y <- sample(0:2, n, replace = TRUE)
  expect_error(
    mcglm(y, z_hat = as.integer(z_hat), x = cbind(1, x1),
          family = "multinomial", method = "cs_akn", Pi = Pi),
    "multinomial"
  )
})

test_that("cs_akn helper reports non-identification clearly", {
  expect_error(
    mismeasured:::.akn_build_Q(diag(2), K = 3),
    "K x K"
  )
  expect_warning(
    expect_error(mismeasured:::.akn_build_Q(matrix(0.5, 2, 2), K = 2),
                 "singular"),
    "near-singular"
  )
})

# ---- weights ----

test_that("cs_akn handles frequency weights", {
  d  <- gen_binary(n = 2000, family = "poisson", seed = 9)
  wt <- runif(nrow(d$df), 0.5, 1.5)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = "cs_akn", weights = wt)
  cf <- coef(fit, method = "cs_akn")
  expect_true(all(is.finite(cf)))
  V <- vcov(fit, method = "cs_akn")
  expect_true(all(is.finite(V)))
})

test_that("cs_akn weighted Gaussian uses the closed-form branch", {
  set.seed(11)
  n <- 200
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  z <- rbinom(n, 1, 0.4)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  x1 <- rnorm(n)
  y <- 0.6 * z - 0.3 + 0.4 * x1 + rnorm(n, sd = 0.5)
  w <- runif(n, 0.5, 1.5)

  fit <- mcglm(y, z_hat = z_hat, x = cbind(1, x1), family = "gaussian",
               method = "cs_akn", Pi = Pi, weights = w)
  expect_true(all(is.finite(coef(fit, method = "cs_akn"))))
})

# ---- agreement with cs as n grows ----

test_that("cs_akn and cs converge to the same point as n grows (Poisson K=2)", {
  d <- gen_binary(n = 20000, family = "poisson", seed = 4242)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = c("cs", "cs_akn"))
  expect_true(max(abs(coef(fit, method = "cs") -
                      coef(fit, method = "cs_akn"))) < 0.01)
})

# ---- Regressed eta and validation tests (Yi et al. 2018 style) ----

test_that("cs_akn works with regressed eta and validation across families", {
  for (fam in c("poisson", "gaussian", "binomial")) {
    for (type in c("internal", "external")) {
      set.seed(101)
      n <- 600
      nv <- 200
      w <- cbind(1, rnorm(n))
      x <- w
      eta_rho <- c(0.2, 0.5)
      eta_01  <- c(-1.4, 0.6)
      eta_10  <- c(-1.3, -0.4)
      rho <- plogis(as.numeric(w %*% eta_rho))
      z   <- rbinom(n, 1, rho)
      p01 <- plogis(as.numeric(w %*% eta_01))
      p10 <- plogis(as.numeric(w %*% eta_10))
      zh  <- ifelse(z == 0L, rbinom(n, 1, p01), 1L - rbinom(n, 1, p10))
      psi_true <- c(0.7, -0.4, 0.5)
      lin_pred <- psi_true[1] * z + psi_true[2] + psi_true[3] * x[, 2]
      y <- switch(fam,
                  poisson  = rpois(n, exp(lin_pred)),
                  binomial = rbinom(n, 1, plogis(lin_pred)),
                  gaussian = rnorm(n, lin_pred, 0.5))
      if (type == "internal") {
        idx <- sample.int(n, nv)
        v <- list(z = z[idx], index = idx)
      } else {
        w_v <- cbind(1, rnorm(nv))
        rho_v <- plogis(as.numeric(w_v %*% eta_rho))
        z_v <- rbinom(nv, 1, rho_v)
        p01_v <- plogis(as.numeric(w_v %*% eta_01))
        p10_v <- plogis(as.numeric(w_v %*% eta_10))
        zh_v <- ifelse(z_v == 0L, rbinom(nv, 1, p01_v), 1L - rbinom(nv, 1, p10_v))
        v <- list(z = z_v, z_hat = zh_v, w = w_v)
      }
      fit <- mcglm(y, z_hat = zh, x = x, family = fam, method = "cs_akn", validation = v)
      expect_s3_class(fit, "mcglm")
      cf <- coef(fit, method = "cs_akn")
      expect_equal(length(cf), 3L)
      expect_true(all(is.finite(cf)))
      V <- vcov(fit, method = "cs_akn")
      expect_true(all(is.finite(V)))
      expect_true(all(diag(V) > 0))
    }
  }
})

test_that("cs_akn internal validation variance matches stacked estimating equation sandwich", {
  set.seed(42)
  n <- 800
  nv <- 300
  w <- cbind(1, rnorm(n))
  x <- w
  q <- ncol(w)
  eta_rho <- c(0.2, 0.6)
  eta_01  <- c(-1.5, 0.8)
  eta_10  <- c(-1.2, -0.6)
  rho <- plogis(as.numeric(w %*% eta_rho))
  z   <- rbinom(n, 1, rho)
  p01 <- plogis(as.numeric(w %*% eta_01))
  p10 <- plogis(as.numeric(w %*% eta_10))
  zh  <- ifelse(z == 0L, rbinom(n, 1, p01), 1L - rbinom(n, 1, p10))
  psi_true <- c(0.8, -0.5, 0.7)
  lin_pred <- psi_true[1] * z + psi_true[2] + psi_true[3] * x[, 2]
  y <- rpois(n, exp(lin_pred))
  v_idx <- sample.int(n, nv)
  v <- list(z = z[v_idx], index = v_idx)
  fit <- mcglm(y, z_hat = zh, x = x, family = "poisson", method = "cs_akn", validation = v)
  V_package <- vcov(fit, method = "cs_akn")

  # Compute stacked sandwich manually
  psi_hat <- coef(fit, method = "cs_akn")
  p01_m <- fit$misclass$p01
  p10_m <- fit$misclass$p10
  x_akn_m <- as.matrix((zh - p01_m) / (1 - p01_m - p10_m))
  Phi <- mismeasured:::.akn_phi(psi_hat, y, x_akn_m, x, 2L, poisson())
  I_akn <- mismeasured:::.akn_information(psi_hat, y, x_akn_m, x, 2L, poisson())
  Ji <- solve(-I_akn)
  p <- length(psi_hat)
  gamma <- psi_hat[1]
  alpha <- psi_hat[-1]
  alpha_lin <- as.numeric(x %*% alpha)
  mu1 <- exp(gamma + alpha_lin)
  mu0 <- exp(alpha_lin)
  delta <- mu1 - mu0
  Delta_m <- 1 - p01_m - p10_m
  v01_m <- p01_m * (1 - p01_m)
  v10_m <- p10_m * (1 - p10_m)
  D_rho <- matrix(0, p, q)
  weight_01 <- (1 - p10_m - zh) / (Delta_m^2) * v01_m
  D_01 <- rbind(-colMeans((y - mu1) * weight_01 * w),
                crossprod(x * (delta * weight_01), w) / n)
  weight_10 <- (zh - p01_m) / (Delta_m^2) * v10_m
  D_10 <- rbind(colMeans((y - mu1) * weight_10 * w),
                -crossprod(x * (delta * weight_10), w) / n)
  D_mat <- cbind(D_rho, D_01, D_10)
  s_mat <- fit$validation$s
  stacked_ef <- matrix(0, n, p + 3 * q)
  stacked_ef[, 1:p] <- Phi
  stacked_ef[v_idx, (p + 1):(p + 3 * q)] <- (n / nv) * s_mat
  top_jacobian <- cbind(Ji, Ji %*% D_mat)
  V_stacked <- top_jacobian %*% crossprod(stacked_ef) %*% t(top_jacobian) / n^2
  expect_equal(unname(V_package), V_stacked, tolerance = 1e-7)
})

test_that("cs_akn formula and mclm interfaces work with validation", {
  set.seed(55)
  n <- 500
  nv <- 200
  w <- cbind(1, rnorm(n))
  x <- w
  eta_rho <- c(0.2, 0.5)
  eta_01  <- c(-1.4, 0.6)
  eta_10  <- c(-1.3, -0.4)
  rho <- plogis(as.numeric(w %*% eta_rho))
  z   <- rbinom(n, 1, rho)
  p01 <- plogis(as.numeric(w %*% eta_01))
  p10 <- plogis(as.numeric(w %*% eta_10))
  zh  <- ifelse(z == 0L, rbinom(n, 1, p01), 1L - rbinom(n, 1, p10))
  y <- rpois(n, exp(0.7 * z - 0.4 + 0.5 * x[, 2]))
  idx <- sample.int(n, nv)
  v <- list(z = z[idx], index = idx)
  df <- data.frame(y = y, z = factor(zh), x1 = x[, 2])
  fit_form <- mcglm(y ~ mc(z) + x1, data = df, family = "poisson", method = "cs_akn", validation = v)
  fit_mat  <- mcglm(y, z_hat = zh, x = x, family = "poisson", method = "cs_akn", validation = v)
  expect_equal(unname(coef(fit_form, method = "cs_akn")), unname(coef(fit_mat, method = "cs_akn")), tolerance = 1e-7)
  expect_equal(unname(vcov(fit_form, method = "cs_akn")), unname(vcov(fit_mat, method = "cs_akn")), tolerance = 1e-7)
  fit_lm <- mclm(y ~ mc(z) + x1, data = df, method = "cs_akn", validation = v)
  expect_equal(length(coef(fit_lm, method = "cs_akn")), 3L)
})
