# Independent equation-level references for the estimated-probability CS
# sandwich. In particular, these do not call the new variance helpers.
.csv_fixture <- function(K = 2L, type = "internal", family = "poisson",
                         n = 700L, nv = 300L, seed = 723L) {
  set.seed(seed)
  Pi <- if (K == 2L) matrix(c(.85, .15, .1, .9), 2L) else
    matrix(c(.85, .1, .05, .1, .85, .05, .05, .1, .85), 3L)
  prevalence <- if (K == 2L) c(.6, .4) else c(.5, .3, .2)
  labels <- function(n) {
    z <- sample.int(K, n, TRUE, prevalence) - 1L
    proxy <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L, 0L)
    list(z = z, z_hat = proxy)
  }
  d <- labels(n)
  x <- cbind(1, rnorm(n))
  eta <- c(0, .65, -.4)[d$z + 1L] - .4 + .5 * x[, 2L]
  y <- switch(family, poisson = rpois(n, exp(eta)),
              binomial = rbinom(n, 1, plogis(eta)), gaussian = rnorm(n, eta))
  if (type == "internal") {
    index <- sample.int(n, nv)
    v <- list(z = d$z[index], index = index)
    proxy_v <- d$z_hat[index]
  } else {
    v <- labels(nv)
    proxy_v <- v$z_hat
  }
  B <- matrix(tabulate(proxy_v + 1L + K * v$z, nbins = K * K), K) / nv
  pi <- colSums(B)
  args <- list(formula = y, z_hat = d$z_hat, x = x, family = family,
               Pi = sweep(B, 2, pi, "/"), pi_z = if (K == 2L) pi[2L] else pi,
               method = c("naive", "bca", "bcm", "cs"))
  list(args = args, v = v, B = B, proxy_v = proxy_v)
}

.csv_jacobian <- function(f, a, h = 1e-5) {
  vapply(seq_along(a), function(j) {
    plus <- minus <- a
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    (f(plus) - f(minus)) / (2 * h)
  }, numeric(length(f(a))))
}

.csv_reference <- function(fit, d) {
  psi <- unname(coef(fit, method = "cs"))
  K <- fit$K
  q <- fit$x
  n <- fit$n
  nv <- length(d$v$z)
  mu <- fit$family$linkinv
  # The multicategory equation directly in joint probabilities, also at K=2.
  score <- function(beta, b) {
    gamma <- c(0, beta[seq_len(K - 1L)])
    eta <- drop(q %*% beta[-seq_len(K - 1L)])
    mus <- vapply(gamma, function(g) mu(eta + g), numeric(n))
    ef <- fit$xi_hat * (fit$y - mu(drop(fit$xi_hat %*% beta)))
    for (l in seq_len(K)) for (j in seq_len(K)) {
      dummy <- numeric(K - 1L)
      if (j > 1L) dummy[j - 1L] <- 1
      xi <- cbind(matrix(rep(dummy, each = n), n), q)
      ef <- ef + b[j + K * (l - 1L)] * xi * (mus[, j] - mus[, l])
    }
    ef
  }
  b <- as.numeric(d$B)
  phi <- score(psi, b)
  J <- .csv_jacobian(function(beta) colMeans(score(beta, b)), psi)
  D <- .csv_jacobian(function(b) colMeans(score(psi, b)), b)
  cell <- d$proxy_v + 1L + K * d$v$z
  s <- diag(K * K)[cell, , drop = FALSE] - rep(b, each = nv)
  Sigma <- crossprod(s) / nv
  meat <- crossprod(phi) / n + n / nv * D %*% Sigma %*% t(D)
  if (!is.null(d$v$index)) {
    pv <- phi[d$v$index, , drop = FALSE]
    Gamma <- crossprod(sweep(pv, 2, colMeans(pv)), s) / nv
    meat <- meat + Gamma %*% t(D) + D %*% t(Gamma)
  }
  Ji <- solve(J)
  list(V = Ji %*% meat %*% t(Ji) / n,
       old = Ji %*% (crossprod(phi) / n) %*% t(Ji) / n,
       phi = phi, J = J, D = D, s = s)
}

test_that("validation sandwich matches independent derivatives for all GLM families", {
  for (K in c(2L, 3L)) for (family in c("poisson", "binomial", "gaussian")) {
    for (type in c("external", "internal")) {
      d <- .csv_fixture(K, type, family)
      old <- do.call(mcglm, d$args)
      fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
      ref <- .csv_reference(fit, d)
      expect_equal(unname(vcov(fit, method = "cs")), ref$V, tolerance = 1e-7)
      expect_equal(unname(vcov(old, method = "cs")), ref$old, tolerance = 1e-7)
      expect_identical(fit$coefficients, old$coefficients)
      for (method in c("naive", "bca", "bcm"))
        expect_identical(vcov(fit, method = method), vcov(old, method = method))
    }
  }
})

test_that("binary eta influence formula equals joint-cell implementation", {
  d <- .csv_fixture()
  fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
  psi <- unname(coef(fit, method = "cs"))
  pi <- d$args$pi_z
  p01 <- d$args$Pi[2, 1]
  p10 <- d$args$Pi[1, 2]
  z <- d$v$z
  zh <- d$proxy_v
  nv <- length(z)
  # Equations s-binary-validation and D-binary-explicit, verbatim algebra.
  s <- cbind(z - pi, (z == 0) * (zh - p01) / (1 - pi),
             (z == 1) * ((1 - zh) - p10) / pi)
  eta <- drop(fit$x %*% psi[-1])
  delta <- fit$family$linkinv(eta + psi[1]) - fit$family$linkinv(eta)
  dq <- colMeans(delta * fit$x)
  D <- rbind(c(-p01 * mean(delta), (1 - pi) * mean(delta), 0),
             cbind(-(p01 + p10) * dq, (1 - pi) * dq, -pi * dq))
  ref <- .csv_reference(fit, d)
  pv <- ref$phi[d$v$index, , drop = FALSE]
  Gamma <- crossprod(sweep(pv, 2, colMeans(pv)), s) / nv
  meat <- crossprod(ref$phi) / fit$n + fit$n / nv * D %*% (crossprod(s) / nv) %*% t(D) +
    Gamma %*% t(D) + D %*% t(Gamma)
  Ji <- solve(ref$J)
  expect_equal(unname(vcov(fit, method = "cs")), Ji %*% meat %*% t(Ji) / fit$n,
               tolerance = 1e-7)
})

test_that("internal covariance agrees with a stacked estimating-equation sandwich", {
  for (nv in c(300L, 700L)) {
    d <- .csv_fixture(nv = nv)
    fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
    ref <- .csv_reference(fit, d)
    nuisance <- matrix(0, fit$n, ncol(ref$s))
    nuisance[d$v$index, ] <- fit$n / nv * ref$s
    ef <- cbind(ref$phi, nuisance)
    # Upper block of inverse stacked Jacobian [[J,D],[0,-I]].
    Ji <- solve(ref$J)
    top <- cbind(Ji, Ji %*% ref$D)
    stacked <- top %*% crossprod(ef) %*% t(top) / fit$n^2
    expect_equal(unname(vcov(fit, method = "cs")), stacked, tolerance = 1e-7)
    ext <- list(z = d$v$z, z_hat = d$proxy_v)
    independent <- do.call(mcglm, c(d$args, list(validation = ext)))
    expect_gt(max(abs(vcov(fit, method = "cs") - vcov(independent, method = "cs"))), 1e-5)
  }
})

test_that("validation SE propagates to summary and confidence intervals", {
  d <- .csv_fixture()
  fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
  se <- sqrt(diag(vcov(fit, method = "cs")))
  expect_equal(fit$se$cs, se)
  ci <- confint(fit, method = "cs")
  expect_equal(unname(ci[, 2] - ci[, 1]), unname(2 * qnorm(.975) * se))
  expect_equal(unname(summary(fit)$coefficients$cs[, "Std. Error"]), unname(se))
})

test_that("validation metadata rejects inconsistent probabilities and unsupported designs", {
  d <- .csv_fixture()
  fit <- function(v = d$v, args = d$args) do.call(mcglm, c(args, list(validation = v)))
  expect_error(fit(list(z = d$v$z, index = rep(1, length(d$v$z)))), "distinct")
  expect_error(fit(list(z = d$v$z, index = d$v$index + 1000L)), "row numbers")
  expect_error(fit(list(z = d$v$z, index = d$v$index, z_hat = 1 - d$proxy_v)), "does not match")
  expect_error(fit(list(z = factor(d$v$z), index = d$v$index)), "codes")
  expect_error(fit(list(z = c(NA, d$v$z[-1]), index = d$v$index)), "codes")
  expect_error(fit(list(z = rep(0, length(d$v$z)), index = d$v$index)), "Every true category")
  expect_error(fit(list(z = d$v$z)), "z_hat")
  wrong <- d$args
  wrong$pi_z <- .8
  expect_error(fit(args = wrong), "empirical estimates")
  wrong <- d$args
  wrong$weights <- rep(2, length(wrong$formula))
  expect_error(fit(args = wrong), "unweighted")
  wrong$weights[] <- 1
  expect_equal(vcov(fit(args = wrong), method = "cs"), vcov(fit(), method = "cs"))
  wrong <- d$args
  wrong$method <- "bcm"
  expect_error(fit(args = wrong), "requires method")
})

test_that("formula, binary c1/c2, and multicategory numerical Jacobian interfaces work", {
  d <- .csv_fixture()
  a <- d$args
  df <- data.frame(y = a$formula, z = a$z_hat, x = a$x[, 2])
  fitted <- mcglm(y ~ mc(z) + x, data = df, family = a$family, method = "cs",
                  Pi = a$Pi, pi_z = a$pi_z, validation = d$v)
  direct <- do.call(mcglm, c(a, list(validation = d$v)))
  expect_equal(unname(vcov(fitted, method = "cs")), unname(vcov(direct, method = "cs")))
  a$Pi <- a$pi_z <- NULL
  a$c1 <- d$B[2, 1]
  a$c2 <- d$B[2, 1] - d$B[1, 2]
  constants <- do.call(mcglm, c(a, list(validation = d$v)))
  expect_equal(vcov(constants, method = "cs"), vcov(direct, method = "cs"))
  d3 <- .csv_fixture(3L)
  analytical <- do.call(mcglm, c(d3$args, list(validation = d3$v)))
  numerical <- do.call(mcglm, c(d3$args, list(validation = d3$v, jacobian = "numerical")))
  expect_equal(vcov(analytical, method = "cs"), vcov(numerical, method = "cs"), tolerance = 1e-5)
})

test_that("external validation uses 1/n_validation scaling even when n_validation > n", {
  d <- .csv_fixture(type = "external", nv = 900L)
  old <- do.call(mcglm, d$args)
  fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
  doubled <- lapply(d$v, rep, times = 2L)
  fit2 <- do.call(mcglm, c(d$args, list(validation = doubled)))
  base <- vcov(old, method = "cs")
  expect_equal(vcov(fit2, method = "cs") - base,
               (vcov(fit, method = "cs") - base) / 2, tolerance = 1e-12)
})

test_that("zero empirical misclassification adds no validation uncertainty", {
  for (K in c(2L, 3L)) {
    d <- .csv_fixture(K, type = "external")
    d$v$z_hat <- d$v$z
    d$args$Pi <- diag(K)
    prevalence <- tabulate(d$v$z + 1L, nbins = K) / length(d$v$z)
    d$args$pi_z <- if (K == 2L) prevalence[2] else prevalence
    old <- do.call(mcglm, d$args)
    fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
    expect_equal(vcov(fit, method = "cs"), vcov(old, method = "cs"), tolerance = 1e-12)
  }
})

# ---- Tests for regressed parameters eta ----

.cs_regressed_fixture <- function(family = "poisson", type = "internal",
                                  n = 600L, nv = 250L, seed = 42L) {
  set.seed(seed)
  w <- cbind(1, rnorm(n))
  x <- w
  q <- ncol(w)
  eta_rho <- c(0.3, -0.5)
  eta_01  <- c(-1.4, 0.4)
  eta_10  <- c(-1.3, -0.3)
  rho <- plogis(as.numeric(w %*% eta_rho))
  z <- rbinom(n, 1, rho)
  p01 <- plogis(as.numeric(w %*% eta_01))
  p10 <- plogis(as.numeric(w %*% eta_10))
  zh <- ifelse(z == 0L, rbinom(n, 1, p01), 1L - rbinom(n, 1, p10))
  psi_true <- c(0.7, -0.4, 0.6)
  lin_pred <- psi_true[1] * z + as.numeric(x %*% psi_true[-1])
  y <- switch(family,
              poisson  = rpois(n, exp(lin_pred)),
              binomial = rbinom(n, 1, plogis(lin_pred)),
              gaussian = rnorm(n, lin_pred, 0.5))
  if (type == "internal") {
    index <- sample.int(n, nv)
    v <- list(z = z[index], index = index)
    proxy_v <- zh[index]
    w_v <- w[index, , drop = FALSE]
  } else {
    w_v <- cbind(1, rnorm(nv))
    rho_v <- plogis(as.numeric(w_v %*% eta_rho))
    z_v <- rbinom(nv, 1, rho_v)
    p01_v <- plogis(as.numeric(w_v %*% eta_01))
    p10_v <- plogis(as.numeric(w_v %*% eta_10))
    proxy_v <- ifelse(z_v == 0L, rbinom(nv, 1, p01_v), 1L - rbinom(nv, 1, p10_v))
    v <- list(z = z_v, z_hat = proxy_v, w = w_v)
  }
  list(y = y, zh = zh, x = x, w = w, v = v, family = family,
       type = type, n = n, nv = nv, w_v = w_v, proxy_v = proxy_v,
       eta_true = list(rho = eta_rho, p01 = eta_01, p10 = eta_10))
}

.cs_regressed_reference <- function(fit, d) {
  psi <- unname(coef(fit, method = "cs"))
  p <- length(psi)
  q <- ncol(d$w)
  n <- d$n
  nv <- d$nv
  fam <- fit$family
  mu_fun <- fam$linkinv
  mu_dot_fun <- fam$mu.eta
  v <- d$v
  z_v <- v$z
  zh_v <- d$proxy_v
  w_v <- d$w_v
  fit_rho <- glm.fit(w_v, z_v, family = binomial())
  fit_01  <- glm.fit(w_v[z_v == 0, , drop = FALSE], zh_v[z_v == 0], family = binomial())
  fit_10  <- glm.fit(w_v[z_v == 1, , drop = FALSE], 1 - zh_v[z_v == 1], family = binomial())
  eta_hat <- list(rho = unname(coef(fit_rho)),
                  p01 = unname(coef(fit_01)),
                  p10 = unname(coef(fit_10)))
  rho_v <- plogis(as.numeric(w_v %*% eta_hat$rho))
  p01_v <- plogis(as.numeric(w_v %*% eta_hat$p01))
  p10_v <- plogis(as.numeric(w_v %*% eta_hat$p10))
  v_rho_v <- rho_v * (1 - rho_v)
  v_01_v  <- p01_v * (1 - p01_v)
  v_10_v  <- p10_v * (1 - p10_v)
  H_rho <- crossprod(w_v * v_rho_v, w_v) / nv
  H_01  <- crossprod(w_v * ((1 - z_v) * v_01_v), w_v) / nv
  H_10  <- crossprod(w_v * (z_v * v_10_v), w_v) / nv
  s_rho <- (w_v * (z_v - rho_v)) %*% solve(H_rho)
  s_01  <- (w_v * ((1 - z_v) * (zh_v - p01_v))) %*% solve(H_01)
  s_10  <- (w_v * (z_v * (1 - zh_v - p10_v))) %*% solve(H_10)
  s_mat <- cbind(s_rho, s_01, s_10)
  rho_m <- plogis(as.numeric(d$w %*% eta_hat$rho))
  p01_m <- plogis(as.numeric(d$w %*% eta_hat$p01))
  p10_m <- plogis(as.numeric(d$w %*% eta_hat$p10))
  c1_m  <- p01_m * (1 - rho_m)
  c2_m  <- c1_m - p10_m * rho_m
  xi_hat <- cbind(d$zh, d$x)
  resid <- d$y - mu_fun(as.numeric(xi_hat %*% psi))
  gamma <- psi[1]
  alpha <- psi[-1]
  eta0 <- as.numeric(d$x %*% alpha)
  delta <- mu_fun(gamma + eta0) - mu_fun(eta0)
  m_mat <- matrix(0, n, p)
  m_mat[, 1]  <- -c1_m * delta
  m_mat[, -1] <- -c2_m * delta * d$x
  phi <- xi_hat * resid - m_mat
  I_hat <- crossprod(xi_hat * mu_dot_fun(as.numeric(xi_hat %*% psi)), xi_hat) / n
  mu_dot1 <- mu_dot_fun(gamma + eta0)
  mu_dot0 <- mu_dot_fun(eta0)
  d_mu_dot <- mu_dot1 - mu_dot0
  M_hat <- matrix(0, p, p)
  M_hat[1, 1]    <- -mean(c1_m * mu_dot1)
  M_hat[1, 2:p]  <- -colMeans(c1_m * d_mu_dot * d$x)
  M_hat[2:p, 1]  <- -colMeans(c2_m * mu_dot1 * d$x)
  M_hat[2:p, 2:p]<- -crossprod(c2_m * d_mu_dot * d$x, d$x) / n
  J <- -(I_hat + M_hat)
  Ji <- solve(J)
  v_rho_m <- rho_m * (1 - rho_m)
  v_01_m  <- p01_m * (1 - p01_m)
  v_10_m  <- p10_m * (1 - p10_m)
  D_rho <- matrix(0, p, q)
  D_rho[1, ]  <- -colMeans(delta * p01_m * v_rho_m * d$w)
  D_rho[-1, ] <- -crossprod(d$x * (delta * (p01_m + p10_m) * v_rho_m), d$w) / n
  D_01 <- matrix(0, p, q)
  D_01[1, ]  <- colMeans(delta * (1 - rho_m) * v_01_m * d$w)
  D_01[-1, ] <- crossprod(d$x * (delta * (1 - rho_m) * v_01_m), d$w) / n
  D_10 <- matrix(0, p, q)
  D_10[1, ]  <- 0
  D_10[-1, ] <- -crossprod(d$x * (delta * rho_m * v_10_m), d$w) / n
  D_mat <- cbind(D_rho, D_01, D_10)
  U <- s_mat %*% t(D_mat)
  S <- crossprod(phi) / n
  meat <- S + (n / nv) * crossprod(U) / nv
  if (d$type == "internal") {
    phi_v <- phi[d$v$index, , drop = FALSE]
    phi_v <- sweep(phi_v, 2, colMeans(phi_v))
    cross <- crossprod(phi_v, U) / nv
    meat <- meat + cross + t(cross)
  }
  V <- Ji %*% meat %*% t(Ji) / n
  list(V = V, J = J, D = D_mat, s = s_mat, phi = phi)
}

test_that("validation sandwich with regressed eta matches independent references for all families", {
  for (type in c("internal", "external")) {
    for (fam in c("poisson", "gaussian", "binomial")) {
      d <- .cs_regressed_fixture(family = fam, type = type)
      fit <- mcglm(d$y, z_hat = d$zh, x = d$x, family = d$family, method = "cs", validation = d$v)
      ref <- .cs_regressed_reference(fit, d)
      expect_equal(unname(vcov(fit, method = "cs")), ref$V, tolerance = 1e-7)
    }
  }
})

test_that("internal covariance with regressed eta agrees with stacked estimating equation sandwich", {
  d <- .cs_regressed_fixture(family = "poisson", type = "internal", nv = 300L)
  fit <- mcglm(d$y, z_hat = d$zh, x = d$x, family = d$family, method = "cs", validation = d$v)
  ref <- .cs_regressed_reference(fit, d)
  p <- length(coef(fit, method = "cs"))
  q <- ncol(d$w)
  stacked_ef <- matrix(0, d$n, p + 3 * q)
  stacked_ef[, 1:p] <- ref$phi
  stacked_ef[d$v$index, (p + 1):(p + 3 * q)] <- (d$n / d$nv) * ref$s
  Ji <- solve(ref$J)
  top_jacobian <- cbind(Ji, Ji %*% ref$D)
  stacked_V <- top_jacobian %*% crossprod(stacked_ef) %*% t(top_jacobian) / d$n^2
  expect_equal(unname(vcov(fit, method = "cs")), stacked_V, tolerance = 1e-7)
})

test_that("regressed eta with w = 1 matches constant probability formula", {
  d <- .cs_regressed_fixture(family = "poisson", type = "internal")
  w_const <- matrix(1, d$n, 1)
  v_const <- list(z = d$v$z, index = d$v$index, w = matrix(1, d$nv, 1))
  fit_const_reg <- mcglm(d$y, z_hat = d$zh, x = cbind(1, d$x[, 2]), family = "poisson", method = "cs",
                         w = w_const, validation = v_const)
  B <- matrix(tabulate(d$zh[d$v$index] + 1L + 2L * d$v$z, nbins = 4), 2) / d$nv
  c1_const <- B[2, 1]
  c2_const <- B[2, 1] - B[1, 2]
  v_old <- list(z = d$v$z, index = d$v$index)
  fit_const_old <- mcglm(d$y, z_hat = d$zh, x = cbind(1, d$x[, 2]), family = "poisson", method = "cs",
                         c1 = c1_const, c2 = c2_const, validation = v_old)
  expect_equal(vcov(fit_const_reg, method = "cs"), vcov(fit_const_old, method = "cs"), tolerance = 1e-8)
})

test_that("eta list and vector inputs give identical point estimates and variances", {
  d <- .cs_regressed_fixture()
  eta_list <- d$eta_true
  eta_vec <- c(eta_list$rho, eta_list$p01, eta_list$p10)
  fit_list <- mcglm(d$y, z_hat = d$zh, x = d$x, family = d$family, method = "cs", eta = eta_list)
  fit_vec  <- mcglm(d$y, z_hat = d$zh, x = d$x, family = d$family, method = "cs", eta = eta_vec)
  expect_identical(coef(fit_list, method = "cs"), coef(fit_vec, method = "cs"))
  expect_identical(vcov(fit_list, method = "cs"), vcov(fit_vec, method = "cs"))
})

test_that("formula and mclm interfaces work with regressed eta validation", {
  d <- .cs_regressed_fixture()
  df <- data.frame(y = d$y, z = factor(d$zh), x1 = d$x[, 2])
  fit_formula <- mcglm(y ~ mc(z) + x1, data = df, family = d$family, method = "cs", validation = d$v)
  fit_matrix  <- mcglm(d$y, z_hat = d$zh, x = d$x, family = d$family, method = "cs", validation = d$v)
  expect_equal(unname(coef(fit_formula, method = "cs")), unname(coef(fit_matrix, method = "cs")))
  expect_equal(unname(vcov(fit_formula, method = "cs")), unname(vcov(fit_matrix, method = "cs")))
  fit_lm <- mclm(y ~ mc(z) + x1, data = df, method = "cs", validation = d$v)
  expect_equal(length(coef(fit_lm, method = "cs")), 3L)
})

test_that("regressed eta validation rejects incompatible inputs", {
  d <- .cs_regressed_fixture()
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs", eta = c(1, 2)), "ncol")
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs", eta = list(rho = 1)), "contain 'rho', 'p01', and 'p10'")
  bad_v <- d$v
  bad_v$w <- matrix(1, d$nv, 5) # mismatch with w in regression
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs", validation = bad_v), "cols")
})

test_that("mcglm handles w supplied as formula, character, and data.frame", {
  d <- .cs_regressed_fixture(family = "poisson", n = 400L, nv = 150L)
  df <- data.frame(y = d$y, z = factor(d$zh), x1 = d$x[, 2])

  # Formula interface with w as formula
  fit_w_form <- mcglm(y ~ mc(z) + x1, data = df, family = "poisson", method = "cs",
                      w = ~ x1, validation = d$v)
  expect_true(all(is.finite(coef(fit_w_form, method = "cs"))))

  # Formula interface with w as character
  fit_w_char <- mcglm(y ~ mc(z) + x1, data = df, family = "poisson", method = "cs",
                      w = "x1", validation = d$v)
  expect_equal(coef(fit_w_form, method = "cs"), coef(fit_w_char, method = "cs"))

  # Formula interface with w as data.frame (including intercept)
  w_df <- data.frame(1, x1 = df$x1)
  fit_w_df <- mcglm(y ~ mc(z) + x1, data = df, family = "poisson", method = "cs",
                    w = w_df, validation = d$v)
  expect_equal(coef(fit_w_form, method = "cs"), coef(fit_w_df, method = "cs"))

  # mclm with w as formula
  fit_mclm_w <- mclm(y ~ mc(z) + x1, data = df, method = "cs", w = ~ x1, validation = d$v)
  expect_true(all(is.finite(coef(fit_mclm_w, method = "cs"))))

  # mclm with cs_akn and validation
  fit_mclm_akn <- mclm(y ~ mc(z) + x1, data = df, method = "cs_akn", w = ~ x1, validation = d$v)
  expect_true(all(is.finite(coef(fit_mclm_akn, method = "cs_akn"))))

  # mclm with eta
  fit_mclm_eta <- mclm(y ~ mc(z) + x1, data = df, method = "cs", eta = d$eta_true)
  expect_true(all(is.finite(coef(fit_mclm_eta, method = "cs"))))
})

test_that("c1 and c2 validation checks catch invalid inputs", {
  d <- .cs_regressed_fixture(n = 100L, nv = 50L)
  df <- data.frame(y = d$y, z = factor(d$zh), x1 = d$x[, 2])

  # In formula interface
  expect_error(mcglm(y ~ mc(z) + x1, data = df, c1 = "bad"), "c1 must be one finite numeric value")
  expect_error(mcglm(y ~ mc(z) + x1, data = df, c1 = c(1, 2)), "c1 must be one finite numeric value")
  expect_error(mcglm(y ~ mc(z) + x1, data = df, c1 = Inf), "c1 must be one finite numeric value")
  expect_error(mcglm(y ~ mc(z) + x1, data = df, c2 = "bad"), "c2 must be one finite numeric value")
  expect_error(mcglm(y ~ mc(z) + x1, data = df, c2 = c(1, 2)), "c2 must be one finite numeric value")
  expect_error(mcglm(y ~ mc(z) + x1, data = df, c2 = NaN), "c2 must be one finite numeric value")

  # In matrix interface (.mcglm_fit)
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, c1 = "bad"), "c1 must be one finite numeric value")
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, c1 = c(1, 2)), "c1 must be one finite numeric value")
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, c1 = Inf), "c1 must be one finite numeric value")
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, c2 = "bad"), "c2 must be one finite numeric value")
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, c2 = c(1, 2)), "c2 must be one finite numeric value")
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, c2 = NaN), "c2 must be one finite numeric value")
})

test_that("validation requires cs or cs_akn method", {
  d <- .cs_regressed_fixture(n = 100L, nv = 50L)
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, method = "naive", validation = d$v),
               "validation requires method = 'cs' or 'cs_akn'")
  expect_error(mcglm(d$y, z_hat = d$zh, x = d$x, method = "bca", validation = d$v),
               "validation requires method = 'cs' or 'cs_akn'")
})

test_that("cs_akn clamping warns when Delta is near zero", {
  d <- .cs_regressed_fixture(n = 100L, nv = 50L)
  expect_warning(
    mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs_akn",
          p01 = 0.5, p10 = 0.5),
    "clamping applied"
  )
  expect_warning(
    mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs_akn",
          p01 = 0.5002, p10 = 0.5002),
    "clamping applied"
  )
})

test_that("cs_akn solver failure returns NA and warns", {
  expect_warning(
    res <- mismeasured:::.mcglm_fit_cs_akn(
      c(30, -30), c(rep(100, 5), rep(0, 5)),
      matrix(c(rep(-10, 5), rep(10, 5)), 10, 1),
      matrix(c(rep(1, 5), rep(-1, 5)), 10, 1),
      K = 2L, family = "poisson"
    ),
    "did not converge"
  )
  expect_true(all(is.na(res)))
})

test_that("unit weights work with cs_akn validation and cs regressed eta", {
  d <- .cs_regressed_fixture(family = "poisson", n = 400L, nv = 150L)
  unit_wt <- rep(1, d$n)

  # cs_akn with validation and weights
  fit_akn_wt <- mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs_akn",
                      validation = d$v, weights = unit_wt)
  expect_true(all(is.finite(coef(fit_akn_wt, method = "cs_akn"))))
  expect_true(all(is.finite(vcov(fit_akn_wt, method = "cs_akn"))))

  # cs with regressed eta and weights
  fit_cs_wt <- mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs",
                     eta = d$eta_true, weights = unit_wt)
  expect_true(all(is.finite(coef(fit_cs_wt, method = "cs"))))
  expect_true(all(is.finite(vcov(fit_cs_wt, method = "cs"))))
})

test_that(".mcglm_normalize_eta handles glm, lm, and error branches", {
  d <- .cs_regressed_fixture(n = 200L, nv = 100L)
  df_val <- data.frame(z = d$v$z, zh = d$zh[d$v$index], x1 = d$x[d$v$index, 2])

  # Fitted glm objects directly into .mcglm_normalize_eta
  g_rho <- stats::glm(z ~ x1, data = df_val, family = stats::binomial)
  g_01  <- stats::glm(zh ~ x1, data = df_val[df_val$z == 0, ], family = stats::binomial)
  g_10  <- stats::glm(I(1 - zh) ~ x1, data = df_val[df_val$z == 1, ], family = stats::binomial)
  norm_glm <- mismeasured:::.mcglm_normalize_eta(list(rho = g_rho, p01 = g_01, p10 = g_10), 2)
  expect_equal(norm_glm$rho, unname(coef(g_rho)))
  expect_equal(norm_glm$p01, unname(coef(g_01)))
  expect_equal(norm_glm$p10, unname(coef(g_10)))

  # Fitted lm objects
  l_rho <- stats::lm(z ~ x1, data = df_val)
  norm_lm <- mismeasured:::.mcglm_normalize_eta(list(rho = l_rho, p01 = g_01, p10 = g_10), 2)
  expect_equal(norm_lm$rho, unname(coef(l_rho)))

  # mcglm with eta list
  fit_eta_norm <- mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs", eta = d$eta_true)
  expect_true(all(is.finite(coef(fit_eta_norm, method = "cs"))))

  # Error: non-numeric or non-finite
  expect_error(mismeasured:::.mcglm_normalize_eta(list(rho = "a", p01 = c(1,1), p10 = c(1,1)), 2),
               "must be numeric")
  expect_error(mismeasured:::.mcglm_normalize_eta(list(rho = c(NaN, 1), p01 = c(1,1), p10 = c(1,1)), 2),
               "must be numeric")

  # Error: wrong length in list
  expect_error(mismeasured:::.mcglm_normalize_eta(list(rho = 1, p01 = c(1,1), p10 = c(1,1)), 2),
               "has length 1, but w has 2")

  # Error: missing elements in list
  expect_error(mismeasured:::.mcglm_normalize_eta(list(rho = c(1,1)), 2),
               "must contain 'rho', 'p01', and 'p10'")

  # Numeric vector of length 3 * q
  norm_vec <- mismeasured:::.mcglm_normalize_eta(rep(0.2, 6), 2)
  expect_equal(length(norm_vec$rho), 2L)
  expect_equal(length(norm_vec$p01), 2L)
  expect_equal(length(norm_vec$p10), 2L)

  # Error: numeric vector wrong length or non-finite
  expect_error(mismeasured:::.mcglm_normalize_eta(c(1, 2, 3), 2), "must have length 3 \\* ncol")
  expect_error(mismeasured:::.mcglm_normalize_eta(c(1, 2, 3, 4, 5, Inf), 2), "must have length 3 \\* ncol")

  # Error: non-numeric non-list
  expect_error(mismeasured:::.mcglm_normalize_eta("bad", 2), "must be a named list")
})

test_that(".mcglm_prepare_cs_validation branches and errors are covered", {
  d <- .cs_regressed_fixture(n = 200L, nv = 100L)

  # Passing eta WITH validation list
  fit_eta_val <- mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs",
                       eta = d$eta_true, validation = d$v)
  expect_true(all(is.finite(coef(fit_eta_val, method = "cs"))))

  # validation$x supplied instead of validation$w
  v_x <- d$v
  v_x$x <- d$w[d$v$index, , drop = FALSE]
  v_x$w <- NULL
  fit_vx <- mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs", validation = v_x)
  expect_true(all(is.finite(coef(fit_vx, method = "cs"))))

  # w specified in mcglm with internal validation
  fit_w_int <- mcglm(d$y, z_hat = d$zh, x = d$x, w = d$w, method = "cs", validation = d$v)
  expect_true(all(is.finite(coef(fit_w_int, method = "cs"))))

  # External validation without w when x has multiple columns
  expect_error(
    mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs",
          validation = list(z = d$v$z, z_hat = d$zh[1:d$nv])),
    "supply w in validation list"
  )

  # External validation without w when x has 1 column
  fit_ext_1col <- mcglm(d$y, z_hat = d$zh, x = matrix(1, d$n, 1), method = "cs",
                        validation = list(z = d$v$z, z_hat = d$zh[1:d$nv]))
  expect_true(all(is.finite(coef(fit_ext_1col, method = "cs"))))

  # validation$z missing category 0 or 1
  expect_error(
    mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs",
          validation = list(z = rep(1L, d$nv), index = d$v$index)),
    "Every true category must occur"
  )
  expect_error(
    mcglm(d$y, z_hat = d$zh, x = d$x, method = "cs",
          validation = list(z = rep(0L, d$nv), index = d$v$index)),
    "Every true category must occur"
  )

  # Ginv fallback in .mcglm_vcov_cs_akn
  V_sing <- mismeasured:::.mcglm_vcov_cs_akn(
    c(0, 0), rep(1, 10), matrix(0, 10, 1), matrix(0, 10, 1),
    K = 2L, family = "gaussian"
  )
  expect_true(is.matrix(V_sing))
})
