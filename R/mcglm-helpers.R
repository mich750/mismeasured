# ---------------------------------------------------------------------------
# Internal helper functions for computing drift, Fisher info, Jacobian
# Used by mcglm() estimators
# ---------------------------------------------------------------------------

#' Estimate pi_z (true prevalence) from observed proxy and Pi
#'
#' Uses Bayesian inversion: if pi_obs = Pi %*% pi_z, then
#' pi_z = solve(Pi) %*% pi_obs. Clamps to [0.01, 0.99] for stability.
#' @param z_hat Integer vector of observed proxy values (0-based).
#' @param Pi K x K misclassification matrix.
#' @return Numeric vector of length K (estimated true prevalences).
#' @keywords internal
.mcglm_estimate_pi_z <- function(z_hat, Pi) {
  K <- nrow(Pi)
  # Coding errors must surface as such, not as empty-category errors
  if (any(z_hat < 0L) || any(z_hat >= K))
    stop("z_hat must be coded 0, ..., K-1 (K = ", K, "); got values in [",
         min(z_hat), ", ", max(z_hat), "].", call. = FALSE)
  # Observed proportions
  tab <- tabulate(z_hat + 1L, nbins = K)
  if (any(tab == 0L))
    stop("Cannot estimate pi_z: proxy category ",
         paste(which(tab == 0L) - 1L, collapse = ", "),
         " has no observations in z_hat (K = ", K, " from Pi). ",
         "Clamping would fabricate a prevalence for it. Supply pi_z ",
         "directly, or collapse Pi to the observed categories.",
         call. = FALSE)
  pi_obs <- tab / sum(tab)
  # Invert: pi_z = Pi^{-1} %*% pi_obs
  pi_raw <- tryCatch(
    as.numeric(solve(Pi, pi_obs)),
    error = function(e)
      stop("Cannot estimate pi_z: the misclassification matrix Pi is ",
           "singular. Supply pi_z directly.", call. = FALSE)
  )
  # Clamp to valid range (two-sided) and renormalize
  pi_z <- pmin(pmax(pi_raw, 0.01), 0.99)
  pi_z <- pi_z / sum(pi_z)
  if (max(abs(pi_z - pi_raw)) > 0.01)
    warning("Estimated pi_z required clamping to [0.01, 0.99]: ",
            "Pi-inversion gave (", paste(round(pi_raw, 4), collapse = ", "),
            "), using (", paste(round(pi_z, 4), collapse = ", "), "). ",
            "The observed proxy frequencies are barely consistent with ",
            "Pi; consider supplying pi_z directly.", call. = FALSE)
  if (K == 2L) return(pi_z[2])  # scalar for binary case
  pi_z
}

#' Build xi_hat design matrix for mcglm
#'
#' For binary Z: xi_hat = cbind(z_hat, x).
#' For multicategory Z: xi_hat = cbind(d_hat, x) where d_hat is the
#' dummy-encoding of z_hat with baseline = 0.
#' @param z_hat Integer proxy covariate vector.
#' @param x Covariate matrix (n x r).
#' @param K Number of Z categories.
#' @return n x p matrix.
#' @keywords internal
.mcglm_build_xi_hat <- function(z_hat, x, K) {
  if (K == 2L) return(cbind(z_hat, x))
  n <- length(z_hat)
  s <- K - 1L
  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  cbind(d_hat, x)
}

# ---- Binary misclassification helpers ----

#' Compute the toggle gap delta_i(psi) for binary misclassification
#' @keywords internal
.mcglm_compute_delta <- function(psi, x, mu_fun) {
  gamma <- psi[1]
  alpha <- psi[-1]
  eta0 <- as.numeric(x %*% alpha)
  mu_fun(gamma + eta0) - mu_fun(eta0)
}

#' Compute per-observation drift m_hat(psi) for binary misclassification
#' @keywords internal
.mcglm_compute_m_bin <- function(psi, x, mu_fun, c1, c2) {
  n <- nrow(x)
  p <- length(psi)
  delta <- .mcglm_compute_delta(psi, x, mu_fun)
  m <- matrix(0, n, p)
  m[, 1]  <- -c1 * delta
  m[, -1] <- -c2 * delta * x
  m
}

#' Compute m_hat(psi) = weighted mean of m_i
#' @keywords internal
.mcglm_compute_mhat_bin <- function(psi, x, mu_fun, c1, c2, wt = NULL) {
  m <- .mcglm_compute_m_bin(psi, x, mu_fun, c1, c2)
  if (is.null(wt)) return(colMeans(m))
  colSums(wt * m) / sum(wt)
}

#' Compute I_hat(psi) from pre-built xi_hat
#' @keywords internal
.mcglm_compute_Ihat <- function(psi, xi_hat, mu_dot_fun, wt = NULL) {
  n <- nrow(xi_hat)
  eta_tilde <- as.numeric(xi_hat %*% psi)
  w <- mu_dot_fun(eta_tilde)
  if (is.null(wt)) {
    crossprod(xi_hat * w, xi_hat) / n
  } else {
    crossprod(xi_hat * (wt * w), xi_hat) / sum(wt)
  }
}

#' Analytical Jacobian M_hat for binary misclassification
#' @keywords internal
.mcglm_compute_Mhat_bin <- function(psi, x, mu_fun, mu_dot_fun, c1, c2,
                                    wt = NULL) {
  gamma <- psi[1]
  alpha <- psi[-1]
  r <- length(alpha)
  p <- 1L + r
  n <- nrow(x)
  N <- if (is.null(wt)) n else sum(wt)

  eta0 <- as.numeric(x %*% alpha)
  mu_dot1 <- mu_dot_fun(gamma + eta0)
  mu_dot0 <- mu_dot_fun(eta0)
  d_mu_dot <- mu_dot1 - mu_dot0

  M <- matrix(0, p, p)

  if (is.null(wt)) {
    M[1, 1]    <- -mean(c1 * mu_dot1)
    M[1, 2:p]  <- -colMeans(c1 * d_mu_dot * x)
    M[2:p, 1]    <- -colMeans(c2 * mu_dot1 * x)
    M[2:p, 2:p]  <- -crossprod(c2 * d_mu_dot * x, x) / n
  } else {
    M[1, 1]    <- -sum(wt * c1 * mu_dot1) / N
    M[1, 2:p]  <- -colSums(wt * c1 * d_mu_dot * x) / N
    M[2:p, 1]    <- -colSums(wt * c2 * mu_dot1 * x) / N
    M[2:p, 2:p]  <- -crossprod(x * (wt * c2 * d_mu_dot), x) / N
  }

  M
}

#' Normalize eta argument into a list of three length-q vectors
#' @keywords internal
.mcglm_normalize_eta <- function(eta, q) {
  extract_coef <- function(val, name) {
    if (inherits(val, "glm") || inherits(val, "lm"))
      val <- stats::coef(val)
    if (!is.numeric(val) || any(!is.finite(val)))
      stop("eta$", name, " must be numeric and finite.", call. = FALSE)
    val <- unname(as.numeric(val))
    if (length(val) != q)
      stop("eta$", name, " has length ", length(val), ", but w has ", q,
           " columns.", call. = FALSE)
    val
  }
  if (is.list(eta)) {
    if (!all(c("rho", "p01", "p10") %in% names(eta)))
      stop("When eta is a list, it must contain 'rho', 'p01', and 'p10'.",
           call. = FALSE)
    list(
      rho = extract_coef(eta$rho, "rho"),
      p01 = extract_coef(eta$p01, "p01"),
      p10 = extract_coef(eta$p10, "p10")
    )
  } else if (is.numeric(eta)) {
    if (length(eta) != 3L * q || any(!is.finite(eta)))
      stop("When eta is a numeric vector, it must have length 3 * ncol(w) = ",
           3L * q, ".", call. = FALSE)
    eta <- unname(as.numeric(eta))
    list(
      rho = eta[seq_len(q)],
      p01 = eta[q + seq_len(q)],
      p10 = eta[2L * q + seq_len(q)]
    )
  } else {
    stop("eta must be a named list with 'rho', 'p01', 'p10' or a numeric vector.",
         call. = FALSE)
  }
}

#' Compute c1 and c2 from regressed parameters eta and covariates w
#'
#' Evaluates rho(x; eta_rho) = expit(w' eta_rho),
#' p01(x; eta_01) = expit(w' eta_01), and
#' p10(x; eta_10) = expit(w' eta_10), and computes
#' c1(x) = p01(x) * (1 - rho(x)) and
#' c2(x) = c1(x) - p10(x) * rho(x).
#'
#' @param eta Named list with elements 'rho', 'p01', 'p10' or numeric vector
#'   of length 3 * ncol(w).
#' @param w Covariate matrix (n x q).
#' @return List with components c1, c2, rho, p01, p10, eta.
#' @keywords internal
.mcglm_compute_c_bin <- function(eta, w) {
  w <- as.matrix(w)
  q <- ncol(w)
  eta_norm <- .mcglm_normalize_eta(eta, q)
  rho <- as.numeric(stats::plogis(w %*% eta_norm$rho))
  p01 <- as.numeric(stats::plogis(w %*% eta_norm$p01))
  p10 <- as.numeric(stats::plogis(w %*% eta_norm$p10))
  c1 <- p01 * (1 - rho)
  c2 <- c1 - p10 * rho
  list(c1 = c1, c2 = c2, rho = rho, p01 = p01, p10 = p10, eta = eta_norm)
}

# ---- Multicategory misclassification helpers ----

#' Compute per-observation drift m_i(psi) for multicategory misclassification
#' @keywords internal
.mcglm_compute_m_multi <- function(psi, x, K, mu_fun, Pi, pi_z) {
  n  <- nrow(x)
  s  <- K - 1
  r  <- ncol(x)
  p  <- s + r

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):p]
  eta_base <- as.numeric(x %*% alpha)

  mu_mat <- sapply(seq_len(K), function(ell) {
    mu_fun(eta_base + gamma[ell])
  })

  m <- matrix(0, n, p)

  for (k in seq_len(s)) {
    for (ell in seq_len(K)) {
      coeff <- pi_z[ell] * Pi[k + 1, ell]
      m[, k] <- m[, k] + coeff * (mu_mat[, ell] - mu_mat[, k + 1])
    }
  }

  for (ell in seq_len(K)) {
    for (j in seq_len(K)) {
      coeff <- pi_z[ell] * Pi[j, ell]
      diff_mu <- mu_mat[, ell] - mu_mat[, j]
      m[, (s + 1):p] <- m[, (s + 1):p] + coeff * diff_mu * x
    }
  }

  m
}

#' @keywords internal
.mcglm_compute_mhat_multi <- function(psi, x, K, mu_fun, Pi, pi_z, wt = NULL) {
  m <- .mcglm_compute_m_multi(psi, x, K, mu_fun, Pi, pi_z)
  if (is.null(wt)) return(colMeans(m))
  colSums(wt * m) / sum(wt)
}

#' Compute I_hat(psi) for multicategory case
#' @keywords internal
.mcglm_compute_Ihat_multi <- function(psi, xi_hat, z_hat, K, mu_dot_fun,
                                      wt = NULL) {
  n <- nrow(xi_hat)
  s <- K - 1
  r <- ncol(xi_hat) - s
  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base <- as.numeric(xi_hat[, (s + 1):(s + r), drop = FALSE] %*% alpha)

  eta_tilde <- eta_base + gamma[z_hat + 1]
  w <- mu_dot_fun(eta_tilde)
  if (is.null(wt)) {
    crossprod(xi_hat * w, xi_hat) / n
  } else {
    crossprod(xi_hat * (wt * w), xi_hat) / sum(wt)
  }
}

#' Analytical Jacobian M_hat for multicategory drift
#'
#' Closed-form block expression for the K-class Jacobian, with
#' B = Pi column-scaled by pi_z and w = rowSums(B) = Pr(Z_hat = j).
#' Replaces the forward-difference numerical Jacobian.
#'
#' @keywords internal
.mcglm_compute_Mhat_multi_analytical <- function(psi, x, K, mu_dot_fun, Pi, pi_z,
                                                  wt = NULL) {
  s <- K - 1L
  r <- ncol(x)
  p <- s + r
  n <- nrow(x)
  N <- if (is.null(wt)) n else sum(wt)

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1L):p]
  eta_base <- as.numeric(x %*% alpha)

  dmu_mat <- vapply(seq_len(K),
                    function(ell) mu_dot_fun(eta_base + gamma[ell]),
                    numeric(n))  # n x K

  B <- sweep(Pi, 2, pi_z, "*")     # B[j, l] = Pi[j, l] * pi_z[l]
  w <- rowSums(B)                   # marginal P(Z_hat = j)

  if (is.null(wt)) {
    mean_dmu <- colMeans(dmu_mat)
  } else {
    mean_dmu <- colSums(wt * dmu_mat) / N
  }

  # Block 1: d m_gamma / d gamma  (s x s)
  M_gg <- matrix(0, s, s)
  for (k in seq_len(s)) {
    for (t in seq_len(s)) {
      if (t != k) {
        M_gg[k, t] <- B[k + 1L, t + 1L] * mean_dmu[t + 1L]
      } else {
        M_gg[k, t] <- (B[k + 1L, k + 1L] - w[k + 1L]) * mean_dmu[k + 1L]
      }
    }
  }

  # Block 2: d m_gamma / d alpha  (s x r)
  dmu_dot <- dmu_mat %*% t(B[2:K, , drop = FALSE])             # n x s
  dmu_k   <- dmu_mat[, 2:K, drop = FALSE]                       # n x s
  s_mat   <- dmu_dot - sweep(dmu_k, 2, w[2:K], "*")
  if (is.null(wt)) {
    M_ga <- crossprod(s_mat, x) / N
  } else {
    M_ga <- crossprod(s_mat, wt * x) / N
  }

  # Block 3: d m_alpha / d gamma  (r x s)
  coeff <- pi_z - w
  M_ag <- matrix(0, r, s)
  for (t in seq_len(s)) {
    dmu_t <- dmu_mat[, t + 1L]
    if (is.null(wt)) {
      M_ag[, t] <- coeff[t + 1L] * colMeans(x * dmu_t)
    } else {
      M_ag[, t] <- coeff[t + 1L] * colSums(wt * x * dmu_t) / N
    }
  }

  # Block 4: d m_alpha / d alpha  (r x r)
  t_vec <- drop(dmu_mat %*% coeff)
  if (is.null(wt)) {
    M_aa <- crossprod(x, x * t_vec) / N
  } else {
    M_aa <- crossprod(x, x * (wt * t_vec)) / N
  }

  M <- matrix(0, p, p)
  M[1L:s, 1L:s]               <- M_gg
  M[1L:s, (s + 1L):p]         <- M_ga
  M[(s + 1L):p, 1L:s]         <- M_ag
  M[(s + 1L):p, (s + 1L):p]   <- M_aa
  M
}

#' Numerical Jacobian M_hat for multicategory drift (forward differences)
#' @keywords internal
.mcglm_compute_Mhat_multi_numerical <- function(psi, x, K, mu_fun, Pi, pi_z,
                                                wt = NULL) {
  p <- length(psi)
  M <- matrix(0, p, p)
  h <- 1e-7
  m0 <- .mcglm_compute_mhat_multi(psi, x, K, mu_fun, Pi, pi_z, wt = wt)
  for (j in seq_len(p)) {
    psi_h <- psi
    psi_h[j] <- psi_h[j] + h
    m1 <- .mcglm_compute_mhat_multi(psi_h, x, K, mu_fun, Pi, pi_z, wt = wt)
    M[, j] <- (m1 - m0) / h
  }
  M
}

#' Dispatch wrapper: analytical (default) or numerical Jacobian for K-class drift
#' @keywords internal
.mcglm_compute_Mhat_multi <- function(psi, x, K, mu_fun, Pi, pi_z, wt = NULL,
                                      jacobian = c("analytical", "numerical"),
                                      mu_dot_fun = NULL) {
  jacobian <- match.arg(jacobian)
  if (jacobian == "analytical") {
    if (is.null(mu_dot_fun))
      stop("Internal error: analytical Jacobian requires mu_dot_fun.")
    .mcglm_compute_Mhat_multi_analytical(psi, x, K, mu_dot_fun, Pi, pi_z,
                                          wt = wt)
  } else {
    .mcglm_compute_Mhat_multi_numerical(psi, x, K, mu_fun, Pi, pi_z, wt = wt)
  }
}
