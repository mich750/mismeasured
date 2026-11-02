# Akazawa-Kinukawa-Nakamura (1998) corrected-score construction for
# misclassified categorical covariates.
#
# Reference: Akazawa, K., Kinukawa, N. and Nakamura, T. (1998).
# A note on the corrected score function adjusting for misclassification.
# J. Japan Statist. Soc., 28(1), 115-123.
#
# Setup matches the package's mcglm() conventions:
#   - Latent Z in {0, ..., K-1}, surrogate Z_hat in {0, ..., K-1}.
#   - Pi: K x K column-stochastic, Pi[j+1, l+1] = P(Z_hat = j | Z = l).
#   - psi = (gamma_1, ..., gamma_{K-1}, alpha_1, ..., alpha_r), gamma_0 := 0.
#   - x: n x r precisely-observed covariates (typically including intercept).
#
# AKN98 builds a (K-1)-vector unbiased surrogate x_akn such that
# E[x_akn | Z = l] = e_l for l = 1, ..., K-1 and = 0 for l = 0. It then
# replaces the dummy-encoded latent in the GLM log-likelihood by x_akn,
# using the identity f*(x_akn) = sum_{l=0}^{K-1} f(e_l) * x_akn_l with
# x_akn_0 := 1 - sum(x_akn). Differentiating gives an unbiased score.
#
# Unlike the package's drift-corrected "cs", the AKN construction does
# NOT require the latent prevalence vector pi_z -- only Pi is needed.

#' Build inverse Q matrix and baseline probability vector
#'
#' Q = [p_1 - p_0, ..., p_(K-1) - p_0] is (K-1) x (K-1) where
#' p_l = (Pi[2, l+1], ..., Pi[K, l+1])^T = E[u | Z = l] and u is the
#' (K-1)-dim dummy encoding of Z_hat with baseline 0.
#'
#' @keywords internal
.akn_build_Q <- function(Pi, K) {
  Pi <- as.matrix(Pi)
  if (nrow(Pi) != K || ncol(Pi) != K)
    stop(".akn_build_Q: Pi must be K x K (got ",
         nrow(Pi), " x ", ncol(Pi), ").")
  P_mat <- Pi[-1, , drop = FALSE]               # (K-1) x K
  p0    <- P_mat[, 1]                            # length K-1
  Q     <- P_mat[, -1, drop = FALSE] - p0        # (K-1) x (K-1)
  cond_Q <- tryCatch(kappa(Q, exact = FALSE),
                     error = function(e) Inf)
  if (!is.finite(cond_Q) || cond_Q > 1e8)
    warning("cs_akn: Q is near-singular (condition number ",
            signif(cond_Q, 3), "); estimates may be unstable.")
  Q_inv <- tryCatch(solve(Q),
                    error = function(e)
                      stop("cs_akn: Q matrix [p_1 - p_0, ..., p_{K-1} - p_0] ",
                           "is singular; AKN98 surrogate is not identified."))
  list(Q_inv = Q_inv, p0 = p0)
}

#' Build AKN98 unbiased surrogate
#'
#' @return n x (K-1) matrix x_akn with E[x_akn[i, ] | Z_i = l] = e_l for
#'   l = 1, ..., K-1 and = 0 for l = 0.
#' @keywords internal
.akn_unbiased_surrogate <- function(z_hat, Q_inv, p0, K) {
  n <- length(z_hat)
  s <- K - 1L
  u <- matrix(0, n, s)
  for (j in seq_len(s)) u[, j] <- as.numeric(z_hat == j)
  centered <- sweep(u, 2, p0, "-")
  centered %*% t(Q_inv)
}

#' Per-observation AKN98 corrected score
#'
#' @return n x p matrix where p = (K-1) + r.
#' @keywords internal
.akn_phi <- function(psi, y, x_akn, x, K, fam) {
  s <- K - 1L
  r <- ncol(x)
  n <- length(y)
  gamma <- psi[seq_len(s)]
  alpha <- psi[s + seq_len(r)]
  alpha_lin <- as.numeric(x %*% alpha)

  mu_mat <- matrix(0, n, s)
  for (i in seq_len(s)) mu_mat[, i] <- fam$linkinv(gamma[i] + alpha_lin)
  mu_0   <- fam$linkinv(alpha_lin)
  x_akn0 <- 1 - rowSums(x_akn)

  resid_gamma <- (y - mu_mat) * x_akn                     # n x (K-1)

  weighted_mu <- mu_0 * x_akn0 + rowSums(mu_mat * x_akn)  # length n
  resid_alpha <- y - weighted_mu
  alpha_score <- resid_alpha * x                          # n x r

  cbind(resid_gamma, alpha_score)
}

#' AKN98 corrected information
#'
#' Returned matrix is the *mean* per observation; the Jacobian of the
#' mean score colMeans(phi) wrt psi is -.akn_information() (negative
#' because the score is residual-based).
#'
#' @keywords internal
.akn_information <- function(psi, y, x_akn, x, K, fam, wt = NULL) {
  s <- K - 1L
  r <- ncol(x)
  n <- length(y)
  N <- if (is.null(wt)) n else sum(wt)
  gamma <- psi[seq_len(s)]
  alpha <- psi[s + seq_len(r)]
  alpha_lin <- as.numeric(x %*% alpha)

  mud_mat <- matrix(0, n, s)
  for (i in seq_len(s)) mud_mat[, i] <- fam$mu.eta(gamma[i] + alpha_lin)
  mud_0   <- fam$mu.eta(alpha_lin)
  x_akn0  <- 1 - rowSums(x_akn)

  p <- s + r
  I <- matrix(0, p, p)

  if (is.null(wt)) {
    diag_gg <- colSums(mud_mat * x_akn) / N
  } else {
    diag_gg <- colSums(mud_mat * x_akn * wt) / N
  }
  for (i in seq_len(s)) I[i, i] <- diag_gg[i]

  for (i in seq_len(s)) {
    if (is.null(wt)) {
      block <- crossprod(mud_mat[, i] * x_akn[, i], x) / N
    } else {
      block <- crossprod((mud_mat[, i] * x_akn[, i]) * wt, x) / N
    }
    I[i, s + seq_len(r)] <- as.numeric(block)
    I[s + seq_len(r), i] <- as.numeric(block)
  }

  W <- mud_0 * x_akn0 + rowSums(mud_mat * x_akn)
  if (is.null(wt)) {
    I[s + seq_len(r), s + seq_len(r)] <- crossprod(x, x * W) / N
  } else {
    I[s + seq_len(r), s + seq_len(r)] <- crossprod(x, x * W * wt) / N
  }

  I
}

#' Closed-form Gaussian AKN98 estimator
#' @keywords internal
.akn_fit_gaussian_closed <- function(y, x_akn, x, K, wt = NULL) {
  s <- K - 1L
  r <- ncol(x)
  if (is.null(wt)) {
    XaX <- crossprod(x, x_akn)
    XX  <- crossprod(x)
    Xay <- crossprod(x_akn, y)
    Xy  <- crossprod(x, y)
    diag_block <- colSums(x_akn)
  } else {
    XaX <- crossprod(x * wt, x_akn)
    XX  <- crossprod(x * wt, x)
    Xay <- crossprod(x_akn * wt, y)
    Xy  <- crossprod(x * wt, y)
    diag_block <- colSums(x_akn * wt)
  }

  A <- matrix(0, s + r, s + r)
  diag(A)[seq_len(s)] <- diag_block
  A[seq_len(s), s + seq_len(r)] <- t(XaX)
  A[s + seq_len(r), seq_len(s)] <- XaX
  A[s + seq_len(r), s + seq_len(r)] <- XX

  as.numeric(solve(A, c(Xay, Xy)))
}

#' Fit AKN98 corrected-score estimator
#'
#' Works for any K >= 2 and for Gaussian, Poisson, Binomial families.
#' Multinomial is unsupported (AKN98 does not cover multinomial responses).
#'
#' @keywords internal
.mcglm_fit_cs_akn <- function(psi_init, y, x_akn, x, K, family, wt = NULL) {
  fam <- .normalize_family(family)
  fam_str <- fam$family

  if (identical(fam_str, "gaussian"))
    return(.akn_fit_gaussian_closed(y, x_akn, x, K, wt = wt))

  score_mean <- function(psi) {
    Phi <- .akn_phi(psi, y, x_akn, x, K, fam)
    if (is.null(wt)) colMeans(Phi) else colSums(Phi * wt) / sum(wt)
  }
  score_jac <- function(psi) {
    -.akn_information(psi, y, x_akn, x, K, fam, wt = wt)
  }
  sol <- nleqslv::nleqslv(psi_init, score_mean, jac = score_jac,
                          control = list(maxit = 500, ftol = 1e-12))
  if (sol$termcd > 2) {
    warning("cs_akn solver did not converge (termcd = ", sol$termcd, ")", call. = FALSE)
    return(rep(NA_real_, length(psi_init)))
  }
  sol$x
}

#' Sandwich variance for cs_akn
#' @keywords internal
.mcglm_vcov_cs_akn <- function(psi, y, x_akn, x, K, family, wt = NULL,
                               validation = NULL, w = NULL,
                               p01 = NULL, p10 = NULL, z_hat = NULL) {
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)

  Phi <- .akn_phi(psi, y, x_akn, x, K, fam)
  if (is.null(wt)) {
    S <- crossprod(Phi) / N
  } else {
    S <- crossprod(Phi * wt, Phi) / N
  }

  if (!is.null(validation)) {
    nv <- validation$n
    if (isTRUE(validation$regressed) && K == 2L) {
      p <- length(psi)
      w_use <- if (!is.null(validation$w)) as.matrix(validation$w) else if (!is.null(w)) as.matrix(w) else as.matrix(x)
      q <- ncol(w_use)

      gamma <- psi[1L]
      alpha <- psi[-1L]
      alpha_lin <- as.numeric(x %*% alpha)
      mu1 <- fam$linkinv(gamma + alpha_lin)
      mu0 <- fam$linkinv(alpha_lin)
      delta <- mu1 - mu0

      p01_v <- if (!is.null(validation$p01)) validation$p01 else p01
      p10_v <- if (!is.null(validation$p10)) validation$p10 else p10
      Delta <- 1 - p01_v - p10_v
      Delta <- ifelse(abs(Delta) < 0.05, ifelse(Delta >= 0, 0.05, -0.05), Delta)
      v01 <- p01_v * (1 - p01_v)
      v10 <- p10_v * (1 - p10_v)

      D_rho <- matrix(0, p, q)
      D_01  <- matrix(0, p, q)
      weight_01 <- (1 - p10_v - z_hat) / (Delta^2) * v01
      if (is.null(wt)) {
        D_01[1L, ]  <- -colMeans((y - mu1) * weight_01 * w_use)
        D_01[-1L, ] <- crossprod(x * (delta * weight_01), w_use) / n
      } else {
        D_01[1L, ]  <- -colSums(wt * (y - mu1) * weight_01 * w_use) / N
        D_01[-1L, ] <- crossprod(x * (wt * delta * weight_01), w_use) / N
      }

      D_10  <- matrix(0, p, q)
      weight_10 <- (z_hat - p01_v) / (Delta^2) * v10
      if (is.null(wt)) {
        D_10[1L, ]  <- colMeans((y - mu1) * weight_10 * w_use)
        D_10[-1L, ] <- -crossprod(x * (delta * weight_10), w_use) / n
      } else {
        D_10[1L, ]  <- colSums(wt * (y - mu1) * weight_10 * w_use) / N
        D_10[-1L, ] <- -crossprod(x * (wt * delta * weight_10), w_use) / N
      }

      D_mat <- cbind(D_rho, D_01, D_10)
      U <- validation$s %*% t(D_mat)

      S <- S + (n / nv) * crossprod(U) / nv
      if (!is.null(validation$index)) {
        phi_v <- Phi[validation$index, , drop = FALSE]
        phi_v <- sweep(phi_v, 2, colMeans(phi_v))
        cross <- crossprod(phi_v, U) / nv
        S <- S + cross + t(cross)
      }
    }
  }

  I_akn <- .akn_information(psi, y, x_akn, x, K, fam, wt = wt)
  J     <- -I_akn
  J_inv <- tryCatch(solve(J), error = function(e) MASS::ginv(J))
  J_inv %*% S %*% t(J_inv) / N
}
