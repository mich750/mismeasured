# ---------------------------------------------------------------------------
# Sandwich variance estimators for mcglm methods
# ---------------------------------------------------------------------------

#' Sandwich variance for the naive estimator
#' @keywords internal
.mcglm_vcov_naive <- function(psi, y, xi_hat, family, wt = NULL) {
  fam    <- .normalize_family(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)
  eta    <- as.numeric(xi_hat %*% psi)
  mu_val <- fam$linkinv(eta)
  w      <- fam$mu.eta(eta)
  eps    <- y - mu_val

  if (is.null(wt)) {
    A <- crossprod(xi_hat * w, xi_hat) / N
    C <- crossprod(xi_hat * eps, xi_hat * eps) / N
  } else {
    A <- crossprod(xi_hat * (wt * w), xi_hat) / N
    C <- crossprod(xi_hat * (wt * eps), xi_hat * eps) / N
  }
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / N
}

#' Sandwich variance for BCA/BCM estimators (binary)
#' @keywords internal
.mcglm_vcov_bc_bin <- function(psi_bc, y, xi_hat, x, family,
                               p01 = NULL, p10 = NULL, pi_z = NULL,
                               c1 = NULL, c2 = NULL,
                               psi_naive = NULL,
                               type = c("bca", "bcm"), corrected = FALSE,
                               wt = NULL) {
  if (!corrected) return(.mcglm_vcov_naive(psi_bc, y, xi_hat, family, wt = wt))

  type   <- match.arg(type)
  fam    <- .normalize_family(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)
  p      <- length(psi_bc)

  if (is.null(psi_naive)) psi_naive <- psi_bc
  eta    <- as.numeric(xi_hat %*% psi_naive)
  w      <- fam$mu.eta(eta)
  resid  <- y - fam$linkinv(eta)

  if (is.null(wt)) {
    A_hat <- crossprod(xi_hat * w, xi_hat) / N
  } else {
    A_hat <- crossprod(xi_hat * (wt * w), xi_hat) / N
  }
  A_inv <- solve(A_hat)
  c1_loc <- if (!is.null(c1)) c1 else p01 * (1 - pi_z)
  c2_loc <- if (!is.null(c2)) c2 else p01 * (1 - pi_z) - p10 * pi_z
  if (is.null(c1_loc) || is.null(c2_loc) ||
      length(c1_loc) == 0L || length(c2_loc) == 0L)
    stop("Corrected variance for 'bca'/'bcm' requires either ",
         "(p01, p10, pi_z) or (c1, c2).")
  M_hat <- .mcglm_compute_Mhat_bin(psi_naive, x, fam$linkinv, fam$mu.eta, c1_loc,
                                    c2_loc, wt = wt)
  m_mat <- .mcglm_compute_m_bin(psi_naive, x, fam$linkinv, c1_loc, c2_loc)
  if (is.null(wt)) {
    m_bar <- colMeans(m_mat)
  } else {
    m_bar <- colSums(wt * m_mat) / N
  }

  if (type == "bca") {
    G     <- diag(p) - A_inv %*% M_hat
    H_inv <- A_inv
  } else {
    H     <- A_hat + M_hat
    H_inv <- solve(H)
    G     <- diag(p) - H_inv %*% M_hat
  }

  score_mat  <- xi_hat * resid
  centered_m <- sweep(m_mat, 2, m_bar)
  L1 <- A_inv %*% t(G)
  L2 <- t(H_inv)

  if (is.null(wt)) {
    S_ss <- crossprod(score_mat) / N^2
    S_mm <- crossprod(centered_m) / N^2
    S_sm <- crossprod(score_mat, centered_m) / N^2
  } else {
    S_ss <- crossprod(score_mat * wt, score_mat) / N^2
    S_mm <- crossprod(centered_m * wt, centered_m) / N^2
    S_sm <- crossprod(score_mat * wt, centered_m) / N^2
  }

  t(L1) %*% S_ss %*% L1 + t(L2) %*% S_mm %*% L2 -
    t(L1) %*% S_sm %*% L2 - t(L2) %*% t(S_sm) %*% L1
}

#' Sandwich variance for corrected-score estimator (binary)
#' @keywords internal
.mcglm_vcov_cs_bin <- function(psi, y, xi_hat, x, family, p01, p10, pi_z,
                               c1 = NULL, c2 = NULL,
                               wt = NULL, validation = NULL) {
  fam    <- .normalize_family(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)

  eta_tilde <- as.numeric(xi_hat %*% psi)
  resid     <- y - fam$linkinv(eta_tilde)
  c1_loc <- if (!is.null(c1)) c1 else p01 * (1 - pi_z)
  c2_loc <- if (!is.null(c2)) c2 else p01 * (1 - pi_z) - p10 * pi_z
  if (is.null(c1_loc) || is.null(c2_loc) ||
      length(c1_loc) == 0L || length(c2_loc) == 0L)
    stop("Variance for 'cs' requires either (p01, p10, pi_z) or (c1, c2).")
  m_mat     <- .mcglm_compute_m_bin(psi, x, fam$linkinv, c1_loc, c2_loc)

  phi_mat <- xi_hat * resid - m_mat

  if (is.null(wt)) {
    S <- crossprod(phi_mat) / N
  } else {
    S <- crossprod(phi_mat * wt, phi_mat) / N
  }

  if (!is.null(validation))
    S <- .mcglm_cs_validation_meat(S, phi_mat, psi, x, 2L,
                                    fam$linkinv, validation)

  I_hat <- .mcglm_compute_Ihat(psi, xi_hat, fam$mu.eta, wt = wt)
  M_hat <- .mcglm_compute_Mhat_bin(psi, x, fam$linkinv, fam$mu.eta, c1_loc, c2_loc,
                                    wt = wt)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}


# Same fall-through for the BCA/BCM corrected sandwich.
# (The earlier definition above already exists; here we extend the binary
#  helper signature to accept c1/c2 the same way.)


# ---- Multicategory variance estimators ----

#' Sandwich variance for naive estimator (multicategory)
#' @keywords internal
.mcglm_vcov_naive_multi <- function(psi, y, xi_hat, z_hat, x, K, family,
                                    wt = NULL) {
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  w   <- fam$mu.eta(eta_tilde)
  eps <- y - fam$linkinv(eta_tilde)

  if (is.null(wt)) {
    A <- crossprod(xi_hat * w, xi_hat) / N
    C <- crossprod(xi_hat * eps, xi_hat * eps) / N
  } else {
    A <- crossprod(xi_hat * (wt * w), xi_hat) / N
    C <- crossprod(xi_hat * (wt * eps), xi_hat * eps) / N
  }
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / N
}

#' Sandwich variance for BCA/BCM estimators (multicategory)
#' @keywords internal
.mcglm_vcov_bc_multi <- function(psi_bc, y, xi_hat, z_hat, x, K, family,
                                 Pi = NULL, pi_z = NULL,
                                 psi_naive = NULL,
                                 type = c("bca", "bcm"), corrected = FALSE,
                                 wt = NULL,
                                 jacobian = c("analytical", "numerical")) {
  jacobian <- match.arg(jacobian)
  if (!corrected)
    return(.mcglm_vcov_naive_multi(psi_bc, y, xi_hat, z_hat, x, K, family,
                                    wt = wt))

  type <- match.arg(type)
  fam  <- .normalize_family(family)
  n    <- length(y)
  N    <- if (is.null(wt)) n else sum(wt)
  s    <- K - 1
  r    <- ncol(x)
  p    <- s + r

  if (is.null(psi_naive)) psi_naive <- psi_bc

  gamma <- c(0, psi_naive[seq_len(s)])
  alpha <- psi_naive[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  w     <- fam$mu.eta(eta_tilde)
  resid <- y - fam$linkinv(eta_tilde)

  if (is.null(wt)) {
    A_hat <- crossprod(xi_hat * w, xi_hat) / N
  } else {
    A_hat <- crossprod(xi_hat * (wt * w), xi_hat) / N
  }
  A_inv <- solve(A_hat)
  M_hat <- .mcglm_compute_Mhat_multi(psi_naive, x, K, fam$linkinv, Pi, pi_z,
                                      wt = wt, jacobian = jacobian,
                                      mu_dot_fun = fam$mu.eta)
  m_mat <- .mcglm_compute_m_multi(psi_naive, x, K, fam$linkinv, Pi, pi_z)
  if (is.null(wt)) {
    m_bar <- colMeans(m_mat)
  } else {
    m_bar <- colSums(wt * m_mat) / N
  }

  if (type == "bca") {
    G     <- diag(p) - A_inv %*% M_hat
    H_inv <- A_inv
  } else {
    H     <- A_hat + M_hat
    H_inv <- solve(H)
    G     <- diag(p) - H_inv %*% M_hat
  }

  score_mat  <- xi_hat * resid
  centered_m <- sweep(m_mat, 2, m_bar)
  L1 <- A_inv %*% t(G)
  L2 <- t(H_inv)

  if (is.null(wt)) {
    S_ss <- crossprod(score_mat) / N^2
    S_mm <- crossprod(centered_m) / N^2
    S_sm <- crossprod(score_mat, centered_m) / N^2
  } else {
    S_ss <- crossprod(score_mat * wt, score_mat) / N^2
    S_mm <- crossprod(centered_m * wt, centered_m) / N^2
    S_sm <- crossprod(score_mat * wt, centered_m) / N^2
  }

  t(L1) %*% S_ss %*% L1 + t(L2) %*% S_mm %*% L2 -
    t(L1) %*% S_sm %*% L2 - t(L2) %*% t(S_sm) %*% L1
}

#' Sandwich variance for corrected-score estimator (multicategory)
#' @keywords internal
.mcglm_vcov_cs_multi <- function(psi, y, xi_hat, z_hat, x, K, family, Pi, pi_z,
                                 wt = NULL,
                                 jacobian = c("analytical", "numerical"),
                                 validation = NULL) {
  jacobian <- match.arg(jacobian)
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  resid <- y - fam$linkinv(eta_tilde)
  m_mat <- .mcglm_compute_m_multi(psi, x, K, fam$linkinv, Pi, pi_z)

  phi_mat <- xi_hat * resid - m_mat
  if (is.null(wt)) {
    S <- crossprod(phi_mat) / N
  } else {
    S <- crossprod(phi_mat * wt, phi_mat) / N
  }

  if (!is.null(validation))
    S <- .mcglm_cs_validation_meat(S, phi_mat, psi, x, K,
                                    fam$linkinv, validation)

  I_hat <- .mcglm_compute_Ihat_multi(psi, xi_hat, z_hat, K, fam$mu.eta,
                                      wt = wt)
  M_hat <- .mcglm_compute_Mhat_multi(psi, x, K, fam$linkinv, Pi, pi_z, wt = wt,
                                      jacobian = jacobian,
                                      mu_dot_fun = fam$mu.eta)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}

# Validate the sampling information before the variance tryCatch in mcglm.
# Joint cells B[j,l] = P(proxy=j, true=l) are sufficient for both K=2 and
# K>2. For K=2, B[2,1]=c1 and B[1,2]=c1-c2; this is exactly the binary
# (pi,p01,p10) influence-function formula by the chain rule.
.mcglm_prepare_cs_validation <- function(validation, z_hat, K, wt,
                                          c1, c2, Pi, pi_z,
                                          x = NULL, w = NULL, eta = NULL) {
  if (!is.null(wt) && any(wt != 1))
    stop("validation currently requires unweighted observations or unit weights.",
         call. = FALSE)
  if (!is.list(validation) || is.null(validation$z) ||
      any(!names(validation) %in% c("z", "z_hat", "index", "w", "x")))
    stop("validation must be a list with z and z_hat (external) or index (internal).",
         call. = FALSE)
  z <- validation$z
  nv <- length(z)
  check_codes <- function(v) {
    is.numeric(v) && is.null(dim(v)) && length(v) == nv &&
      all(is.finite(v)) && all(v == floor(v)) && all(v >= 0 & v < K)
  }
  if (nv < 2L || !check_codes(z))
    stop("validation$z must contain at least two codes in 0, ..., K-1.",
         call. = FALSE)
  index <- validation$index
  proxy <- validation$z_hat
  if (!is.null(index)) {
    if (!is.numeric(index) || !is.null(dim(index)) || length(index) != nv ||
        any(!is.finite(index)) || any(index != floor(index)) ||
        any(index < 1 | index > length(z_hat)) || anyDuplicated(index))
      stop("validation$index must contain distinct regression row numbers, one per z.",
           call. = FALSE)
    if (!is.null(proxy) &&
        (!check_codes(proxy) || any(proxy != z_hat[index])))
      stop("validation$z_hat does not match z_hat at validation$index.",
           call. = FALSE)
    proxy <- z_hat[index]
  }
  if (!check_codes(proxy))
    stop("validation$z_hat must contain one code in 0, ..., K-1 per z.",
         call. = FALSE)

  is_regressed <- (K == 2L) && (!is.null(eta) || !is.null(w) ||
    !is.null(validation$w) || !is.null(validation$x) ||
    (is.null(c1) && is.null(Pi) && is.null(pi_z)) ||
    (length(c1) > 1L || length(c2) > 1L))

  if (!is_regressed || K > 2L) {
    cells <- proxy + 1L + K * z  # column-major vec(B), as in the paper
    b <- tabulate(cells, nbins = K * K) / nv
    B <- matrix(b, K, K)
    if (any(colSums(B) == 0))
      stop("Every true category must occur in validation.", call. = FALSE)
    expected <- if (K == 2L) c(B[2, 1], B[2, 1] - B[1, 2]) else b
    supplied <- if (K == 2L) c(c1, c2) else as.numeric(sweep(Pi, 2, pi_z, "*"))
    if (length(supplied) != length(expected) || any(!is.finite(supplied)) ||
        max(abs(supplied - expected)) > 1e-8)
      stop("CS probabilities must be empirical estimates from the supplied validation sample; ",
           "supply pi_z explicitly (or c1/c2 for K=2).", call. = FALSE)
    return(list(b = b, cells = cells, index = index, n = nv, regressed = FALSE,
                type = if (is.null(index)) "external" else "internal"))
  }

  # Regressed parameters eta (K = 2)
  idx0 <- which(z == 0L)
  idx1 <- which(z == 1L)
  if (length(idx0) == 0L || length(idx1) == 0L)
    stop("Every true category must occur in validation.", call. = FALSE)

  if (!is.null(validation$w)) {
    w_val <- as.matrix(validation$w)
  } else if (!is.null(validation$x)) {
    w_val <- as.matrix(validation$x)
  } else if (!is.null(index)) {
    if (!is.null(w)) {
      w_val <- as.matrix(w)[index, , drop = FALSE]
    } else {
      w_val <- as.matrix(x)[index, , drop = FALSE]
    }
  } else {
    if (!is.null(w) || (!is.null(x) && ncol(as.matrix(x)) > 1L)) {
      stop("For external validation with regressed eta, supply w in validation list.",
           call. = FALSE)
    } else {
      w_val <- matrix(1, nv, 1)
    }
  }

  q <- ncol(w_val)
  n <- length(z_hat)
  if (!is.null(w)) {
    w_main <- as.matrix(w)
  } else if (!is.null(x)) {
    w_main <- as.matrix(x)
  } else {
    w_main <- matrix(1, n, 1)
  }

  if (ncol(w_main) != q)
    stop("w in validation (", q, " cols) does not match w in regression (",
         ncol(w_main), " cols).", call. = FALSE)

  if (is.null(eta)) {
    fit_rho <- stats::glm.fit(w_val, z, family = stats::binomial())
    if (!fit_rho$converged) warning("prevalence model did not converge", call. = FALSE)

    fit_01 <- stats::glm.fit(w_val[idx0, , drop = FALSE], proxy[idx0], family = stats::binomial())
    if (!fit_01$converged) warning("misclassification model p01 did not converge", call. = FALSE)

    fit_10 <- stats::glm.fit(w_val[idx1, , drop = FALSE], 1L - proxy[idx1], family = stats::binomial())
    if (!fit_10$converged) warning("misclassification model p10 did not converge", call. = FALSE)

    eta_rho <- unname(stats::coef(fit_rho))
    eta_01  <- unname(stats::coef(fit_01))
    eta_10  <- unname(stats::coef(fit_10))
    eta <- list(rho = eta_rho, p01 = eta_01, p10 = eta_10)
  } else {
    eta <- .mcglm_normalize_eta(eta, q)
    eta_rho <- eta$rho
    eta_01  <- eta$p01
    eta_10  <- eta$p10
  }

  rho_val <- as.numeric(stats::plogis(w_val %*% eta_rho))
  p01_val <- as.numeric(stats::plogis(w_val %*% eta_01))
  p10_val <- as.numeric(stats::plogis(w_val %*% eta_10))

  v_rho_val <- rho_val * (1 - rho_val)
  v_01_val  <- p01_val * (1 - p01_val)
  v_10_val  <- p10_val * (1 - p10_val)

  H_rho <- crossprod(w_val * v_rho_val, w_val) / nv
  H_01  <- crossprod(w_val * ((1 - z) * v_01_val), w_val) / nv
  H_10  <- crossprod(w_val * (z * v_10_val), w_val) / nv

  s_rho <- (w_val * (z - rho_val)) %*% solve(H_rho)
  s_01  <- (w_val * ((1 - z) * (proxy - p01_val))) %*% solve(H_01)
  s_10  <- (w_val * (z * (1 - proxy - p10_val))) %*% solve(H_10)

  s_mat <- cbind(s_rho, s_01, s_10)

  rho_main <- as.numeric(stats::plogis(w_main %*% eta_rho))
  p01_main <- as.numeric(stats::plogis(w_main %*% eta_01))
  p10_main <- as.numeric(stats::plogis(w_main %*% eta_10))

  c1_main <- p01_main * (1 - rho_main)
  c2_main <- c1_main - p10_main * rho_main

  list(
    type = if (is.null(index)) "external" else "internal",
    n = nv,
    index = index,
    regressed = TRUE,
    eta = eta,
    w = w_main,
    w_val = w_val,
    s = s_mat,
    c1 = c1_main,
    c2 = c2_main,
    rho = rho_main,
    p01 = p01_main,
    p10 = p10_main
  )
}

# Equations Sigmahat-cs-external/internal and D-multicategory-explicit in
# GLM bias correction.tex and Estimated probability corrected score.txt.
# All quantities here are on the sqrt(n) scale;
# the existing sandwich callers divide the resulting covariance by n.
.mcglm_cs_validation_meat <- function(S, phi_mat, psi, x, K, mu_fun,
                                       validation) {
  n <- nrow(x)
  nv <- validation$n

  if (isTRUE(validation$regressed)) {
    w <- validation$w
    q <- ncol(w)
    p <- length(psi)

    gamma <- psi[1L]
    alpha <- psi[-1L]
    eta0 <- as.numeric(x %*% alpha)
    delta <- mu_fun(gamma + eta0) - mu_fun(eta0)

    rho <- validation$rho
    p01 <- validation$p01
    p10 <- validation$p10

    v_rho <- rho * (1 - rho)
    v_01  <- p01 * (1 - p01)
    v_10  <- p10 * (1 - p10)

    D_rho <- matrix(0, p, q)
    D_rho[1L, ] <- -colMeans(delta * p01 * v_rho * w)
    D_rho[-1L, ] <- -crossprod(x * (delta * (p01 + p10) * v_rho), w) / n

    D_01 <- matrix(0, p, q)
    D_01[1L, ] <- colMeans(delta * (1 - rho) * v_01 * w)
    D_01[-1L, ] <- crossprod(x * (delta * (1 - rho) * v_01), w) / n

    D_10 <- matrix(0, p, q)
    D_10[1L, ] <- 0
    D_10[-1L, ] <- -crossprod(x * (delta * rho * v_10), w) / n

    D_mat <- cbind(D_rho, D_01, D_10)

    U <- validation$s %*% t(D_mat)

    S <- S + (n / nv) * crossprod(U) / nv
    if (!is.null(validation$index)) {
      phi_v <- phi_mat[validation$index, , drop = FALSE]
      phi_v <- sweep(phi_v, 2, colMeans(phi_v))
      cross <- crossprod(phi_v, U) / nv
      S <- S + cross + t(cross)
    }
    return(S)
  }

  s <- K - 1L
  gamma <- c(0, psi[seq_len(s)])
  eta <- as.numeric(x %*% psi[-seq_len(s)])
  mu <- vapply(gamma, function(g) mu_fun(eta + g), numeric(n))
  D <- matrix(0, length(psi), K * K)
  for (ell in seq_len(K)) {
    for (j in seq_len(K)) {
      column <- j + K * (ell - 1L)
      delta <- mu[, j] - mu[, ell]
      if (j > 1L) D[j - 1L, column] <- mean(delta)
      D[-seq_len(s), column] <- colMeans(x * delta)
    }
  }

  # Row i equals D s_i. Computing in coefficient space avoids allocating
  # an n_V by K^2 matrix and never inverts the singular cell covariance.
  influence <- sweep(t(D[, validation$cells, drop = FALSE]), 2,
                     as.numeric(D %*% validation$b))
  S <- S + (n / nv) * crossprod(influence) / nv
  if (!is.null(validation$index)) {
    phi_v <- phi_mat[validation$index, , drop = FALSE]
    phi_v <- sweep(phi_v, 2, colMeans(phi_v))
    cross <- crossprod(phi_v, influence) / nv
    S <- S + cross + t(cross)
  }
  S
}
