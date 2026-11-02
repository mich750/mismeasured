#' Bias-Corrected GLM Estimation with Misclassified Covariates
#'
#' Fits a GLM in which one covariate is observed through a misclassified
#' proxy and applies the bias-correction estimators of Battaglia,
#' Christensen, Hansen and Sacher (2025). Five estimators are available:
#' \emph{naive} (uncorrected proxy-score), \emph{BCA} (additive bias
#' correction), \emph{BCM} (multiplicative bias correction), \emph{CS}
#' (corrected score), and \emph{one-step} (joint mixture likelihood via
#' \pkg{RTMB}). Both binary and multicategory latent regressors are
#' supported.
#'
#' @section Model and notation:
#' Let \eqn{Y_i} denote a scalar response, \eqn{x_i \in \mathbb{R}^r} a
#' vector of cleanly observed covariates (typically including an
#' intercept), and \eqn{Z_i \in \{0, 1, \dots, K-1\}} a latent categorical
#' regressor that is observed only through a proxy \eqn{\hat Z_i}.
#' For the binary case (\eqn{K = 2}) we let
#' \eqn{\xi_i = (Z_i, x_i^\top)^\top} and
#' \eqn{\hat\xi_i = (\hat Z_i, x_i^\top)^\top}; for \eqn{K > 2},
#' \eqn{Z_i} is dummy-encoded with baseline category 0 and
#' \eqn{\xi_i = (d_i^\top, x_i^\top)^\top}.
#' The parameter is \eqn{\psi = (\gamma^\top, \alpha^\top)^\top}, where
#' \eqn{\gamma \in \mathbb{R}^{K-1}} are the coefficients on the
#' (dummy-coded) latent regressor and \eqn{\alpha \in \mathbb{R}^r}
#' the coefficients on \eqn{x_i}. The conditional mean is
#' \eqn{E[Y_i | Z_i, x_i] = \mu(\psi^\top \xi_i)} for
#' \eqn{\mu = g^{-1}} the inverse canonical link.
#'
#' \emph{Misclassification mechanism.}
#' For \eqn{K = 2}, the mechanism is summarized by the false-positive
#' rate \eqn{p_{01} = \Pr(\hat Z = 1 \mid Z = 0)}, the false-negative
#' rate \eqn{p_{10} = \Pr(\hat Z = 0 \mid Z = 1)} and the prevalence
#' \eqn{\pi = \Pr(Z = 1)}. For \eqn{K > 2}, by the \eqn{K \times K}
#' matrix \eqn{\Pi_{j\ell} = \Pr(\hat Z = j - 1 \mid Z = \ell - 1)}
#' (column-stochastic) and the prevalence vector
#' \eqn{(\pi_0, \dots, \pi_{K-1})}. Misclassification is assumed
#' nondifferential: \eqn{\hat Z \perp\!\!\!\perp (Y, x) \mid Z}.
#'
#' \emph{Estimators.}
#' Let \eqn{\hat\psi} solve the proxy score
#' \eqn{\hat U_n(\psi) = n^{-1} \sum_i \hat\xi_i \{Y_i - \mu(\psi^\top \hat\xi_i)\} = 0}
#' (this is the naive GLM estimator). Define the per-observation drift
#' \eqn{m_i(\psi) = (-c_1 \delta_i(\psi),\, -c_2 \delta_i(\psi)\, x_i)^\top}
#' with \eqn{\delta_i(\psi) = \mu(\gamma + \alpha^\top x_i) - \mu(\alpha^\top x_i)},
#' \eqn{c_1 = p_{01}(1 - \pi)} and
#' \eqn{c_2 = p_{01}(1 - \pi) - p_{10} \pi}, and let
#' \eqn{\hat I(\psi) = n^{-1} \sum_i \dot\mu(\psi^\top \hat\xi_i) \hat\xi_i \hat\xi_i^\top}
#' and \eqn{\hat M(\psi) = \partial \hat m / \partial \psi^\top}. Then:
#'
#' \describe{
#'   \item{BCA:}{\eqn{\hat\psi_{\mathrm{bca}} = \hat\psi - \hat I(\hat\psi)^{-1} \hat m(\hat\psi)}}
#'   \item{BCM:}{\eqn{\hat\psi_{\mathrm{bcm}} = \hat\psi - (\hat I(\hat\psi) + \hat M(\hat\psi))^{-1} \hat m(\hat\psi)}}
#'   \item{CS:}{solves \eqn{n^{-1} \sum_i \phi_i(\psi) = 0} with
#'     \eqn{\phi_i(\psi) = \hat\xi_i \{Y_i - \mu(\psi^\top \hat\xi_i)\} - m_i(\psi)}}
#'   \item{one-step:}{maximizes the integrated mixture likelihood
#'     \eqn{\prod_i \sum_\ell \pi_\ell\, \Pi_{\hat z_i+1, \ell+1}\, f(Y_i \mid Z = \ell, \psi)}
#'     by automatic differentiation (\pkg{RTMB}).}
#' }
#'
#' Multicategory analogues replace \eqn{(p_{01}, p_{10}, \pi)} with
#' \eqn{(\Pi, \pi_z)} throughout; see Appendix B of Battaglia et al.
#' (2025) for explicit formulas.
#'
#' \emph{Inference.}
#' All five estimators report an asymptotic variance through
#' \code{\link{vcov.mcglm}}: \eqn{A^{-1} C A^{-1}} for naive / BCA / BCM
#' (Theorems on the two-step expansion and on bias-corrected estimators
#' under drifting misclassification), \eqn{J^{-1} S J^{-\top}} for CS
#' (Z-estimator sandwich, with \eqn{J = -(\hat I + \hat M)},
#' \eqn{S = E[\phi_i \phi_i^\top]}), and the inverse Hessian of the
#' integrated log-likelihood for one-step. See \code{vcov_corrected}
#' below for an alternative, more conservative BCA/BCM sandwich.
#'
#' @section Required inputs by method:
#' Only the chosen \code{method}s' inputs are validated; you do not have
#' to supply parameters for methods you are not fitting.
#'
#' \tabular{lll}{
#'   \strong{Method}    \tab \strong{Binary (K = 2)}                    \tab \strong{Multicategory (K > 2)} \cr
#'   \code{"naive"}    \tab \emph{none}                                 \tab \emph{none} \cr
#'   \code{"bca"}, \code{"bcm"}, \code{"cs"} \tab \code{c1} and \code{c2}, or \code{p01}/\code{p10}/\code{pi_z}, or \code{Pi} \tab \code{Pi} and \code{pi_z} \cr
#'   \code{"onestep"} (\code{fix_omega = FALSE}) \tab \emph{none} (mixture weights estimated) \tab \emph{none} \cr
#'   \code{"onestep"} (\code{fix_omega = TRUE})  \tab \code{p01}, \code{p10}, \code{pi_z} \tab \code{Pi}, \code{pi_z}
#' }
#'
#' \emph{Auto-derivations performed by} \code{mcglm}\emph{:}
#' \itemize{
#'   \item If \code{Pi} is supplied with \code{nrow(Pi) == 2},
#'     \code{p01 = Pi[2, 1]} and \code{p10 = Pi[1, 2]} are extracted.
#'   \item If \code{Pi} and \code{z_hat} are available but \code{pi_z}
#'     is not, \code{pi_z} is estimated by Bayesian inversion of the
#'     observed proxy frequencies (\code{.mcglm_estimate_pi_z}).
#'   \item If \code{p01}, \code{p10}, \code{pi_z} are available but
#'     \code{c1}, \code{c2} are not, the latter are computed from the
#'     formulas above.
#' }
#' Supplying \code{c1} and \code{c2} directly is an "advanced" path
#' that bypasses these checks. The formula interface
#' (\code{y ~ mc(z, Pi) + ...}) extracts \code{Pi} automatically and is
#' the recommended entry point.
#'
#' @param formula Either a formula of the form
#'   \code{y ~ mc(z, Pi) + x1 + x2}, where \code{\link{mc}()} marks the
#'   misclassified covariate \eqn{\hat Z_i} together with its
#'   misclassification matrix \eqn{\Pi}; or, when using the matrix
#'   interface, a numeric response vector of length \eqn{n}.
#' @param data A data frame (required for the formula interface). The
#'   variable named in \code{mc()} should be a \code{\link{factor}} or
#'   integer-valued vector with levels coded as \eqn{0, 1, \dots, K-1}.
#' @param family A \code{\link[stats]{family}} object or one of
#'   the strings \code{"poisson"}, \code{"binomial"}, \code{"gaussian"}
#'   (any \eqn{K}, all five methods), or \code{"multinomial"}
#'   (\code{naive} and \code{onestep} only).
#' @param method Character vector of estimators to fit. Any subset of
#'   \code{c("naive", "bca", "bcm", "cs", "cs_akn", "onestep")}; the
#'   default is the four analytical estimators
#'   \code{c("naive", "bca", "bcm", "cs")}. The \code{"cs_akn"} entry
#'   selects the Akazawa--Kinukawa--Nakamura (1998) corrected-score
#'   construction (an alternative formulation of \code{"cs"} based on
#'   the unbiased-surrogate transform \eqn{x = Q^{-1}(u - p_0)}); it
#'   needs only \eqn{\Pi} (not \eqn{\pi_z}) and is unsupported for
#'   \code{family = "multinomial"}.
#' @param p01 False-positive rate \eqn{p_{01} = \Pr(\hat Z = 1 \mid Z = 0)}
#'   (\code{K = 2} only). Auto-extracted from \code{Pi} when supplied.
#' @param p10 False-negative rate \eqn{p_{10} = \Pr(\hat Z = 0 \mid Z = 1)}
#'   (\code{K = 2} only). Auto-extracted from \code{Pi} when supplied.
#' @param pi_z Latent prevalence: a scalar
#'   \eqn{\pi = \Pr(Z = 1)} (\code{K = 2}) or a \eqn{K}-vector
#'   \eqn{(\pi_0, \dots, \pi_{K-1})} of class probabilities. If
#'   \code{NULL}, \code{pi_z} is estimated from \code{z_hat} and
#'   \code{Pi} via \eqn{\hat\pi = \Pi^{-1} \hat\pi_{\text{obs}}}
#'   (clamped to \eqn{[0.01, 0.99]}).
#' @param Pi The \eqn{K \times K} misclassification matrix
#'   \eqn{\Pi_{j\ell} = \Pr(\hat Z = j - 1 \mid Z = \ell - 1)}; columns
#'   must sum to 1. Extracted automatically from the \code{mc()} term
#'   when using the formula interface.
#' @param K Number of latent categories. Auto-detected from \code{Pi}
#'   or, failing that, from the levels of \code{z_hat}.
#' @param eta Optional regressed misclassification/prevalence parameters
#'   \eqn{\hat\eta = (\hat\eta_\rho^\top, \hat\eta_{01}^\top, \hat\eta_{10}^\top)^\top}
#'   for binary misclassification. Either a named list with components
#'   \code{rho}, \code{p01}, \code{p10} (each a coefficient vector or fitted
#'   \code{glm} object), or a numeric vector of length \eqn{3 \times \text{ncol}(w)}.
#' @param w Optional covariate matrix (or one-sided formula) for the
#'   misclassification and prevalence models \eqn{\expit(w^\top \eta)}.
#'   Defaults to the regression covariates \eqn{x}.
#' @param c1 Pre-computed misclassification constant
#'   \eqn{c_1 = p_{01}(1 - \pi)} (binary case, advanced use). Computed
#'   from \code{p01} and \code{pi_z} when \code{NULL}.
#' @param c2 Pre-computed misclassification constant
#'   \eqn{c_2 = p_{01}(1 - \pi) - p_{10} \pi} (binary case, advanced
#'   use). Computed from \code{p01}, \code{p10} and \code{pi_z} when
#'   \code{NULL}.
#' @param iterate Logical. If \code{TRUE}, iterate the BCA/BCM updates
#'   to convergence (the iterated BCM solves the corrected score
#'   equation; iterated BCA performs Newton steps on the additive
#'   correction).
#' @param jacobian Character. Either \code{"analytical"} (default) or
#'   \code{"numerical"}. Controls how the drift Jacobian
#'   \eqn{M = \partial m / \partial \psi} is evaluated for the
#'   multicategory (\eqn{K \ge 3}) BCM and CS estimators. The analytical
#'   form is exact and faster; the numerical form (forward differences)
#'   is provided as a fallback for cross-checking. The binary case
#'   (\eqn{K = 2}) always uses the analytical Jacobian.
#' @param fix_omega Logical. If \code{TRUE}, fix the mixture weights
#'   in the one-step estimator at the values implied by the supplied
#'   misclassification parameters; if \code{FALSE} (default), they are
#'   estimated jointly with \eqn{\psi} via a softmax parameterization.
#' @param vcov_corrected Logical. If \code{TRUE}, the BCA/BCM variance
#'   accounts for the additional uncertainty due to plug-in estimation
#'   of \eqn{\hat m(\hat\psi)} (joint score-and-drift sandwich); if
#'   \code{FALSE} (default), the drifting-regime asymptotic variance
#'   \eqn{A^{-1} C A^{-1}} is used. Both share the same point estimate.
#' @param weights Optional positive frequency weights of length \eqn{n}.
#' @param validation Optional list describing the sample used to estimate the
#'   probabilities for \code{method = "cs"}. For independent external
#'   validation use \code{list(z = true_codes, z_hat = proxy_codes)}. For an
#'   internal simple random subsample use
#'   \code{list(z = true_codes, index = regression_row_numbers)}; proxy codes
#'   are taken from those rows, which remain in the regression. Codes must be
#'   numeric integers in \code{0, ..., K-1}, in the model's category order.
#'   Supply \code{pi_z} and \code{Pi} (or binary \code{p01/p10}, or
#'   \code{c1/c2}) computed from these same validation observations using
#'   empirical proportions. This argument changes only the CS covariance,
#'   adding validation uncertainty and, for internal validation, the overlap
#'   covariance. Requires unweighted observations (or all weights equal to
#'   one) and every true category present in validation. With \code{NULL},
#'   the existing variance conditional on the supplied probabilities is used,
#'   including when \code{pi_z} is inferred from proxy frequencies.
#' @param J Number of response categories (multinomial family only;
#'   auto-detected).
#' @param homoskedastic Logical. For one-step Gaussian fits, assume a
#'   single residual variance across mixture components; if
#'   \code{FALSE}, binary fits estimate one variance per latent class,
#'   while \eqn{K > 2} fits estimate one baseline variance and one
#'   pooled nonbaseline variance.
#' @param optim_control List of control parameters passed to
#'   \code{\link[stats]{nlminb}} for the one-step estimator.
#' @param z_hat Integer vector of observed proxy values
#'   \eqn{\hat Z_i \in \{0, \dots, K-1\}} (matrix interface only).
#' @param x Covariate matrix \eqn{[x_1, \dots, x_n]^\top} of dimension
#'   \eqn{n \times r}; should include an intercept column when one is
#'   desired (matrix interface only).
#'
#' @return An object of class \code{"mcglm"} with components:
#'   \describe{
#'     \item{coefficients}{Named list of coefficient vectors
#'       \eqn{\hat\psi}, one per method.}
#'     \item{vcov}{Named list of asymptotic variance-covariance
#'       matrices, one per method (sandwich estimators from the
#'       paper's theorems).}
#'     \item{se}{Named list of standard errors per method.}
#'     \item{naive_fit}{The \code{\link[stats]{glm}} object from the
#'       naive fit (\code{NULL} for multinomial).}
#'     \item{method, family, K, n, p}{Metadata.}
#'     \item{weights}{Frequency weights used (\code{NULL} if
#'       unweighted).}
#'     \item{loglik_onestep, npar_onestep, vcov_onestep}{One-step
#'       log-likelihood, full optimized parameter count, and variance
#'       (when \code{onestep} is fit).}
#'     \item{call, formula}{The matched call and the formula
#'       (\code{NULL} if the matrix interface was used).}
#'   }
#'   See \code{\link{summary.mcglm}}, \code{\link{vcov.mcglm}},
#'   \code{\link{confint.mcglm}}, \code{\link{predict.mcglm}}, and
#'   \code{\link{logLik.mcglm}} for downstream methods.
#'
#' @references
#' Nakamura, T. (1990). Corrected score function for errors-in-variables
#' models: methodology and application to generalized linear models.
#' \emph{Biometrika}, 77(1), 127--137. \doi{10.1093/biomet/77.1.127}.
#' Originating reference for the corrected-score concept.
#'
#' Akazawa, K., Kinukawa, N. and Nakamura, T. (1998). A note on the
#' corrected score function adjusting for misclassification.
#' \emph{Journal of the Japan Statistical Society}, 28(1), 115--123.
#' Establishes existence of a corrected score for misclassified discrete
#' covariates and gives explicit forms for normal, Poisson, and binary
#' logistic GLMs; foundational reference for the \code{"cs"} method.
#'
#' Battaglia, L., Christensen, T., Hansen, S. and Sacher, S. (2025).
#' Inference for regression with variables generated by AI or machine
#' learning. \emph{arXiv preprint arXiv:2402.15585}. Drifting-asymptotics
#' framework for linear models with AI/ML-generated covariates;
#' originating reference for the BCA/BCM and one-step estimators.
#'
#' Warmuz, M. and Beresewicz, M. (2026). Bias corrections for GLMs with
#' misclassified covariates. Working paper. Extends the drifting-regime
#' BCA/BCM to GLMs, derives the corrected-score estimator under fixed
#' misclassification probabilities for binary and multicategory latent
#' regressors with explicit drift formulas, and gives the plug-in
#' robustness conditions used here when \eqn{(\Pi, \pi)} are estimated
#' from auxiliary data.
#'
#' @seealso \code{\link{simex}} for SIMEX and MC-SIMEX corrections;
#'   \code{\link{mc}} for the formula-interface marker;
#'   \code{\link{summary.mcglm}}, \code{\link{vcov.mcglm}},
#'   \code{\link{confint.mcglm}}.
#'
#' @examples
#' # --- Formula interface (binary misclassification, Poisson) ---
#' set.seed(42)
#' n  <- 2000
#' Pi <- matrix(c(0.90, 0.10,    # Pi[2,1] = p01 = 0.10
#'                0.15, 0.85),   # Pi[1,2] = p10 = 0.15
#'              nrow = 2, byrow = FALSE)
#' z      <- rbinom(n, 1, 0.4)               # latent Z, prevalence pi = 0.4
#' z_hat  <- z
#' z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)        # flip 0 -> 1
#' z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)    # flip 1 -> 0
#' x1 <- rnorm(n)
#' y  <- rpois(n, exp(0.8 * z + -0.5 + 0.7 * x1))
#' df <- data.frame(y = y, z = factor(z_hat), x1 = x1)
#'
#' fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
#'              method = c("naive", "bca", "bcm", "cs"))
#' fit
#' summary(fit)
#' confint(fit, method = "cs")
#'
#' # Same problem with the AKN98 unbiased-surrogate corrected score
#' # (needs only Pi, not pi_z). The two CS variants give population-
#' # equivalent estimators and agree closely in large samples.
#' fit_akn <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
#'                  method = c("cs", "cs_akn"))
#' coef(fit_akn, method = "cs_akn")
#'
#' # --- Matrix interface, supplying p01/p10/pi_z directly ---
#' x_mat <- cbind(1, x1)
#' fit2 <- mcglm(y, z_hat = as.integer(z_hat), x = x_mat,
#'               family = "poisson", method = c("naive", "cs"),
#'               p01 = 0.10, p10 = 0.15, pi_z = 0.4)
#'
#' # --- Multicategory misclassification (K = 3) ---
#' \donttest{
#' set.seed(3)
#' n   <- 1000
#' Pi3 <- matrix(c(0.8, 0.1, 0.1,
#'                 0.1, 0.8, 0.1,
#'                 0.1, 0.1, 0.8), 3, 3, byrow = TRUE)
#' pi3 <- c(0.5, 0.3, 0.2)
#' Z   <- sample(0:2, n, replace = TRUE, prob = pi3)
#' Zh  <- vapply(Z, function(z) sample(0:2, 1, prob = Pi3[, z + 1]), 0L)
#' x1  <- rnorm(n)
#' y   <- rpois(n, exp(c(0, 1.0, -0.9)[Z + 1] + 0.8 - 0.7 * x1))
#' fit3 <- mcglm(y, z_hat = Zh, x = cbind(1, x1), family = "poisson",
#'               method = c("naive", "cs"), Pi = Pi3, pi_z = pi3)
#' coef(fit3, method = "cs")
#' }
#'
#' @export
mcglm <- function(formula, data = NULL, family = "poisson",
                  method = c("naive", "bca", "bcm", "cs"),
                  p01 = NULL, p10 = NULL, pi_z = NULL,
                  Pi = NULL, K = NULL,
                  c1 = NULL, c2 = NULL, eta = NULL, w = NULL,
                  iterate = FALSE,
                  jacobian = c("analytical", "numerical"),
                  fix_omega = FALSE,
                  vcov_corrected = FALSE,
                  weights = NULL,
                  J = NULL,
                  homoskedastic = TRUE,
                  optim_control = list(),
                  z_hat = NULL, x = NULL, validation = NULL) {

  cl <- match.call()
  jacobian <- match.arg(jacobian)
  formula_obj <- NULL

  # Process w if supplied as formula
  if (!is.null(w)) {
    if (inherits(w, "formula")) {
      w <- stats::model.matrix(w, data)
    } else if (is.character(w)) {
      w <- stats::model.matrix(stats::as.formula(paste("~", paste(w, collapse = " + "))), data)
    } else if (is.data.frame(w)) {
      w <- as.matrix(w)
    }
  }

  # --- Dispatch: formula vs matrix interface ---
  if (inherits(formula, "formula")) {
    formula_obj <- formula
    caller_env <- parent.frame()
    # Offsets are never plumbed into the corrected estimators; reject
    # rather than silently fit without them.
    if (!is.null(attr(terms(formula, data = data), "offset")))
      stop("offset() terms are not supported by mcglm(); the corrected ",
           "estimators would ignore the offset and give biased estimates.",
           call. = FALSE)
    parsed <- .mcglm_parse_formula(formula, data, caller_env)
    y     <- parsed$y
    z_hat <- parsed$z_hat
    x     <- parsed$x
    if (is.null(Pi) && !is.null(parsed$Pi)) Pi <- parsed$Pi
    if (is.null(K))  K  <- parsed$K
    z_levels <- parsed$z_levels
    x_levels <- parsed$x_levels
    x_names  <- colnames(parsed$x)
  } else {
    # Matrix interface: first arg is y. Coefficients keep the positional
    # gamma/alpha names; no factor-level bookkeeping applies.
    y <- formula
    z_levels <- NULL
    x_levels <- NULL
    x_names  <- NULL
  }

  if (!is.null(Pi)) Pi <- .validate_mc_matrix(Pi, "Pi")
  .check_probability_scalar <- function(value, name) {
    if (is.null(value)) return(invisible(NULL))
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 1)
      stop(name, " must be one finite probability in [0, 1].", call. = FALSE)
    invisible(NULL)
  }
  .check_probability_scalar(p01, "p01")
  .check_probability_scalar(p10, "p10")
  if (!is.null(pi_z) &&
      (!is.numeric(pi_z) || !length(pi_z) || any(!is.finite(pi_z)) ||
       any(pi_z < 0 | pi_z > 1)))
    stop("pi_z must contain finite probabilities in [0, 1].", call. = FALSE)
  for (nm in c("c1", "c2")) {
    value <- get(nm)
    if (!is.null(value) &&
        (!is.numeric(value) || any(!is.finite(value)) || (length(value) != 1L && !is.null(data) && length(value) != nrow(data))))
      stop(nm, " must be one finite numeric value or a vector of length n.", call. = FALSE)
  }

  # --- Derive p01/p10 from 2x2 Pi ---
  if (!is.null(Pi) && nrow(Pi) == 2L) {
    if (is.null(p01)) p01 <- Pi[2, 1]
    if (is.null(p10)) p10 <- Pi[1, 2]
  }

  # --- Estimate pi_z from observed data + Pi if not supplied ---
  if (is.null(pi_z) && !is.null(Pi) && !is.null(z_hat)) {
    .mcglm_check_z_hat(z_hat)
    pi_z <- .mcglm_estimate_pi_z(as.integer(z_hat), Pi)
  }

  # --- Derive c1/c2 from (p01, p10, pi_z) if not supplied directly ---
  if (is.null(c1) && !is.null(p01) && !is.null(pi_z))
    c1 <- p01 * (1 - pi_z)
  if (is.null(c2) && !is.null(p01) && !is.null(p10) && !is.null(pi_z))
    c2 <- p01 * (1 - pi_z) - p10 * pi_z

  # --- Core fitting ---
  out <- .mcglm_fit(y = y, z_hat = z_hat, x = x, family = family,
                    method = method, p01 = p01, p10 = p10, pi_z = pi_z,
                    Pi = Pi, K = K, c1 = c1, c2 = c2, eta = eta, w = w, iterate = iterate,
                    jacobian = jacobian,
                    fix_omega = fix_omega,
                    vcov_corrected = vcov_corrected,
                    weights = weights, J = J,
                    homoskedastic = homoskedastic,
                    optim_control = optim_control,
                    z_levels = z_levels, x_names = x_names,
                    validation = validation)
  out$call     <- cl
  out$formula  <- formula_obj
  out$z_levels <- z_levels
  out$x_levels <- x_levels
  out
}


#' Validate z_hat coding for the mcglm matrix interface
#' @keywords internal
.mcglm_check_z_hat <- function(z_hat) {
  if (is.factor(z_hat))
    stop("z_hat must be an integer vector coded 0, ..., K-1, not a factor ",
         "(as.integer(factor) is 1-based and would distort the estimates).",
         call. = FALSE)
  if (any(!is.finite(z_hat)) || any(z_hat != round(z_hat)))
    stop("z_hat must contain finite integer codes 0, ..., K-1.",
         call. = FALSE)
  invisible(TRUE)
}

#' Parse mcglm formula extracting mc() term
#' @keywords internal
.mcglm_parse_formula <- function(formula, data, env, require_y = TRUE,
                                 z_levels = NULL, x_levels = NULL) {
  if (is.null(data))
    stop("'data' argument is required when using the formula interface.")

  lhs <- formula[[2]]
  # Extract response (optional for prediction on new data)
  y_name <- deparse(lhs)
  y <- data[[y_name]]
  if (is.null(y) && require_y)
    stop("Response variable '", y_name, "' not found in data.")

  # Let R's formula machinery identify mc() and preserve ordinary formula
  # semantics (intercept removal, transformations, interactions, contrasts).
  tt <- terms(formula, specials = "mc", data = data)
  special_idx <- attr(tt, "specials")$mc
  if (length(special_idx) != 1L)
    stop("Formula must contain exactly one mc() term, e.g., y ~ mc(z, Pi) + x")

  variables <- as.list(attr(tt, "variables"))[-1L]
  mc_call <- variables[[special_idx]]
  if (!is.call(mc_call) || !identical(mc_call[[1L]], as.name("mc")) ||
      length(mc_call) < 2L || !is.symbol(mc_call[[2L]]))
    stop("mc() must contain a bare variable name, e.g. mc(z, Pi).",
         call. = FALSE)

  factors <- attr(tt, "factors")
  mc_term_idx <- which(factors[special_idx, , drop = TRUE] != 0)
  if (length(mc_term_idx) != 1L || attr(tt, "order")[mc_term_idx] != 1L)
    stop("mc() must appear as one bare main-effect term; interactions and ",
         "transformations involving mc() are not supported.", call. = FALSE)

  var_name <- as.character(mc_call[[2L]])
  # Pi may instead be supplied to mcglm()/a wrapper, and prediction does not
  # need it. Preserve the existing behavior that an unresolvable expression
  # yields NULL and is validated later by the fitting caller when required.
  mat_val <- if (length(mc_call) >= 3L) {
    tryCatch(eval(mc_call[[3L]], data, env), error = function(e) NULL)
  } else {
    NULL
  }
  mc_info <- list(variable = var_name,
                  Pi = if (is.null(mat_val)) NULL else as.matrix(mat_val))

  # Remove the mc() term without reconstructing the formula from strings.
  # This retains the original intercept, contrast, and environment metadata.
  x_terms <- stats::drop.terms(tt, mc_term_idx, keep.response = FALSE)

  # Get z_hat from data
  z_var <- data[[mc_info$variable]]
  if (is.null(z_var))
    stop("mc() variable '", mc_info$variable, "' not found in data.")

  Pi <- mc_info$Pi

  # Convert z to 0-based integer. When reference levels are supplied
  # (prediction on new data), recode AGAINST THEM: a factor in new data
  # carries its own level set/order, and coding by that would silently
  # mis-map categories.
  if (is.factor(z_var) || is.character(z_var)) {
    if (is.null(z_levels)) {
      if (is.character(z_var)) z_var <- factor(z_var)
      z_levels <- levels(z_var)
      z_hat <- as.integer(z_var) - 1L
    } else {
      idx <- match(as.character(z_var), z_levels)
      if (anyNA(idx))
        stop("mc() variable '", mc_info$variable, "' has level(s) ",
             paste0("'", unique(as.character(z_var)[is.na(idx)]), "'",
                    collapse = ", "),
             " not present in the training levels (",
             paste(z_levels, collapse = ", "), ").", call. = FALSE)
      z_hat <- idx - 1L
    }
  } else {
    if (!is.null(z_levels))
      z_levels <- NULL  # integer coding: levels do not apply
    .mcglm_check_z_hat(z_var)
    z_hat <- as.integer(z_var)
  }

  # Determine K (reference levels, when supplied, fix K regardless of
  # which categories happen to appear in this particular data set)
  K <- if (!is.null(Pi)) nrow(Pi)
       else if (!is.null(z_levels)) length(z_levels)
       else length(unique(z_hat))

  # Validate Pi when provided
  if (!is.null(Pi)) {
    Pi <- .validate_mc_matrix(Pi, "Pi")
  }

  # Build x matrix from remaining terms; x_levels (when supplied) recode
  # factor covariates in new data against the training levels, exactly
  # like predict.lm()'s xlev mechanism
  mf <- model.frame(x_terms, data = data, xlev = x_levels,
                    na.action = na.pass)
  if (anyNA(mf))
    stop("Covariates contain missing values after applying the ",
         "training factor levels; check for NA values or factor ",
         "levels not seen during fitting.", call. = FALSE)
  x <- model.matrix(x_terms, data = mf)
  x_levels_out <- .getXlevels(attr(mf, "terms"), mf)
  if (is.null(x_levels_out)) x_levels_out <- list()

  list(y = if (is.null(y)) NULL else as.numeric(y),
       z_hat = z_hat, x = x, Pi = Pi, K = K,
       z_levels = z_levels,
       x_levels = if (is.null(x_levels)) x_levels_out else x_levels)
}


#' Core mcglm fitting (internal, called by mcglm after dispatch)
#' @keywords internal
.mcglm_fit <- function(y, z_hat, x, family = "poisson",
                       method = c("naive", "bca", "bcm", "cs"),
                       p01 = NULL, p10 = NULL, pi_z = NULL,
                       Pi = NULL, K = NULL,
                       c1 = NULL, c2 = NULL, eta = NULL, w = NULL,
                       iterate = FALSE,
                       jacobian = c("analytical", "numerical"),
                       fix_omega = FALSE,
                       vcov_corrected = FALSE,
                       weights = NULL,
                       J = NULL,
                       homoskedastic = TRUE,
                       optim_control = list(),
                       z_levels = NULL,
                       x_names = NULL, validation = NULL) {

  # --- input validation ---
  y     <- as.numeric(y)
  .mcglm_check_z_hat(z_hat)
  z_hat <- as.integer(z_hat)
  x     <- as.matrix(x)
  n     <- length(y)
  jacobian <- match.arg(jacobian)
  method <- match.arg(method,
                      c("naive", "bca", "bcm", "cs", "cs_akn", "onestep"),
                      several.ok = TRUE)

  stopifnot(nrow(x) == n, length(z_hat) == n)

  # --- validate frequency weights ---
  wt <- NULL
  if (!is.null(weights)) {
    wt <- as.numeric(weights)
    if (length(wt) != n)
      stop("weights must have length n (", n, "), got ", length(wt))
    if (any(!is.finite(wt)) || any(wt <= 0))
      stop("weights must contain only finite positive values")
  }

  # --- multinomial family ---
  is_multinomial <- is.character(family) && family == "multinomial"
  if (is_multinomial) {
    unsupported <- setdiff(method, c("naive", "onestep"))
    if (length(unsupported) > 0)
      stop("For multinomial family, only 'naive' and 'onestep' methods are ",
           "supported. Unsupported: ", paste(unsupported, collapse = ", "))
    if ("cs_akn" %in% method)
      stop("'cs_akn' is not supported for multinomial responses (AKN98 ",
           "covers univariate-response GLMs only).")
    y <- as.integer(y)
    if (is.null(J)) J <- length(unique(y))
    if (any(y < 0L) || any(y >= J))
      stop("For multinomial, y must be integers in {0, ..., J-1}")
  } else {
    family <- .normalize_family(family)
  }

  # Determine binary vs multicategory. Prefer Pi: length(unique(z_hat))
  # undercounts K when a category is unobserved in the sample.
  if (is.null(K)) {
    if (!is.null(Pi)) {
      K <- nrow(Pi)
    } else {
      unique_z <- sort(unique(z_hat))
      K <- if (all(unique_z %in% c(0L, 1L))) 2L else length(unique_z)
    }
  }
  is_binary <- (K == 2L)

  if (!is.null(pi_z)) {
    expected <- if (is_binary) 1L else K
    if (length(pi_z) != expected)
      stop("pi_z must have length ", expected, " for K = ", K, ".",
           call. = FALSE)
    if (!is_binary && abs(sum(pi_z) - 1) > 1e-6)
      stop("For K > 2, pi_z must sum to 1.", call. = FALSE)
  }

  if (any(z_hat < 0L) || any(z_hat >= K))
    stop("z_hat must be coded 0, ..., K-1 (K = ", K, "); got values in [",
         min(z_hat), ", ", max(z_hat), "]. 1-based coding silently ",
         "distorts binary estimates and indexes out of bounds for K > 2.",
         call. = FALSE)

  # A category declared in Pi/K but absent from z_hat makes its gamma
  # non-identifiable (empty dummy column -> NA), and the drift terms of
  # bca/bcm/cs depend on that gamma, so the whole correction is
  # infeasible. Fail here with the actual cause instead of the
  # downstream singular-matrix / non-finite errors.
  z_counts <- tabulate(z_hat + 1L, nbins = K)
  if (any(z_counts == 0L)) {
    empty <- which(z_counts == 0L) - 1L
    lab <- if (!is.null(z_levels) && length(z_levels) == K) {
      paste0("'", z_levels[empty + 1L], "' (code ", empty, ")")
    } else {
      paste0("code ", empty)
    }
    stop("Category ", paste(lab, collapse = ", "), " of z_hat has no ",
         "observations, but K = ", K, " (from Pi or the K argument). ",
         "Its gamma coefficient is not identifiable and the ",
         "misclassification correction depends on it. Either collapse ",
         "Pi to the observed categories (drop the row and column, then ",
         "renormalize the columns to sum to 1) or check the factor ",
         "level coding of z_hat.", call. = FALSE)
  }

  # Check c1 and c2
  for (nm in c("c1", "c2")) {
    value <- get(nm)
    if (!is.null(value) &&
        (!is.numeric(value) || any(!is.finite(value)) || (length(value) != 1L && length(value) != n)))
      stop(nm, " must be one finite numeric value or a vector of length n.", call. = FALSE)
  }

  # If eta is supplied, compute c1 and c2 from eta and w
  if (!is.null(eta)) {
    w_use <- if (!is.null(w)) as.matrix(w) else x
    c_pack <- .mcglm_compute_c_bin(eta, w_use)
    c1   <- c_pack$c1
    c2   <- c_pack$c2
    p01  <- c_pack$p01
    p10  <- c_pack$p10
    pi_z <- c_pack$rho
    eta  <- c_pack$eta
    w    <- w_use
  }

  # Validation affects CS covariance (and estimates eta if needed)
  if (!is.null(validation)) {
    if (!any(c("cs", "cs_akn") %in% method))
      stop("validation requires method = 'cs' or 'cs_akn'.", call. = FALSE)
    validation <- .mcglm_prepare_cs_validation(validation, z_hat, K, wt,
                                                c1, c2, Pi, pi_z,
                                                x = x, w = w, eta = eta)
    if (isTRUE(validation$regressed)) {
      c1   <- validation$c1
      c2   <- validation$c2
      p01  <- validation$p01
      p10  <- validation$p10
      pi_z <- validation$rho
      eta  <- validation$eta
      w    <- validation$w
    }
  }

  # --- validate misclassification parameters ---
  needs_correction <- any(method %in% c("bca", "bcm", "cs", "cs_akn", "onestep"))
  needs_cs_akn     <- "cs_akn" %in% method
  if (needs_correction) {
    if (is_binary) {
      needs_bca_bcm_cs <- any(method %in% c("bca", "bcm", "cs"))
      needs_onestep    <- "onestep" %in% method
      if (needs_bca_bcm_cs && (is.null(c1) || is.null(c2)))
        stop("For binary bca/bcm/cs, supply c1 and c2 (or Pi, or p01/p10/pi_z).")
      if (needs_onestep && fix_omega && (is.null(p01) || is.null(p10) || is.null(pi_z)))
        stop("For onestep with fix_omega=TRUE, supply p01, p10, and pi_z.")
      if (needs_cs_akn) {
        if (is.null(Pi)) {
          if (is.null(p01) || is.null(p10))
            stop("For binary cs_akn, supply Pi or both p01 and p10.")
          if (length(p01) == 1L && length(p10) == 1L)
            Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2L, 2L)
        }
      }
    } else {
      needs_bca_bcm_cs <- any(method %in% c("bca", "bcm", "cs"))
      needs_onestep    <- "onestep" %in% method
      if (needs_bca_bcm_cs) {
        if (is.null(Pi) || is.null(pi_z))
          stop("For multicategory bca/bcm/cs, supply Pi and pi_z.")
        Pi   <- as.matrix(Pi)
        pi_z <- as.numeric(pi_z)
        stopifnot(nrow(Pi) == K, ncol(Pi) == K, length(pi_z) == K)
      }
      if (needs_onestep && fix_omega && (is.null(Pi) || is.null(pi_z)))
        stop("For multicategory onestep with fix_omega=TRUE, supply Pi and pi_z.")
      if (needs_cs_akn) {
        if (is.null(Pi))
          stop("For multicategory cs_akn, supply Pi.")
        Pi <- as.matrix(Pi)
        stopifnot(nrow(Pi) == K, ncol(Pi) == K)
      }
    }
  }


  # --- build xi_hat once ---
  if (!is_multinomial) {
    xi_hat <- .mcglm_build_xi_hat(z_hat, x, K)
  }

  # --- AKN98 unbiased surrogate (built once, reused for fit + variance) ---
  x_akn <- NULL
  if (needs_cs_akn) {
    if (is_binary && (!is.null(p01) && !is.null(p10))) {
      Delta <- 1 - p01 - p10
      if (any(abs(Delta) < 0.001))
        warning("cs_akn: 1 - p01 - p10 is near zero (< 0.001); clamping applied for numerical stability.", call. = FALSE)
      Delta <- ifelse(abs(Delta) < 0.001, ifelse(Delta >= 0, 0.001, -0.001), Delta)
      x_akn <- as.matrix((z_hat - p01) / Delta)
    } else {
      Q_pack <- .akn_build_Q(Pi, K)
      x_akn  <- .akn_unbiased_surrogate(z_hat, Q_pack$Q_inv, Q_pack$p0, K)
    }
  }

  # --- fit ---
  results <- list()
  onestep_vcov   <- NULL
  onestep_loglik <- NULL
  onestep_npar   <- NULL
  naive_glm_fit  <- NULL

  if (is_multinomial) {
    naive_mn <- .mcglm_fit_naive_multinomial(y, z_hat, x, J, K, wt = wt)
    psi      <- naive_mn$coefficients
    results$naive <- psi
    p_total  <- length(psi)

    if ("onestep" %in% method) {
      os <- .mcglm_fit_onestep_multinomial(y, z_hat, x, J, K,
                                           p01 = p01, p10 = p10, pi_z = pi_z,
                                           Pi = Pi, fix_omega = fix_omega,
                                           optim_control = optim_control,
                                           wt = wt)
      results$onestep <- os$coefficients
      onestep_vcov    <- os$vcov
      onestep_loglik  <- os$loglik
      onestep_npar    <- os$npar
    }

  } else if (is_binary) {
    naive <- .mcglm_fit_naive_bin(y, xi_hat, family, wt = wt)
    psi   <- naive$coefficients
    results$naive <- psi
    p_total <- length(psi)
    naive_glm_fit <- naive$glm_fit

    if ("bca" %in% method)
      results$bca <- .mcglm_fit_bca_bin(psi, y, xi_hat, x, family, c1, c2,
                                         iterate = iterate, wt = wt)
    if ("bcm" %in% method)
      results$bcm <- .mcglm_fit_bcm_bin(psi, y, xi_hat, x, family, c1, c2,
                                         iterate = iterate, wt = wt)
    if ("cs"  %in% method)
      results$cs  <- .mcglm_fit_cs_bin(psi, y, xi_hat, x, family, c1, c2,
                                        wt = wt)
    if ("cs_akn" %in% method)
      results$cs_akn <- .mcglm_fit_cs_akn(psi, y, x_akn, x, K, family, wt = wt)

    if ("onestep" %in% method) {
      os <- .mcglm_fit_onestep_bin(y, z_hat, x, family, p01, p10, pi_z,
                                   fix_omega = fix_omega,
                                   homoskedastic = homoskedastic,
                                   optim_control = optim_control, wt = wt)
      results$onestep <- os$coefficients
      onestep_vcov    <- os$vcov
      onestep_loglik  <- os$loglik
      onestep_npar    <- os$npar
    }

  } else {
    naive <- .mcglm_fit_naive_multi(y, xi_hat, family, wt = wt)
    psi   <- naive$coefficients
    results$naive <- psi
    p_total <- length(psi)
    naive_glm_fit <- naive$glm_fit

    if ("bca" %in% method)
      results$bca <- .mcglm_fit_bca_multi(psi, y, xi_hat, z_hat, x, K, family,
                                           Pi, pi_z, iterate = iterate, wt = wt)
    if ("bcm" %in% method)
      results$bcm <- .mcglm_fit_bcm_multi(psi, y, xi_hat, z_hat, x, K, family,
                                           Pi, pi_z, iterate = iterate, wt = wt,
                                           jacobian = jacobian)
    if ("cs"  %in% method)
      results$cs  <- .mcglm_fit_cs_multi(psi, y, xi_hat, z_hat, x, K, family,
                                          Pi, pi_z, wt = wt,
                                          jacobian = jacobian)
    if ("cs_akn" %in% method)
      results$cs_akn <- .mcglm_fit_cs_akn(psi, y, x_akn, x, K, family, wt = wt)

    if ("onestep" %in% method) {
      os <- .mcglm_fit_onestep_multi(y, z_hat, x, K, family, Pi, pi_z,
                                     fix_omega = fix_omega,
                                     homoskedastic = homoskedastic,
                                     optim_control = optim_control, wt = wt)
      results$onestep <- os$coefficients
      onestep_vcov    <- os$vcov
      onestep_loglik  <- os$loglik
      onestep_npar    <- os$npar
    }
  }

  # --- parameter names ---
  if (is_multinomial) {
    s <- K - 1
    r <- ncol(x)
    block_size <- s + r
    nms <- character(p_total)
    for (jj in seq_len(J - 1)) {
      offset <- (jj - 1) * block_size
      if (s == 1) {
        nms[offset + 1] <- paste0("y", jj, ":gamma")
      } else if (s > 1) {
        nms[offset + seq_len(s)] <- paste0("y", jj, ":gamma", seq_len(s))
      }
      nms[offset + s + seq_len(r)] <- paste0("y", jj, ":alpha", seq_len(r) - 1)
    }
  } else {
    # Formula fits with factor inputs get informative names: gammas
    # labelled by the (non-baseline) factor levels of the mc() variable,
    # alphas by the model.matrix column names. Integer z / matrix
    # interface keeps the positional gamma/alpha names.
    alpha_nms <- if (!is.null(x_names) && length(x_names) == ncol(x)) {
      x_names
    } else {
      paste0("alpha", seq_len(ncol(x)) - 1)
    }
    gamma_nms <- if (!is_binary && !is.null(z_levels) &&
                     length(z_levels) == K) {
      paste0("gamma_", z_levels[-1])
    } else {
      if (is_binary) "gamma" else paste0("gamma", seq_len(K - 1))
    }
    nms <- c(gamma_nms, alpha_nms)
  }
  results <- lapply(results, function(v) {
    v <- as.numeric(v)
    names(v) <- nms
    v
  })

  # --- per-method asymptotic variance ---
  vcov_list <- list()
  if (!is_multinomial) {
    for (nm in names(results)) {
      psi_nm <- unname(results[[nm]])
      V <- tryCatch(
        if (nm == "naive") {
          if (is_binary)
            .mcglm_vcov_naive(psi_nm, y, xi_hat, family, wt = wt)
          else
            .mcglm_vcov_naive_multi(psi_nm, y, xi_hat, z_hat, x, K, family,
                                    wt = wt)
        } else if (nm %in% c("bca", "bcm")) {
          have_probs <- (!is.null(p01) && !is.null(p10) && !is.null(pi_z)) ||
                        (!is.null(c1) && !is.null(c2))
          if (is_binary) {
            .mcglm_vcov_bc_bin(psi_nm, y, xi_hat, x, family,
                               p01 = p01, p10 = p10, pi_z = pi_z,
                               c1 = c1, c2 = c2,
                               psi_naive = unname(results$naive),
                               type = nm,
                               corrected = vcov_corrected && have_probs,
                               wt = wt)
          } else {
            have_multi <- !is.null(Pi) && !is.null(pi_z)
            .mcglm_vcov_bc_multi(psi_nm, y, xi_hat, z_hat, x, K, family,
                                 Pi = Pi, pi_z = pi_z,
                                 psi_naive = unname(results$naive),
                                 type = nm,
                                 corrected = vcov_corrected && have_multi,
                                 wt = wt,
                                 jacobian = jacobian)
          }
        } else if (nm == "cs") {
          if (is_binary)
            .mcglm_vcov_cs_bin(psi_nm, y, xi_hat, x, family, p01, p10, pi_z,
                               c1 = c1, c2 = c2,
                               wt = wt, validation = validation)
          else
            .mcglm_vcov_cs_multi(psi_nm, y, xi_hat, z_hat, x, K, family,
                                 Pi, pi_z, wt = wt,
                                 jacobian = jacobian, validation = validation)
        } else if (nm == "cs_akn") {
          .mcglm_vcov_cs_akn(psi_nm, y, x_akn, x, K, family, wt = wt,
                             validation = validation, w = w,
                             p01 = p01, p10 = p10, z_hat = z_hat)
        } else if (nm == "onestep") {
          onestep_vcov
        },
        error = function(e) {
          warning("Variance computation failed for method '", nm, "': ",
                  conditionMessage(e))
          NULL
        }
      )
      if (!is.null(V)) {
        dimnames(V) <- list(nms, nms)
      }
      vcov_list[[nm]] <- V
    }
  } else if (!is.null(onestep_vcov)) {
    vcov_list$onestep <- onestep_vcov
  }

  se_list <- lapply(vcov_list, function(V) {
    if (is.null(V)) return(NULL)
    s <- sqrt(pmax(diag(V), 0))
    names(s) <- rownames(V)
    s
  })

  out <- list(
    coefficients = results,
    vcov         = vcov_list,
    se           = se_list,
    naive_fit    = naive_glm_fit,
    method       = method,
    family       = family,
    K            = K,
    n            = n,
    p            = p_total,
    weights      = wt,
    y            = y,
    z_hat        = z_hat,
    x            = x,
    is_multinomial = is_multinomial,
    misclass     = list(c1 = c1, c2 = c2, p01 = p01, p10 = p10, pi_z = pi_z, Pi = Pi, eta = eta, w = w)
  )
  if (!is_multinomial) out$xi_hat <- xi_hat
  if (!is.null(validation)) out$validation <- validation
  if (is_multinomial)  out$J     <- J

  if (!is.null(onestep_vcov)) {
    out$vcov_onestep   <- onestep_vcov  # kept for backwards compatibility
    out$loglik_onestep <- onestep_loglik
    out$npar_onestep   <- onestep_npar
  }

  structure(out, class = "mcglm")
}


#' Default method for an \code{mcglm} object
#' @keywords internal
.mcglm_default_method <- function(object) {
  if (!is.null(object$method) && length(object$method))
    return(object$method[length(object$method)])
  names(object$coefficients)[1]
}

.mcglm_family_string <- function(object) {
  if (is.character(object$family)) return(object$family)
  if (inherits(object$family, "family")) return(object$family$family)
  as.character(object$family)
}


#' Print method for \code{mcglm} fits
#'
#' Mirrors the layout of \code{print.glm}: shows the call, the family,
#' a coefficient table (one column per estimation method), the residual
#' degrees of freedom and (when available) the deviance and AIC of the
#' naive GLM fit.
#' @param x An object of class \code{"mcglm"}.
#' @param digits Number of digits to display.
#' @param ... Unused.
#' @export
print.mcglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  if (!is.null(x$call)) print(x$call) else cat("mcglm(...)\n")
  cat("\nFamily: ", .mcglm_family_string(x),
      "  |  n = ", x$n, ", K = ", x$K,
      ", p = ", x$p, "\n", sep = "")
  cat("Methods: ", paste(x$method, collapse = ", "), "\n\n", sep = "")

  cat("Coefficients:\n")
  coefs <- do.call(cbind, x$coefficients)
  colnames(coefs) <- toupper(names(x$coefficients))
  print.default(format(coefs, digits = digits), print.gap = 2L,
                quote = FALSE)

  if (!is.null(x$naive_fit)) {
    rdf <- x$naive_fit$df.residual
    cat("\nDegrees of Freedom:", x$naive_fit$df.null, "Total (i.e. Null); ",
        rdf, "Residual\n")
    cat("Null Deviance:    ",
        format(signif(x$naive_fit$null.deviance, digits)), "\n")
    cat("Residual Deviance:",
        format(signif(x$naive_fit$deviance, digits)),
        " | AIC (naive):", format(signif(stats::AIC(x$naive_fit), digits)),
        "\n")
  }
  if (!is.null(x$loglik_onestep)) {
    cat("logLik (onestep): ", format(signif(x$loglik_onestep, digits)), "\n")
  }
  invisible(x)
}


#' Extract coefficients from an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Optional character: which estimator to extract
#'   (\code{"naive"}, \code{"bca"}, \code{"bcm"}, \code{"cs"},
#'   \code{"onestep"}). If \code{NULL}, returns the coefficients of the
#'   default method (the last one fit), consistently with
#'   \code{\link{vcov.mcglm}} and \code{\link{predict.mcglm}} -- and as
#'   required by downstream tools such as \code{lmtest::coeftest} and
#'   \pkg{marginaleffects}. Use \code{method = "all"} for the named list
#'   of all fitted coefficient vectors (the pre-0.6.0 default).
#' @param ... Unused.
#' @export
coef.mcglm <- function(object, method = NULL, ...) {
  if (identical(method, "all")) return(object$coefficients)
  if (is.null(method)) method <- .mcglm_default_method(object)
  if (!method %in% names(object$coefficients))
    stop("Method '", method, "' was not fit. Available: ",
         paste(names(object$coefficients), collapse = ", "))
  object$coefficients[[method]]
}


#' Variance-covariance matrix of an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Optional character: which estimator to use. Defaults to the
#'   last method fit (e.g. \code{"onestep"} when present, else the most
#'   refined correction available).
#' @param ... Unused.
#' @return A \eqn{p \times p} sandwich variance-covariance matrix.
#' @export
vcov.mcglm <- function(object, method = NULL, ...) {
  vcov_list <- object[["vcov", exact = TRUE]]
  if (is.null(vcov_list)) {
    # backwards-compat (older fits stored only vcov_onestep)
    return(object[["vcov_onestep", exact = TRUE]])
  }
  if (is.null(method)) method <- .mcglm_default_method(object)
  if (!method %in% names(vcov_list))
    stop("No variance available for method '", method, "'.")
  vcov_list[[method]]
}


#' Summary of an \code{mcglm} fit
#'
#' Returns a Wald-style coefficient table (estimate, standard error,
#' z-value, p-value) for each fitted method, using the asymptotic
#' sandwich variance derived in Battaglia et al. (2025).
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
summary.mcglm <- function(object, ...) {
  fam_str <- .mcglm_family_string(object)
  ans <- list(
    call    = object$call,
    family  = fam_str,
    n       = object$n,
    p       = object$p,
    K       = object$K,
    method  = object$method,
    z_levels = object$z_levels,
    coefficients = list(),
    deviance     = if (!is.null(object$naive_fit)) object$naive_fit$deviance,
    df.residual  = if (!is.null(object$naive_fit)) object$naive_fit$df.residual,
    aic_naive    = if (!is.null(object$naive_fit)) stats::AIC(object$naive_fit),
    loglik_onestep = object$loglik_onestep
  )
  for (nm in names(object$coefficients)) {
    est <- object$coefficients[[nm]]
    V   <- if (!is.null(object$vcov)) object$vcov[[nm]] else NULL
    if (!is.null(V)) {
      se <- sqrt(pmax(diag(V), 0))
      z  <- est / se
      pv <- 2 * stats::pnorm(-abs(z))
      tab <- cbind(Estimate = est,
                   `Std. Error` = se,
                   `z value` = z,
                   `Pr(>|z|)` = pv)
    } else {
      tab <- cbind(Estimate = est)
    }
    ans$coefficients[[nm]] <- tab
  }
  class(ans) <- "summary.mcglm"
  ans
}


#' @export
print.summary.mcglm <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\nCall:\n")
  if (!is.null(x$call)) print(x$call) else cat("mcglm(...)\n")
  cat("\nFamily: ", x$family,
      "  |  n = ", x$n, ", K = ", x$K, ", p = ", x$p, "\n", sep = "")
  cat("Methods: ", paste(x$method, collapse = ", "), "\n", sep = "")
  if (!is.null(x$z_levels))
    cat("z categories (Pi assumed in this order): ",
        paste0(x$z_levels, c(" (baseline)", rep("", length(x$z_levels) - 1L)),
               collapse = ", "), "\n", sep = "")

  for (nm in x$method) {
    cat("\n--- ", toupper(nm), " ---\n", sep = "")
    tab <- x$coefficients[[nm]]
    if (is.null(tab)) next
    has_p <- ncol(tab) >= 4L
    stats::printCoefmat(tab, digits = digits,
                        has.Pvalue = has_p, P.values = has_p,
                        signif.stars = signif.stars && has_p,
                        cs.ind = if (has_p) 1:2 else 1L,
                        tst.ind = if (has_p) 3L else integer(0),
                        na.print = "NA")
  }

  if (!is.null(x$deviance)) {
    cat("\nResidual deviance (naive): ",
        format(signif(x$deviance, digits)),
        "  on ", x$df.residual, " degrees of freedom\n", sep = "")
    cat("AIC (naive): ", format(signif(x$aic_naive, digits)), "\n", sep = "")
  }
  if (!is.null(x$loglik_onestep)) {
    cat("logLik (onestep): ", format(signif(x$loglik_onestep, digits)), "\n",
        sep = "")
  }
  if (length(x$coefficients) > 1L) {
    cat("\nBias correction (difference from naive):\n")
    naive <- x$coefficients[["naive"]][, "Estimate"]
    others <- setdiff(names(x$coefficients), "naive")
    diffs <- sapply(others, function(o)
      x$coefficients[[o]][, "Estimate"] - naive)
    if (!is.matrix(diffs)) diffs <- matrix(diffs, nrow = 1L,
                                           dimnames = list(names(naive), others))
    print.default(format(diffs, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  invisible(x)
}


#' Standard errors of coefficient estimates
#'
#' @param object An \code{mcglm} object.
#' @param method Optional method name.
#' @return Named numeric vector of standard errors.
#' @export
se.mcglm <- function(object, method = NULL) {
  if (is.null(object$se)) return(sqrt(diag(vcov(object, method = method))))
  if (is.null(method)) method <- .mcglm_default_method(object)
  object$se[[method]]
}


#' Wald confidence intervals for an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param parm Optional names or indices of parameters to return.
#' @param level Coverage level.
#' @param method Optional estimation method (defaults to the last one fit,
#'   e.g. \code{"onestep"} if present).
#' @param ... Unused.
#' @export
confint.mcglm <- function(object, parm, level = 0.95, method = NULL, ...) {
  if (is.null(method)) method <- .mcglm_default_method(object)
  est <- coef(object, method = method)
  V   <- vcov(object, method = method)
  if (is.null(V))
    stop("No variance available for method '", method, "'.")
  se   <- sqrt(pmax(diag(V), 0))
  alpha <- (1 - level) / 2
  q    <- stats::qnorm(1 - alpha)
  ci   <- cbind(est - q * se, est + q * se)
  pct  <- paste0(format(c(alpha, 1 - alpha) * 100,
                        trim = TRUE, scientific = FALSE,
                        digits = 3), " %")
  colnames(ci) <- pct
  rownames(ci) <- names(est)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}


#' Fitted values for an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Estimation method.
#' @param ... Unused.
#' @export
fitted.mcglm <- function(object, method = NULL, ...) {
  if (isTRUE(object$is_multinomial))
    stop("fitted() not yet implemented for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(object)
  psi <- coef(object, method = method)
  fam <- .normalize_family(object$family)
  eta <- as.numeric(object$xi_hat %*% unname(psi))
  setNames(fam$linkinv(eta), seq_along(eta))
}


#' Linear-predictor or response predictions
#'
#' Currently supports in-sample prediction only.
#' @param object An \code{mcglm} object.
#' @param newdata Optional new data. For formula fits, a data frame
#'   containing the \code{mc()} variable (interpreted as the covariate at
#'   which to predict -- pass the \emph{true} category when available,
#'   e.g. when predicting onto a validation or probability sample) and
#'   the other covariates. For matrix-interface fits, a list with
#'   components \code{z_hat} (integer, coded \code{0, ..., K-1}) and
#'   \code{x} (covariate matrix as in the fit). \code{NULL} (default)
#'   predicts at the in-sample proxy design \eqn{\hat\xi_i}.
#' @param type One of \code{"link"}, \code{"response"}.
#' @param method Estimation method.
#' @param ... Unused.
#' @export
predict.mcglm <- function(object, newdata = NULL,
                          type = c("link", "response"),
                          method = NULL, ...) {
  type <- match.arg(type)
  if (isTRUE(object$is_multinomial))
    stop("predict() not yet implemented for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(object)
  psi <- unname(coef(object, method = method))

  if (is.null(newdata)) {
    xi <- object$xi_hat
  } else if (is.data.frame(newdata)) {
    if (is.null(object$formula))
      stop("predict() with a data-frame 'newdata' requires a fit from the ",
           "formula interface; for matrix-interface fits pass ",
           "newdata = list(z_hat = , x = ).")
    parts <- .mcglm_parse_formula(object$formula, newdata,
                                  environment(object$formula),
                                  require_y = FALSE,
                                  z_levels = object$z_levels,
                                  x_levels = object$x_levels)
    if (any(parts$z_hat < 0L) || any(parts$z_hat >= object$K))
      stop("newdata contains categories outside the fitted model ",
           "(K = ", object$K, ").")
    xi <- .mcglm_build_xi_hat(parts$z_hat, parts$x, object$K)
  } else if (is.list(newdata) &&
             all(c("z_hat", "x") %in% names(newdata))) {
    .mcglm_check_z_hat(newdata$z_hat)
    z_new <- as.integer(newdata$z_hat)
    x_new <- as.matrix(newdata$x)
    if (any(z_new < 0L) || any(z_new >= object$K))
      stop("newdata$z_hat must be coded 0, ..., K-1 (K = ", object$K, ").")
    xi <- .mcglm_build_xi_hat(z_new, x_new, object$K)
  } else {
    stop("'newdata' must be NULL, a data frame (formula fits), or a list ",
         "with components 'z_hat' and 'x' (matrix fits).")
  }

  if (ncol(xi) != object$p)
    stop("newdata produces ", ncol(xi), " columns; the model has ",
         object$p, " coefficients.")
  eta <- as.numeric(xi %*% psi)
  if (type == "link") return(eta)
  fam <- .normalize_family(object$family)
  fam$linkinv(eta)
}


#' Residuals for an \code{mcglm} fit
#'
#' @param object An \code{mcglm} object.
#' @param method Estimation method.
#' @param type One of \code{"response"}, \code{"pearson"}, \code{"deviance"}.
#' @param ... Unused.
#' @export
residuals.mcglm <- function(object, method = NULL,
                            type = c("response", "pearson", "deviance"),
                            ...) {
  type <- match.arg(type)
  if (isTRUE(object$is_multinomial))
    stop("residuals() not yet implemented for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(object)
  fam      <- .normalize_family(object$family)
  mu_val   <- fitted(object, method = method)
  r        <- object$y - mu_val
  if (type == "response") return(r)
  if (type == "pearson") return(r / sqrt(fam$variance(mu_val)))
  wt <- if (is.null(object$weights)) rep.int(1, length(r)) else object$weights
  d  <- fam$dev.resids(object$y, mu_val, wt = wt)
  sign(r) * sqrt(d)
}


#' Number of observations
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
nobs.mcglm <- function(object, ...) object$n


#' Family of an \code{mcglm} fit
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
family.mcglm <- function(object, ...) {
  .normalize_family(object$family)
}


#' Formula of an \code{mcglm} fit (NULL for matrix interface)
#' @param x An \code{mcglm} object.
#' @param ... Unused.
#' @export
formula.mcglm <- function(x, ...) x$formula


#' Proxy design matrix used by \code{mcglm}
#'
#' Returns the assembled \eqn{\hat\xi_i = (\hat d_i, x_i)} design matrix
#' (proxy-encoded misclassified covariate plus other covariates).
#' @param object An \code{mcglm} object.
#' @param ... Unused.
#' @export
model.matrix.mcglm <- function(object, ...) object$xi_hat


#' Log-likelihood of an \code{mcglm} fit
#'
#' Returns the naive-GLM log-likelihood for \code{method = "naive"} (the
#' only estimator with a closed-form likelihood for proxy data) or the
#' integrated mixture log-likelihood for \code{method = "onestep"}. Other
#' methods do not correspond to a likelihood; an error is thrown.
#' @param object An \code{mcglm} object.
#' @param method Estimation method.
#' @param ... Unused.
#' @export
logLik.mcglm <- function(object, method = NULL, ...) {
  if (is.null(method))
    method <- if ("onestep" %in% object$method) "onestep" else "naive"
  if (method == "onestep" && !is.null(object$loglik_onestep)) {
    val <- object$loglik_onestep
    df  <- if (!is.null(object$npar_onestep)) {
      object$npar_onestep
    } else {
      length(coef(object, method = "onestep"))
    }
    attr(val, "df")   <- df
    attr(val, "nobs") <- object$n
    class(val) <- "logLik"
    return(val)
  }
  if (method == "naive" && !is.null(object$naive_fit)) {
    return(stats::logLik(object$naive_fit))
  }
  stop("logLik not defined for method '", method, "'. ",
       "Available: ", paste(intersect(c("naive", "onestep"), object$method),
                            collapse = ", "))
}


#' AIC of an \code{mcglm} fit
#' @param object An \code{mcglm} object.
#' @param ... Additional arguments (ignored).
#' @param k Penalty per parameter (default 2).
#' @export
AIC.mcglm <- function(object, ..., k = 2) {
  ll <- logLik(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}


#' BIC of an \code{mcglm} fit
#' @param object An \code{mcglm} object.
#' @param ... Additional arguments (ignored).
#' @export
BIC.mcglm <- function(object, ...) {
  ll <- logLik(object)
  -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df")
}
