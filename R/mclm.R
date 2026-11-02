#' Bias-Corrected Linear Model with Misclassified Covariates
#'
#' Convenience wrapper around \code{\link{mcglm}} that fixes the family
#' to \code{"gaussian"} (identity link), giving a misclassification-corrected
#' linear-regression interface analogous to \code{\link[stats]{lm}}. All
#' arguments other than \code{family} are forwarded to \code{mcglm} and
#' carry the same meaning.
#'
#' @inheritParams mcglm
#' @inherit mcglm return
#'
#' @seealso \code{\link{mcglm}} for the general GLM interface and the
#'   underlying methodology; \code{\link{simex}} for SIMEX/MC-SIMEX
#'   alternatives.
#'
#' @examples
#' set.seed(33)
#' n  <- 1000
#' Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
#' z  <- rbinom(n, 1, 0.4)
#' z_hat <- z
#' z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
#' z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.1)
#' x1 <- rnorm(n)
#' y  <- -0.5 + 0.8 * z + 0.7 * x1 + rnorm(n, 0, 0.5)
#' df <- data.frame(y = y, z = factor(z_hat), x1 = x1)
#'
#' fit <- mclm(y ~ mc(z, Pi) + x1, data = df)
#' summary(fit)
#'
#' @export
mclm <- function(formula, data = NULL,
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
                 optim_control = list(), validation = NULL,
                 z_hat = NULL, x = NULL) {

  cl <- match.call()
  mcall <- cl
  mcall[[1L]] <- quote(mismeasured::mcglm)
  mcall$family <- "gaussian"
  out <- eval(mcall, parent.frame())
  out$call <- cl
  out
}
