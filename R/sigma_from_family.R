#' Compute sigma as a function of eta for a GLM family
#'
#' For a GLM family with variance function V(mu), link function g, and
#' dispersion parameter phi, the sampling standard deviation of the MLE
#' on the link scale is:
#'   sigma(eta) = sqrt(phi * V(mu)) / |dmu/deta|
#' where mu = linkinv(eta + offset).
#'
#' @param eta Numeric vector of values on the link scale (parameter of interest).
#' @param family An R family object (e.g., binomial(), MASS::negative.binomial(0.5)).
#' @param offset Numeric scalar. Shift applied before evaluating the variance
#'   function. For NB/Poisson with log link, this is log(exposure). Default 0.
#' @param phi Numeric scalar. Dispersion parameter. For gaussian family this is
#'   sigma^2 (user-supplied). For binomial/NB/Poisson this is 1 (the default).
#' @return Numeric vector of sigma values (same length as eta).
#'
#' @keywords internal
sigma_from_family <- function(eta, family, offset = 0, phi = 1) {
  eta_full <- eta + offset
  mu <- family$linkinv(eta_full)
  mu_eta <- pmax(abs(family$mu.eta(eta_full)), .Machine$double.eps)
  sqrt(phi * family$variance(mu)) / mu_eta
}

#' Resolve family/sigma arguments into a sigma_fun closure
#'
#' Validates the family/sigma combination and returns a list with:
#' \itemize{
#'   \item \code{sigma_fun}: NULL (fixed-sigma path) or a closure
#'         \code{function(eta)} returning sigma at each eta.
#'   \item \code{family}: the family object (possibly set to NULL for
#'         gaussian(identity) short-circuit).
#' }
#'
#' @param family A family object or NULL.
#' @param sigma_missing Logical: was sigma missing in the calling function?
#' @param sigma Numeric sigma value (only used when sigma_missing is FALSE).
#' @param offset Numeric scalar offset.
#' @return A list with components \code{sigma_fun} and \code{family}.
#'
#' @keywords internal
resolve_sigma_family <- function(family, sigma_missing, sigma = NULL, offset = 0) {
  sigma_fun <- NULL
  if (!is.null(family)) {
    if (!inherits(family, "family")) {
      stop("'family' must be an R family object (e.g., binomial(), gaussian()).")
    }
    assert_number(offset)
    if (family$family == "gaussian") {
      if (sigma_missing) {
        stop("For gaussian() family, 'sigma' must be specified (dispersion parameter).")
      }
      if (family$link == "identity") {
        family <- NULL
      } else {
        phi <- sigma^2
        sigma_fun <- function(eta) sigma_from_family(eta, family, offset, phi)
      }
    } else {
      if (!sigma_missing) {
        stop("Cannot specify 'sigma' with a non-gaussian family. Sigma is determined by the family.")
      }
      sigma_fun <- function(eta) sigma_from_family(eta, family, offset)
    }
  }
  list(sigma_fun = sigma_fun, family = family)
}
