#' Gaussian quadrature for mixture density integration
#'
#' Provides S3 methods for \code{integrate_density} and
#' \code{integrate_density_log} that use Gaussian quadrature matched to
#' the mixture family:
#' \itemize{
#'   \item \code{normMix}: Gauss-Hermite quadrature
#'   \item \code{betaMix}: Gauss-Jacobi quadrature
#'   \item \code{gammaMix}: Gauss-Laguerre quadrature
#' }
#'
#' The integration method is controlled by the global option
#' \code{RBesT.integrate_method} (default \code{"GQ"}). Setting it to
#' \code{"adaptive"} falls through to the default adaptive method.
#'
#' Tolerance control is enabled by default: starting from
#' \code{RBesT.GQ_nodes} nodes, the node count is repeatedly increased
#' (by \code{RBesT.GQ_node_growth}, default doubling) until successive
#' estimates agree within \code{max(RBesT.GQ_abs_tol, RBesT.GQ_rel_tol * |I|)}.
#' Setting \code{RBesT.GQ_rel_tol} to a non-finite or non-positive value
#' (e.g. \code{Inf}) restores the legacy single-shot evaluation at
#' \code{RBesT.GQ_nodes} nodes.
#'
#' Relevant options (all consulted only on the \code{"GQ"} path):
#' \itemize{
#'   \item \code{RBesT.GQ_nodes} (default \code{20L}) -- starting node count
#'   \item \code{RBesT.GQ_rel_tol} (default \code{1e-4}) -- relative tolerance;
#'     a non-finite or non-positive value (e.g. \code{Inf}) disables refinement
#'     (single evaluation at \code{RBesT.GQ_nodes})
#'   \item \code{RBesT.GQ_abs_tol} (default \code{1e-6}) -- absolute tolerance floor
#'   \item \code{RBesT.GQ_max_nodes} (default \code{240L}) -- refinement cap
#'   \item \code{RBesT.GQ_node_growth} (default \code{2}) -- node growth factor
#'   \item \code{RBesT.GQ_on_nonconvergence} (default \code{"adaptive"}) -- one of
#'     \code{"adaptive"} (fall through to adaptive integration), \code{"warn"},
#'     \code{"error"}, \code{"silent"}
#' }
#'
#' @name integrate_gauss_quadrature
#' @keywords internal
NULL

# --- Tolerance-controlled refinement driver ---------------------------------

#' Refine a Gaussian-quadrature estimate to a requested tolerance
#'
#' Repeatedly evaluates \code{eval_at_n(n)} at a geometrically growing node
#' count until successive estimates agree within
#' \code{max(abs_tol, rel_tol * |I|)}. When \code{rel_tol} is \code{NULL} or
#' non-finite/non-positive (e.g. \code{Inf}) a single evaluation at \code{n0}
#' nodes is returned (legacy behaviour).
#'
#' @param eval_at_n function mapping an integer node count to a scalar estimate
#'   of the (summed-over-components) integral.
#' @param fallback_fn optional thunk returning the adaptive \code{.default}
#'   result for the same integral; used when \code{on_fail == "adaptive"}.
#' @details User-controllable arguments (sourced from \code{RBesT.GQ_*}
#'   options) are validated: \code{n0} must be a single integerish value
#'   \code{>= 10} and \code{< max_n}, \code{rel_tol} and \code{abs_tol} single
#'   non-negative doubles, \code{growth} a single double \code{> 1}, and
#'   \code{on_fail} one of \code{"adaptive"}, \code{"warn"}, \code{"error"},
#'   \code{"silent"}.
#' @keywords internal
#' @noRd
.gq_integrate_to_tol <- function(
  eval_at_n,
  fallback_fn = NULL,
  n0 = getOption("RBesT.GQ_nodes", 20L),
  rel_tol = getOption("RBesT.GQ_rel_tol", 1e-4),
  abs_tol = getOption("RBesT.GQ_abs_tol", 1e-6),
  max_n = getOption("RBesT.GQ_max_nodes", 240L),
  growth = getOption("RBesT.GQ_node_growth", 2),
  on_fail = getOption("RBesT.GQ_on_nonconvergence", "adaptive")
) {
  assert_int(n0, lower = 10L)
  n0 <- as.integer(n0)
  assert_number(rel_tol, lower = 0)
  assert_number(abs_tol, lower = 0)
  if (!is.finite(rel_tol)) {
    return(eval_at_n(n0))
  }
  assert_int(max_n, lower = n0 + 1L)
  assert_number(growth, lower = 1)
  stopifnot(growth > 1)
  assert_choice(on_fail, c("adaptive", "warn", "error", "silent"))
  max_n <- as.integer(max_n)
  n <- n0
  prev <- eval_at_n(n)
  repeat {
    n_next <- min(max_n, as.integer(ceiling(n * growth)))
    cur <- eval_at_n(n_next)
    err <- abs(cur - prev)
    if (err <= max(abs_tol, rel_tol * abs(cur))) {
      return(cur)
    }
    if (n_next >= max_n) {
      return(.gq_handle_nonconvergence(cur, err, on_fail, fallback_fn))
    }
    prev <- cur
    n <- n_next
  }
}

#' Apply the non-convergence policy for the GQ refinement driver
#' @keywords internal
#' @noRd
.gq_handle_nonconvergence <- function(best, err, policy, fallback_fn = NULL) {
  msg <- sprintf(
    "Gauss-quadrature refinement did not reach the requested tolerance (achieved error %.3g at the node cap).",
    err
  )
  switch(
    policy,
    silent = best,
    warn = {
      warning(msg, call. = FALSE)
      best
    },
    error = stop(msg, call. = FALSE),
    adaptive = {
      if (is.null(fallback_fn)) {
        warning(msg, call. = FALSE)
        best
      } else {
        fallback_fn()
      }
    },
    {
      warning(msg, call. = FALSE)
      best
    }
  )
}

# --- GH node/weight cache ---------------------------------------------------

.GH_cache <- new.env(parent = emptyenv())

.get_GH_rule <- function(n) {
  key <- as.character(n)
  if (!exists(key, envir = .GH_cache)) {
    rule <- statmod::gauss.quad(n, kind = "hermite")
    assign(key, rule, envir = .GH_cache)
  }
  get(key, envir = .GH_cache)
}

# --- normMix methods ---------------------------------------------------------

#' @exportS3Method
#' @describeIn integrate_density Gauss-Hermite method for normMix
#' @keywords internal
integrate_density.normMix <- function(mix, integrand, ...) {
  method <- getOption("RBesT.integrate_method", "GQ")
  if (method == "adaptive") {
    return(NextMethod())
  }

  eval_at_n <- function(n) {
    gh <- .get_GH_rule(n)
    t_nodes <- gh$nodes
    weights <- gh$weights
    Nc <- ncol(mix)
    total <- 0
    for (k in seq_len(Nc)) {
      x_k <- mix[2, k] + sqrt(2) * mix[3, k] * t_nodes
      g_vals <- integrand(x_k)
      total <- total + mix[1, k] / sqrt(pi) * sum(weights * g_vals)
    }
    total
  }
  fallback_fn <- function() integrate_density.default(mix, integrand, ...)
  .gq_integrate_to_tol(eval_at_n, fallback_fn = fallback_fn)
}

#' @exportS3Method
#' @describeIn integrate_density_log Gauss-Hermite method for normMix
#' @keywords internal
integrate_density_log.normMix <- function(mix, log_integrand, ...) {
  method <- getOption("RBesT.integrate_method", "GQ")
  if (method == "adaptive") {
    return(NextMethod())
  }

  integrate_density.normMix(mix, function(x) exp(log_integrand(x)))
}

# --- betaMix methods (Gauss-Jacobi) ------------------------------------------

#' @exportS3Method
#' @describeIn integrate_density_log Gauss-Jacobi method for betaMix
#' @keywords internal
integrate_density_log.betaMix <- function(mix, log_integrand, ...) {
  method <- getOption("RBesT.integrate_method", "GQ")
  if (method == "adaptive") {
    return(NextMethod())
  }
  if (dlink(mix)$name != "identity") {
    return(NextMethod())
  }

  eval_at_n <- function(n) {
    Nc <- ncol(mix)
    total <- 0
    for (k in seq_len(Nc)) {
      a_k <- mix[2, k]
      b_k <- mix[3, k]
      rule <- statmod::gauss.quad.prob(
        n,
        dist = "beta",
        alpha = a_k,
        beta = b_k
      )
      total <- total +
        mix[1, k] * sum(rule$weights * exp(log_integrand(rule$nodes)))
    }
    total
  }
  fallback_fn <- function() {
    integrate_density_log.default(mix, log_integrand, ...)
  }
  .gq_integrate_to_tol(eval_at_n, fallback_fn = fallback_fn)
}

# --- gammaMix methods (Gauss-Laguerre) ----------------------------------------

#' @exportS3Method
#' @describeIn integrate_density_log Gauss-Laguerre method for gammaMix
#' @keywords internal
integrate_density_log.gammaMix <- function(mix, log_integrand, ...) {
  method <- getOption("RBesT.integrate_method", "GQ")
  if (method == "adaptive") {
    return(NextMethod())
  }
  if (dlink(mix)$name != "identity") {
    return(NextMethod())
  }

  eval_at_n <- function(n) {
    Nc <- ncol(mix)
    total <- 0
    for (k in seq_len(Nc)) {
      shape_k <- mix[2, k]
      rate_k <- mix[3, k]
      rule <- statmod::gauss.quad.prob(
        n,
        dist = "gamma",
        alpha = shape_k,
        beta = 1 / rate_k
      )
      total <- total +
        mix[1, k] * sum(rule$weights * exp(log_integrand(rule$nodes)))
    }
    total
  }
  fallback_fn <- function() {
    integrate_density_log.default(mix, log_integrand, ...)
  }
  .gq_integrate_to_tol(eval_at_n, fallback_fn = fallback_fn)
}
