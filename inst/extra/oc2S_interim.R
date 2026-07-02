#' OC for Two-Stage Normal-Endpoint Designs with Interim Analysis
#'
#' Computes operating characteristics for a two-stage design with an
#' interim analysis (IA) and a user-defined continuation rule.
#'
#' The approach avoids brute-force Monte Carlo by:
#' 1. Computing the full-trial decision boundary once (the expensive step).
#' 2. Using Gauss-Hermite quadrature over the interim data distribution
#'    to integrate out the IA outcomes deterministically.
#' 3. At each quadrature point, evaluating the conditional OC via the
#'    pre-computed boundary — just pnorm + boundary lookup.
#'
#' @param prior1 normMix prior for arm 1 (treatment) used at final analysis
#' @param prior2 normMix prior for arm 2 (control) used at final analysis
#' @param n1 total sample size arm 1
#' @param n2 total sample size arm 2
#' @param n1_ia sample size arm 1 at interim
#' @param n2_ia sample size arm 2 at interim
#' @param decision decision2S object for the final analysis
#' @param ia_rule function receiving interim posteriors and returning
#'   TRUE to continue or FALSE to stop. Called as
#'   ia_rule(post1_ia, post2_ia) when no IA priors are given, or
#'   ia_rule(post1_ia, post2_ia, post1_ia_info, post2_ia_info)
#'   when prior1_ia / prior2_ia are provided. The first pair are
#'   posteriors from the analysis priors; the second from the IA priors.
#' @param sigma1 reference scale for arm 1 (default: from prior)
#' @param sigma2 reference scale for arm 2 (default: from prior)
#' @param family optional GLM family (e.g., negative.binomial(theta))
#' @param offset1 offset for arm 1 (log-exposure for NB/Poisson)
#' @param offset2 offset for arm 2
#' @param prior1_ia optional normMix prior for arm 1 at IA. When
#'   provided, an additional set of interim posteriors is computed
#'   from this prior and passed to ia_rule as the third argument.
#' @param prior2_ia optional normMix prior for arm 2 at IA. Same
#'   as prior1_ia but for arm 2 (fourth argument to ia_rule).
#' @param Ngrid_ia number of Gauss-Hermite nodes per arm at IA
#' @param eps tail-probability cutoff for boundary computation
#' @param Ngrid number of grid points for decision boundary
#'
#' @return A function f(theta1, theta2) returning a named numeric
#'   vector (scalar input) or data.frame (vector input) with
#'   components: power, stop_prob.
oc2S_interim <- function(
    prior1,
    prior2,
    n1,
    n2,
    n1_ia,
    n2_ia,
    decision,
    ia_rule,
    sigma1,
    sigma2,
    family = NULL,
    offset1 = 0,
    offset2 = offset1,
    prior1_ia = NULL,
    prior2_ia = NULL,
    Ngrid_ia = 21L,
    eps = 1e-6,
    Ngrid = 10
) {
  stopifnot(inherits(prior1, "normMix"))
  stopifnot(inherits(prior2, "normMix"))
  stopifnot(n1_ia < n1, n2_ia < n2)
  stopifnot(is.function(ia_rule))

  # --- Resolve sigma / family ---
  # Build sigma function from family if provided (public API only)
  sigma_fun1 <- NULL
  sigma_fun2 <- NULL
  if (!is.null(family)) {
    stopifnot(inherits(family, "family"))
    .sigma_from_family <- function(eta, fam, offset) {
      eta_full <- eta + offset
      mu <- fam$linkinv(eta_full)
      mu_eta <- pmax(abs(fam$mu.eta(eta_full)), .Machine$double.eps)
      sqrt(fam$variance(mu)) / mu_eta
    }
    sigma_fun1 <- function(eta) .sigma_from_family(eta, family, offset1)
    sigma_fun2 <- function(eta) .sigma_from_family(eta, family, offset2)
  }

  if (is.null(sigma_fun1)) {
    if (missing(sigma1)) sigma1 <- RBesT::sigma(prior1)
    if (missing(sigma2)) sigma2 <- RBesT::sigma(prior2)
  } else {
    sigma1 <- RBesT::sigma(prior1)
    sigma2 <- RBesT::sigma(prior2)
  }
  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  # --- IA-specific priors (for informative interim decisions) ---
  has_ia_priors <- !is.null(prior1_ia) || !is.null(prior2_ia)
  if (is.null(prior1_ia)) prior1_ia <- prior1
  if (is.null(prior2_ia)) prior2_ia <- prior2
  sigma(prior1_ia) <- sigma1
  sigma(prior2_ia) <- sigma2

  # --- Compute full-trial decision boundary ONCE ---
  if (is.null(sigma_fun1)) {
    crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision,
                                   sigma1, sigma2, eps, Ngrid)
  } else {
    crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision,
                                   family = family, offset1 = offset1,
                                   offset2 = offset2, eps = eps, Ngrid = Ngrid)
  }

  lower.tail <- attr(decision, "lower.tail")
  stopifnot(is(decision, "decision2S_1sided"))

  # --- Gauss-Hermite rule ---
  gh <- statmod::gauss.quad(Ngrid_ia, kind = "hermite")
  gh_nodes <- gh$nodes
  gh_weights <- gh$weights / sqrt(pi)  # normalised for N(0,1) measure

  # --- Remaining sample sizes ---
  n1_rem <- n1 - n1_ia
  n2_rem <- n2 - n2_ia

  # --- Returned OC function ---
  oc_fun <- function(theta1, theta2) {
    stopifnot(length(theta1) == length(theta2))

    results <- vapply(seq_along(theta1), function(idx) {
      t1 <- theta1[idx]
      t2 <- theta2[idx]

      # SE of interim MLE under true parameters
      se1_ia <- if (is.null(sigma_fun1)) {
        sigma1 / sqrt(n1_ia)
      } else {
        sigma_fun1(t1) / sqrt(n1_ia)
      }
      se2_ia <- if (is.null(sigma_fun2)) {
        sigma2 / sqrt(n2_ia)
      } else {
        sigma_fun2(t2) / sqrt(n2_ia)
      }

      # GH grid for interim outcomes centred on true values
      y1_grid <- t1 + sqrt(2) * se1_ia * gh_nodes
      y2_grid <- t2 + sqrt(2) * se2_ia * gh_nodes

      # Conditional SE of full-sample mean given interim
      # (only the remaining-data component contributes variance)
      se1_full_cond <- if (is.null(sigma_fun1)) {
        sigma1 * sqrt(n1_rem) / n1
      } else {
        sigma_fun1(t1) * sqrt(n1_rem) / n1
      }
      se2_full_cond <- if (is.null(sigma_fun2)) {
        sigma2 * sqrt(n2_rem) / n2
      } else {
        sigma_fun2(t2) * sqrt(n2_rem) / n2
      }

      # Pre-cache boundary over the range we'll need
      y1_full_range <- c(
        (n1_ia * min(y1_grid) + n1_rem * t1) / n1 - 4 * se1_full_cond,
        (n1_ia * max(y1_grid) + n1_rem * t1) / n1 + 4 * se1_full_cond
      )
      y2_full_range <- c(
        (n2_ia * min(y2_grid) + n2_rem * t2) / n2 - 4 * se2_full_cond,
        (n2_ia * max(y2_grid) + n2_rem * t2) / n2 + 4 * se2_full_cond
      )
      crit_y1(y2_full_range, lim1 = y1_full_range)

      power <- 0
      stop_prob <- 0

      for (j in seq_along(y2_grid)) {
        for (i in seq_along(y1_grid)) {
          wt <- gh_weights[i] * gh_weights[j]

          # Interim posteriors (analysis priors)
          post1_ia <- postmix(prior1, m = y1_grid[i], se = se1_ia)
          post2_ia <- postmix(prior2, m = y2_grid[j], se = se2_ia)

          # IA rule: TRUE = continue, FALSE = stop
          if (has_ia_priors) {
            post1_ia_info <- postmix(prior1_ia, m = y1_grid[i], se = se1_ia)
            post2_ia_info <- postmix(prior2_ia, m = y2_grid[j], se = se2_ia)
            ia_continue <- ia_rule(post1_ia, post2_ia, post1_ia_info, post2_ia_info)
          } else {
            ia_continue <- ia_rule(post1_ia, post2_ia)
          }
          if (!ia_continue) {
            stop_prob <- stop_prob + wt
            next
          }

          # Conditional mean of full-sample statistic
          mu1_cond <- (n1_ia * y1_grid[i] + n1_rem * t1) / n1
          mu2_cond <- (n2_ia * y2_grid[j] + n2_rem * t2) / n2

          # Evaluate conditional OC using pre-computed boundary:
          # integrate P(y1_full < B(y2_full)) * f(y2_full) dy2_full
          # Direct GH quadrature over the single Normal component
          x_inner <- mu2_cond + sqrt(2) * se2_full_cond * gh_nodes
          log_vals <- pnorm(crit_y1(x_inner, lim1 = y1_full_range),
                            mu1_cond, se1_full_cond,
                            lower.tail = lower.tail, log.p = TRUE)
          oc_val <- sum(gh_weights * exp(log_vals))
          power <- power + wt * oc_val
        }
      }

      c(power = power, stop_prob = stop_prob)
    }, numeric(2))

    # Return as data.frame for vector inputs
    if (length(theta1) == 1) {
      return(results[, 1])
    }
    data.frame(
      theta1 = theta1,
      theta2 = theta2,
      power = results["power", ],
      stop_prob = results["stop_prob", ]
    )
  }

  oc_fun
}
