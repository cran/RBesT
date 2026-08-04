#' Decision Boundary for 2 Sample Designs
#'
#' The `decision2S_boundary` function defines a 2 sample design
#' (priors, sample sizes, decision function) for the calculation of
#' the decision boundary. A function is returned which calculates the
#' critical value of the first sample \eqn{y_{1,c}} as a function of
#' the outcome in the second sample \eqn{y_2}. At the decision
#' boundary, the decision function will change between 0 (failure) and
#' 1 (success) for the respective outcomes.
#'
#' @template args-boundary2S
#'
#' @details For a 2 sample design the specification of the priors, the
#' sample sizes and the decision function, \eqn{D(y_1,y_2)}, uniquely
#' defines the decision boundary
#'
#' \deqn{D_1(y_2) = \max_{y_1}\{D(y_1,y_2) = 1\},}{D_1(y_2) = max_{y_1}{D(y_1,y_2) = 1},}
#'
#' which is the critical value of \eqn{y_{1,c}} conditional on the
#' value of \eqn{y_2} whenever the decision \eqn{D(y_1,y_2)} function
#' changes its value from 0 to 1 for a decision function with
#' `lower.tail=TRUE` (otherwise the definition is \eqn{D_1(y_2) =
#' \max_{y_1}\{D(y_1,y_2) = 0\}}{D_1(y_2) = max_{y_1}{D(y_1,y_2) =
#' 0}}). The decision function may change at most at a single critical
#' value for given \eqn{y_{2}} as only one-sided decision functions
#' are supported. Here, \eqn{y_2} is defined for binary and Poisson
#' endpoints as the sufficient statistic \eqn{y_2 = \sum_{i=1}^{n_2}
#' y_{2,i}} and for the normal case as the mean \eqn{\bar{y}_2 = 1/n_2
#' \sum_{i=1}^{n_2} y_{2,i}}.
#'
#' @return For one-sided decision functions, returns a function with a
#' single argument. This function calculates in dependence of the
#' outcome \eqn{y_2} in sample 2 the critical value \eqn{y_{1,c}} for
#' which the defined design will change the decision from 0 to 1 (or
#' vice versa, depending on the decision function).
#' For two-sided decision functions, returns a list with components
#' `lower_or_equal_than` and `higher_than`, containing the critical
#' value functions for the lower and upper one-sided decision
#' boundary components.
#'
#' @family design2S
#'
#' @examples
#'
#' # see ?decision2S for details of example
#' priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
#' priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
#' # the success criteria is for delta which are larger than some
#' # threshold value which is why we set lower.tail=FALSE
#' successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
#' # the futility criterion acts in the opposite direction
#' futilityCrit <- decision2S(c(0.90), c(40), TRUE)
#'
#' # success criterion boundary
#' successBoundary <- decision2S_boundary(priorP, priorT, 10, 20, successCrit)
#'
#' # futility criterion boundary
#' futilityBoundary <- decision2S_boundary(priorP, priorT, 10, 20, futilityCrit)
#'
#' curve(successBoundary(x), -25:25 - 49, xlab = "y2", ylab = "critical y1")
#' curve(futilityBoundary(x), lty = 2, add = TRUE)
#'
#' # hence, for mean in sample 2 of 10, the critical value for y1 is
#' y1c <- futilityBoundary(-10)
#'
#' # around the critical value the decision for futility changes
#' futilityCrit(postmix(priorP, m = y1c + 1E-3, n = 10), postmix(priorT, m = -10, n = 20))
#' futilityCrit(postmix(priorP, m = y1c - 1E-3, n = 10), postmix(priorT, m = -10, n = 20))
#'
#' @export
decision2S_boundary <- function(prior1, prior2, n1, n2, decision, ...) {
  UseMethod("decision2S_boundary")
}
#' @export
decision2S_boundary.default <- function(prior1, prior2, n1, n2, decision, ...) {
  "Unknown density"
}

#' @templateVar fun decision2S_boundary
#' @template design2S-binomial
# If the \code{eps} argument is specificed, then the
# returned function will use the additional \code{lim2}
# argument to limit the search for the critical value.
#' @export
decision2S_boundary.betaMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  ## only n2=0 is supported
  assert_number(n1, lower = 1, finite = TRUE)
  assert_number(n2, lower = 0, finite = TRUE)

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_betaMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  } else {
    decision2S_boundary_betaMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  }
}

decision2S_boundary_betaMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  crit_lower <- decision2S_boundary_betaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    eps,
    ...
  )
  crit_upper <- decision2S_boundary_betaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    eps,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}

decision2S_boundary_betaMix_1sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  decision <- if (has_lower(decision)) {
    lower(decision)
  } else {
    upper(decision)
  }
  decision2S_boundary_betaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    eps,
    ...
  )
}

decision2S_boundary_betaMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  if (!missing(eps)) {
    assert_number(eps, lower = 0, upper = 0.1, finite = TRUE)
  }

  cond_decisionDist <- function(post2cond) {
    function(m1) {
      ## Note: Subtracting from the decision 0.25 leads to
      ## negative decisions being at -0.25 while positives are
      ## at 0.75; since uniroot_int *always* returns the x which
      ## has lowest absolute value we are guaranteed that y2crit
      ## is just before the jump
      ## decision(post1cond, post2[[m2+1]]) - 0.25
      decision(postmix(prior1, r = m1, n = n1), post2cond) - 0.25
      ## decision(postmix(prior2, r=m2, n=n2), post1cond) - 0.25
    }
  }

  ## saves the decision boundary conditional on the outcome of the
  ## second variable
  clim1 <- c(Inf, -Inf)
  clim2 <- c(Inf, -Inf)
  boundary <- c()
  full_boundary <- missing(eps)

  lower.tail <- attr(decision, "lower.tail")

  update_boundary <- function(lim1, lim2) {
    boundary <<- rep(NA, diff(lim2) + 1)
    clim2 <<- lim2
    clim1 <<- lim1
    ## The step objective is monotone in m1 in the direction set by
    ## lower.tail; the critical y1 itself is monotone in y2 for the
    ## supported one-sided decisions.
    extend_dir <- if (lower.tail) "downX" else "upX"

    decFun_for <- function(y2) {
      if (n2 == 0) {
        cond_decisionDist(prior2)
      } else {
        cond_decisionDist(postmix(prior2, r = y2, n = n2))
      }
    }
    ## raw integer critical value at y2 within a warm bracket;
    ## uniroot_int extends it (bounded by lim1) if the boundary moved
    ## further and reports +/-Inf where the decision is constant.
    solve_raw <- function(y2, bracket) {
      uniroot_int(decFun_for(y2), bracket, extendInt = extend_dir, clamp = lim1)
    }
    ## map the constant-decision sentinels into the ordered integer
    ## domain [-1, n1 + 1] so they participate in the monotone fill:
    ## never true (-Inf -> -1), always true (Inf -> n1 + 1).
    map_root <- function(root) {
      if (is.finite(root)) {
        root
      } else if (root < 0) {
        -1
      } else {
        n1 + 1
      }
    }
    set_val <- function(y2, val) {
      boundary[y2 - lim2[1] + 1] <<- val
    }

    if (getOption("RBesT.decision2S_boundary", "adaptive") == "adaptive") {
      ## Adaptive solve replacement. Because the critical value is
      ## monotone in y2, if the roots at the two ends of a y2 subinterval
      ## are equal then every y2 between them shares that value -- the run
      ## is filled without any (expensive pmixdiff) uniroot_int calls.
      ## Otherwise the midpoint is solved, tightly bracketed by the two
      ## neighbour roots (crit(m) is guaranteed to lie between them), and
      ## the two halves recurse. Steep boundaries end up solving every
      ## point; flat runs cost O(log length).
      bracket_for <- function(ra, rb) {
        lo <- max(min(ra, rb), lim1[1])
        hi <- min(max(ra, rb), lim1[2])
        if (lo >= hi) lim1 else c(lo, hi)
      }
      fill <- function(a, ra, b, rb) {
        if (ra == rb) {
          for (y2 in a:b) {
            set_val(y2, ra)
          }
          return(invisible())
        }
        if (b - a <= 1L) {
          set_val(a, ra)
          set_val(b, rb)
          return(invisible())
        }
        m <- (a + b) %/% 2L
        rm <- map_root(solve_raw(m, bracket_for(ra, rb)))
        fill(a, ra, m, rm)
        fill(m, rm, b, rb)
      }
      ra <- map_root(solve_raw(lim2[1], lim1))
      rb <- map_root(solve_raw(lim2[2], lim1))
      fill(lim2[1], ra, lim2[2], rb)
    } else {
      ## Legacy sequential sweep: warm-start each y2 point's integer
      ## search with a tight bracket around the previous boundary.
      prev <- NA
      for (y2 in lim2[1]:lim2[2]) {
        bracket <- if (is.na(prev)) {
          lim1
        } else {
          c(max(prev - 2, lim1[1]), min(prev + 2, lim1[2]))
        }
        root <- solve_raw(y2, bracket)
        set_val(y2, map_root(root))
        prev <- if (is.finite(root)) root else NA
      }
    }

    if (lower.tail) {
      ## if lower.tail==TRUE, then the condition becomes true when
      ## going from large to small values, hence we need to integrate from
      ## 0 to boundary
      boundary <<- pmax(boundary - 1, -1)
    }
    return()
  }

  if (full_boundary) {
    update_boundary(c(0, n1), c(0, n2))
  }

  decision_boundary <- function(y2, lim1) {
    ## check if we need to recalculate the decision grid for the
    ## case of enabled approximate method
    assert_integerish(y2, lower = 0, upper = n2, any.missing = FALSE)

    if (!full_boundary) {
      if (missing(lim1)) {
        ## if not hint is given we search the full sample
        ## space which should be OK, as the complexity is
        ## log(N)
        lim1 <- c(0, n1)
      } else {
        assert_integerish(lim1, lower = 0, upper = n1, any.missing = FALSE)
      }
      lim2 <- c(min(y2), max(y2))
      ## check if the decision grid needs to be recomputed
      if (
        lim1[1] < clim1[1] |
          lim1[2] > clim1[2] |
          lim2[1] < clim2[1] |
          lim2[2] > clim2[2]
      ) {
        ## ensure that lim1 never shrinks
        lim1[1] <- min(lim1[1], clim1[1])
        lim1[2] <- max(lim1[2], clim1[2])
        update_boundary(lim1, lim2)
      }
    }

    ## make sure y2 is an integer which is the value of
    ## the second read-out for which we return the decision
    ## boundary
    ## TODO: handle case with eps with care
    crit <- boundary[(y2 - clim2[1]) + 1]
    if (!full_boundary) {
      ## in case the lower boundary of the searched grid is not
      ## zero, then we cannot say anything about cases when the
      ## decision is always negative
      if (!lower.tail) {
        ## in this case the decision changes from negative to
        ## positive when going from small to large
        ## values. Hence, if the decision is always negative,
        ## then we can be sure of that we can never be sure,
        ## but should the decision be negative at all values,
        ## it can change at larger values.
        crit[crit == n1 + 1] <- NA
      } else {
        ## now the decision changes from positive to negative
        ## when going from small to large => should the
        ## decision not change in the clim1 domain then we do
        ## not know if it happens later
        if (clim1[1] > 0) {
          crit[crit == -1] <- NA
        }
        ## however, if crit==Inf then we can be sure that the
        ## decision is indeed always positive
      }
    }
    return(crit)
  }
  decision_boundary
}


## Adaptive boundary tracer for the normal endpoint.
##
## Recursively bisects the y2 interval instead of sweeping a fixed grid.
## Each new point is seeded from its two bracketing solved neighbours and
## root-found within their two roots as bracket -- guaranteed to contain
## the root because y1_c(y2) is monotone in y2 for the supported one-sided
## decisions; the solver's extendInt remains the safety net. delta2 is the
## resolution floor.
##
## Acceptance uses a nested (quarter-point) error estimate rather than a
## single-midpoint linearity check. For a subinterval [a, b] with solved
## endpoints and midpoint (a, ra), (m, rm), (b, rb) we form the quadratic
## Lagrange interpolant Q through those three points, solve the two quarter
## points q1 = (a+m)/2, q3 = (m+b)/2, and accept the subinterval when
## max(|rq1 - Q(q1)|, |rq3 - Q(q3)|) < lin_tol; otherwise we subdivide.
## Because acceptance is decided at points Q never saw, a boundary that is
## straight at the midpoint but bent between samples cannot alias as
## linear -- this removes the need for an explicit minimum-subdivision
## floor. The quarter points q1, q3 of a parent are exactly the midpoints
## of its two children, so with memoization no y2 is ever solved twice.
##
## The exact solved roots are fed to the downstream cubic splinefun. This
## is interpolation-consistent: recursive bisection produces a smoothly
## graded (dyadic) point density -- dense where the boundary bends, coarse
## where it is linear -- and the quadratic-at-quarter-points accept test is
## a faithful proxy for the cubic consumer, so lin_tol caps the residual
## curvature over any accepted subinterval. The cost saving is fewer root
## solves; accuracy matches or beats the uniform grid at the same delta2.
##
## Returns cbind(grid, crit), sorted and de-duplicated for splinefun, or
## NULL if a constant/no-root (NA) region is hit -- the caller then falls
## back to the robust uniform grid sweep.
trace_boundary2S_adaptive <- function(
  solve_at,
  lim1,
  lim2,
  delta2,
  scale1,
  lin_tol
) {
  ## memoized solver: one uniroot per distinct y2, ever -- a quarter point
  ## of a parent interval is reused as the midpoint of a child interval.
  cache_y <- numeric(0)
  cache_r <- numeric(0)
  bad <- FALSE

  solve_memo <- function(y2, bracket) {
    hit <- which(abs(cache_y - y2) < 1e-9)
    if (length(hit)) {
      return(cache_r[hit[1L]])
    }
    r <- solve_at(y2, bracket)
    if (is.na(r)) {
      bad <<- TRUE
    }
    cache_y[[length(cache_y) + 1L]] <<- y2
    cache_r[[length(cache_r) + 1L]] <<- r
    r
  }

  ## bracket for a point known to lie between two solved neighbours;
  ## monotonicity of y1_c(y2) guarantees the root is inside [r_lo, r_hi].
  ## A degenerate (equal-root) bracket falls back to a symmetric scale1
  ## window so uniroot's extendInt can recover.
  bracket_between <- function(r_lo, r_hi) {
    br <- sort(c(r_lo, r_hi))
    if (diff(br) < 1e-9) {
      mean(br) + c(-2, 2) * scale1
    } else {
      br
    }
  }

  ## quadratic Lagrange interpolant through (a, ra), (m, rm), (b, rb)
  quad <- function(x, a, ra, m, rm, b, rb) {
    ra *
      ((x - m) * (x - b)) /
      ((a - m) * (a - b)) +
      rm * ((x - a) * (x - b)) / ((m - a) * (m - b)) +
      rb * ((x - a) * (x - m)) / ((b - a) * (b - m))
  }

  ## recurse over [a, b]; the midpoint m is already solved and passed in
  ## (it was a quarter point of the parent -- reused, not re-solved).
  recurse <- function(a, ra, m, rm, b, rb) {
    if (bad || (b - a) <= delta2) {
      return(invisible())
    }
    q1 <- (a + m) / 2
    q3 <- (m + b) / 2
    rq1 <- solve_memo(q1, bracket_between(ra, rm))
    rq3 <- solve_memo(q3, bracket_between(rm, rb))
    if (bad) {
      return(invisible())
    }
    ## nested error estimate: solved roots vs the quadratic they did not build
    err <- max(
      abs(rq1 - quad(q1, a, ra, m, rm, b, rb)),
      abs(rq3 - quad(q3, a, ra, m, rm, b, rb))
    )
    if (err > lin_tol) {
      ## q1, q3 become the midpoints of the two children.
      recurse(a, ra, q1, rq1, m, rm)
      recurse(m, rm, q3, rq3, b, rb)
    }
  }

  a <- lim2[1]
  b <- lim2[2]
  ra <- solve_memo(a, lim1)
  rb <- solve_memo(b, lim1)
  if (!bad) {
    m <- (a + b) / 2
    rm <- solve_memo(m, bracket_between(ra, rb))
    if (!bad) {
      recurse(a, ra, m, rm, b, rb)
    }
  }

  if (bad) {
    return(NULL)
  }
  ord <- order(cache_y)
  y <- cache_y[ord]
  r <- cache_r[ord]
  keep <- c(TRUE, diff(y) > 1e-9)
  cbind(y[keep], r[keep])
}

## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object
solve_boundary2S_normMix <- function(
  decision,
  mix1,
  mix2,
  n1,
  n2,
  lim1,
  lim2,
  delta2,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  y2_at = NULL
) {
  assert_class(decision, "decision2S_atomic")

  sigma2 <- sigma(mix2)

  if (is.null(sigma_fun1)) {
    ## sigma is known and fixed
    sigma1 <- sigma(mix1)
    sem1 <- sigma1 / sqrt(n1)
    scale1 <- sigma1 / (n1^0.25)

    cond_decisionStep <-
      function(post2) {
        function(m1) {
          min(decision(postmix(mix1, m = m1, se = sem1), post2, dist = TRUE))
        }
      }
  } else {
    sigma1 <- sigma_fun1(summary(mix1)["mean"])
    sem1 <- sigma1 / sqrt(n1)
    scale1 <- sigma1 / (n1^0.25)

    ## Use the reference data sigma (sigma_fun evaluated at the prior
    ## mean) to set the search domain.  The boundary can only exist
    ## where the data is informative, i.e. where sigma_fun1(m1) is of
    ## similar magnitude to sigma1.  A factor of 8 is generous enough
    ## to never clip a genuine boundary.
    max_half_width <- max(8 * sigma1, diff(lim1))

    cond_decisionStep <- function(post2) {
      ## Family path: bound uniroot's interval extension to the region
      ## where the data is informative. Returning NA outside that window
      ## (centred on the current bracket) stops the extension cheaply; a
      ## genuinely constant decision then yields no root and the grid
      ## point is dropped. center is fixed per grid point because
      ## cond_decisionStep() is called with the current lim1.
      center <- mean(lim1)
      function(m1) {
        if (abs(m1 - center) > max_half_width) {
          return(NA_real_)
        }
        se_m1 <- sigma_fun1(m1) / sqrt(n1)
        min(decision(postmix(mix1, m = m1, se = se_m1), post2, dist = TRUE))
      }
    }
  }

  tol <- min(delta2 / 100, .Machine$double.eps^0.25)

  ## Direction of the (continuous) decision objective as a function of
  ## m1. min(decision(..., dist = TRUE)) is monotone in m1 and shares
  ## the atomic's lower.tail: for lower.tail = TRUE the objective
  ## decreases through the root ("downX"), for lower.tail = FALSE it
  ## increases ("upX"). We provide good starting limits and let uniroot
  ## extend the bracket itself instead of expanding lim1 by hand.
  extend_dir <- if (isTRUE(attr(decision, "lower.tail"))) "downX" else "upX"

  ## Root-find the critical y1 for a second-sample outcome y2 within a
  ## warm bracket (extendInt widens it if the boundary moved out).
  ## Returns NA where the decision is constant (no sign change).
  solve_at <- function(y2, bracket) {
    if (n2 == 0) {
      post2 <- mix2
    } else {
      se2 <- if (is.null(sigma_fun2)) {
        sigma2 / sqrt(n2)
      } else {
        sigma_fun2(y2) / sqrt(n2)
      }
      post2 <- postmix(mix2, m = y2, se = se2)
    }
    ind_fun <- cond_decisionStep(post2)
    tryCatch(
      uniroot(ind_fun, bracket, extendInt = extend_dir, tol = tol)$root,
      error = function(e) NA_real_
    )
  }

  ## Direct solve at caller-specified y2 points (used by the exact
  ## linear-boundary fast path, which needs only two points to define a
  ## provably-straight boundary). Bypasses both the adaptive tracer and the
  ## uniform grid. Points with no root (NA) are dropped by the caller.
  if (!is.null(y2_at)) {
    crit <- vapply(y2_at, solve_at, numeric(1), bracket = lim1)
    keep <- !is.na(crit)
    return(cbind(y2_at[keep], crit[keep]))
  }

  ## Adaptive tracing is the default; the uniform grid sweep is kept as a
  ## fallback (constant/no-root regions) and as an opt-out via
  ## options(RBesT.decision2S_boundary = "grid"). lin_tol is the accepted
  ## deviation of the boundary from its quadratic fit at the quarter points
  ## of a subinterval; it is set relative to sem1, the posterior scale of
  ## the arm-1 mean that governs the boundary's own resolution. The nested
  ## quarter-point accept test is aliasing-resistant, so sem1/100
  ## ("boundary located to ~1% of its natural sem1 scale") caps the
  ## worst-case error at grid level without the over-tight sem1/1000 the
  ## curvature-blind single-midpoint test required.
  if (getOption("RBesT.decision2S_boundary", "adaptive") == "adaptive") {
    adaptive <- trace_boundary2S_adaptive(
      solve_at,
      lim1,
      lim2,
      delta2,
      scale1,
      lin_tol = sem1 / 100
    )
    if (!is.null(adaptive)) {
      return(adaptive)
    }
  }

  ## Uniform grid sweep: a fixed grid of y2 points, warm-started with a
  ## symmetric +/-2*scale1 bracket around the previous root.
  grid <- seq(lim2[1], lim2[2], length = diff(lim2) / delta2)
  crit <- rep(NA_real_, length(grid))
  bracket <- lim1
  for (i in seq_along(grid)) {
    y1c <- solve_at(grid[i], bracket)
    if (is.na(y1c)) {
      bracket <- lim1
      next
    }
    crit[i] <- y1c
    bracket <- c(y1c - 2 * scale1, y1c + 2 * scale1)
  }

  ## Drop grid points where no boundary was found (crit is NA).
  ## The spline / approxfun covers only the range with real boundaries;
  ## extrapolation at the tails produces extreme values that give the
  ## correct pnorm result (~ 0 or 1).
  keep <- !is.na(crit)
  cbind(grid[keep], crit[keep])
}

#' @templateVar fun decision2S_boundary
#' @template design2S-normal
#' @export
decision2S_boundary.normMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps = 1e-6,
  Ngrid = 5,
  family = NULL,
  offset1 = 0,
  offset2 = offset1,
  ...
) {
  # Determine data sigma for each arm
  resolved1 <- resolve_sigma_family(family, missing(sigma1), sigma1, offset1)
  resolved2 <- resolve_sigma_family(family, missing(sigma2), sigma2, offset2)
  sigma_fun1 <- resolved1$sigma_fun
  sigma_fun2 <- resolved2$sigma_fun
  family <- resolved1$family

  # Resolve sigma for the fixed-sigma path
  if (is.null(sigma_fun1)) {
    if (missing(sigma1)) {
      sigma1 <- RBesT::sigma(prior1)
      message("Using default prior 1 reference scale ", sigma1)
    }
    if (missing(sigma2)) {
      sigma2 <- RBesT::sigma(prior2)
      message("Using default prior 2 reference scale ", sigma2)
    }
  } else {
    sigma1 <- RBesT::sigma(prior1)
    sigma2 <- RBesT::sigma(prior2)
  }

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_normMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      sigma_fun1 = sigma_fun1,
      sigma_fun2 = sigma_fun2,
      ...
    )
  } else {
    decision2S_boundary_normMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      sigma_fun1 = sigma_fun1,
      sigma_fun2 = sigma_fun2,
      ...
    )
  }
}

decision2S_boundary_normMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  ...
) {
  crit_lower <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    sigma_fun1 = sigma_fun1,
    sigma_fun2 = sigma_fun2,
    ...
  )
  crit_upper <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    sigma_fun1 = sigma_fun1,
    sigma_fun2 = sigma_fun2,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}

decision2S_boundary_normMix_1sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  ...
) {
  decision <- if (has_lower(decision)) {
    lower(decision)
  } else {
    upper(decision)
  }
  decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    sigma1,
    sigma2,
    eps,
    Ngrid,
    sigma_fun1 = sigma_fun1,
    sigma_fun2 = sigma_fun2,
    ...
  )
}

## Build the interpolating boundary function from the solved (y2, y1c)
## points produced by solve_boundary2S_normMix. This is the single seam
## through which the boundary consumer is chosen:
##   - empty input  -> constant extreme-value function (decision never/always
##     satisfied); the sign follows the decision's lower.tail so pnorm ~ 0/1.
##   - `linear`     -> piecewise-linear approxfun for the provably-linear
##     single-component / constant-sigma case.
##   - otherwise    -> smooth spline through the (non-uniformly graded) points.
fit_boundary_consumer <- function(boundary_discrete, linear, decision) {
  if (nrow(boundary_discrete) == 0) {
    lower.tail <- attr(decision, "lower.tail")
    const_val <- if (lower.tail) -1e100 else 1e100
    return(function(y2) rep(const_val, length(y2)))
  }
  if (linear) {
    approxfun(boundary_discrete[, 1], boundary_discrete[, 2], rule = 2)
  } else {
    ## Shape-preserving monotone Hermite (Fritsch-Carlson): matches the
    ## adaptive tracer's curvature-controlled point placement, avoids the
    ## overshoot a global cubic ("fmm") exhibits over the non-uniformly
    ## graded points, and respects the boundary's monotonicity where it
    ## holds without imposing it where it does not.
    splinefun(
      boundary_discrete[, 1],
      boundary_discrete[, 2],
      method = "monoH.FC"
    )
  }
}

decision2S_boundary_normMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  ...
) {
  assert_class(decision, "decision2S_atomic")

  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)

  sem1 <- sigma1 / sqrt(n1)
  sem2 <- sigma2 / sqrt(n2)

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  ## only n2 can be zero
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  if (n2 == 0) {
    sem2 <- sigma(prior2) / sqrt(1E-1)
  }

  ## change the reference scale of the prior such that the prior
  ## represents the distribution of the respective means
  mean_prior1 <- prior1
  sigma(mean_prior1) <- sem1

  ## discretization step-size
  delta2 <- sem2 / Ngrid

  ## for the case of mix1 and mix2 having just 1 component, then one
  ## can prove that the decision boundary is a linear function.
  ## Hence we only calculate a very rough grid and apply linear
  ## interpolation.
  ## With family path, linearity no longer holds even for single-component priors.

  linear_boundary <- FALSE
  if (ncol(prior1) == 1 && ncol(prior2) == 1 && is.null(sigma_fun1)) {
    linear_boundary <- TRUE
    ## we could relax this even further
    delta2 <- sigma2 / Ngrid
  }

  ## the boundary function depends only on the samples sizes n1, n2,
  ## the priors and the decision, but not the assumed truths

  clim2 <- c(Inf, -Inf)

  ## the boundary function which gives conditional on the second
  ## variable the critical value where the decision changes
  boundary <- NA
  boundary_discrete <- matrix(NA, nrow = 0, ncol = 2)
  exact_line_built <- FALSE
  linear_fallback <- FALSE

  decision_boundary <- function(y2, lim1) {
    ## Exact linear fast path (C1/C2/C3): for single-component priors with
    ## constant sigma the boundary y1_c(y2) is provably a straight line, and
    ## for multi-criteria decisions the per-criterion lines are parallel, so
    ## their min/max envelope is still a single line (no kink). We therefore
    ## solve at just two anchor points, fit the exact line, and cache it as a
    ## closed-form function valid for all y2 -- bypassing the grid sweep and
    ## the lazy caching/extension machinery entirely. n2 == 0 collapses to a
    ## horizontal line (two equal anchors -> slope 0) automatically.
    if (exact_line_built) {
      return(boundary(y2))
    }
    if (linear_boundary && !linear_fallback) {
      if (missing(lim1)) {
        lim1 <- qmix(mean_prior1, c(eps / 2, 1 - eps / 2))
      }
      anchor2 <- range(y2)
      if (diff(anchor2) == 0) {
        anchor2 <- anchor2 + c(-1, 1) * max(delta2, abs(anchor2[1]), 1)
      }
      pts <- solve_boundary2S_normMix(
        decision,
        prior1,
        prior2,
        n1,
        n2,
        lim1,
        anchor2,
        delta2,
        sigma_fun1 = sigma_fun1,
        sigma_fun2 = sigma_fun2,
        y2_at = anchor2
      )
      if (nrow(pts) == 2 && pts[2, 1] != pts[1, 1]) {
        boundary_discrete <<- pts
        slope <- (pts[2, 2] - pts[1, 2]) / (pts[2, 1] - pts[1, 1])
        intercept <- pts[1, 2] - slope * pts[1, 1]
        boundary <<- local({
          m <- slope
          b <- intercept
          function(y2) b + m * y2
        })
        clim2 <<- c(-Inf, Inf)
        exact_line_built <<- TRUE
        return(boundary(y2))
      }
      ## Two-anchor solve failed (e.g. a root did not exist at an anchor);
      ## fall through to the general grid/adaptive path, which stays linear
      ## via the approxfun consumer. Latch so later calls skip the retry.
      linear_fallback <<- TRUE
    }

    lim2 <- range(y2)

    ## check if boundary function must be recomputed
    if (lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
      new_lim2 <- clim2
      ## note: the <<- assignment is needed to set the variable in the enclosure
      if (missing(lim1)) {
        lim1 <- qmix(mean_prior1, c(eps / 2, 1 - eps / 2))
      }
      if (nrow(boundary_discrete) == 0) {
        ## boundary hasn't been calculated before, do it all
        boundary_discrete <<- solve_boundary2S_normMix(
          decision,
          prior1,
          prior2,
          n1,
          n2,
          lim1,
          lim2,
          delta2,
          sigma_fun1 = sigma_fun1,
          sigma_fun2 = sigma_fun2
        )
        new_lim2 <- lim2
      } else {
        if (lim2[1] < clim2[1]) {
          ## the lower bound is not low enough... only add the region which is missing
          new_left_lim2 <- min(lim2[1], clim2[1] - 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(
            decision,
            prior1,
            prior2,
            n1,
            n2,
            lim1,
            c(new_left_lim2, clim2[1] - delta2),
            delta2,
            sigma_fun1 = sigma_fun1,
            sigma_fun2 = sigma_fun2
          )
          new_lim2[1] <- new_left_lim2
          boundary_discrete <<- rbind(boundary_extra, boundary_discrete)
        }
        if (lim2[2] > clim2[2]) {
          ## the upper bound is not large enough.. again only add what's missing
          new_right_lim2 <- max(lim2[2], clim2[2] + 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(
            decision,
            prior1,
            prior2,
            n1,
            n2,
            lim1,
            c(clim2[2] + delta2, new_right_lim2),
            delta2,
            sigma_fun1 = sigma_fun1,
            sigma_fun2 = sigma_fun2
          )
          new_lim2[2] <- new_right_lim2
          boundary_discrete <<- rbind(boundary_discrete, boundary_extra)
        }
      }
      ## only for debugging
      ## assert_that(all(order(boundary_discrete[,1]) == 1:nrow(boundary_discrete)), msg="x grid must stay ordered!")
      boundary <<- fit_boundary_consumer(
        boundary_discrete,
        linear_boundary,
        decision
      )
      clim2 <<- new_lim2
    }

    return(boundary(y2))
  }

  decision_boundary
}

#' @templateVar fun decision2S_boundary
#' @template design2S-poisson
#' @export
decision2S_boundary.gammaMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  assert_that(likelihood(prior1) == "poisson")
  assert_that(likelihood(prior2) == "poisson")

  # only the second n2 argument may be 0
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_gammaMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  } else {
    decision2S_boundary_gammaMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  }
}

#' @keywords internal
decision2S_boundary_gammaMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  crit_lower <- decision2S_boundary_gammaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    eps,
    ...
  )
  crit_upper <- decision2S_boundary_gammaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    eps,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}

#' @keywords internal
decision2S_boundary_gammaMix_1sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  decision <- if (has_lower(decision)) {
    lower(decision)
  } else {
    upper(decision)
  }
  decision2S_boundary_gammaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    eps,
    ...
  )
}

#' @keywords internal
decision2S_boundary_gammaMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  if (!missing(eps)) {
    assert_number(eps, lower = 0, upper = 0.1, finite = TRUE)
  }

  cond_decisionStep <- function(post2) {
    function(m1) {
      decision(postmix(prior1, n = n1, m = m1 / n1), post2) - 0.25
    }
  }

  clim1 <- c(Inf, -Inf)
  clim2 <- c(Inf, -Inf)
  boundary <- NA
  grid <- NA
  lower.tail <- attr(decision, "lower.tail")

  decision_boundary <- function(y2, lim1) {
    if (missing(lim1)) {
      lambda1 <- summary(prior1, probs = c())["mean"] * n1
      lim1 <- qpois(c(eps / 2, 1 - eps / 2), lambda1)
    }

    lim2 <- range(y2)

    assert_number(lim1[1], lower = 0, finite = TRUE)
    assert_number(lim1[2], lower = 0, finite = TRUE)
    assert_number(lim2[1], lower = 0, finite = TRUE)
    assert_number(lim2[2], lower = 0, finite = TRUE)

    ## check if the boundary needs to be recomputed
    if (
      lim1[1] < clim1[1] |
        lim1[2] > clim1[2] |
        lim2[1] < clim2[1] |
        lim2[2] > clim2[2]
    ) {
      ## ensure that lim1 never shrinks in size
      lim1[1] <- min(lim1[1], clim1[1])
      lim1[2] <- max(lim1[2], clim1[2])
      grid <<- lim2[1]:lim2[2]
      Neval <- length(grid)
      boundary <<- rep(NA, Neval)
      ## The step objective is monotone in m1 in the direction set by
      ## lower.tail; the critical y1 itself is monotone in y2 for the
      ## supported one-sided decisions.
      extend_dir <- if (lower.tail) "downX" else "upX"

      decFun_for <- function(y2) {
        if (n2 == 0) {
          cond_decisionStep(prior2)
        } else {
          cond_decisionStep(postmix(prior2, n = n2, m = y2 / n2))
        }
      }
      solve_raw <- function(y2, bracket) {
        uniroot_int(
          decFun_for(y2),
          bracket,
          extendInt = extend_dir,
          clamp = lim1
        )
      }
      ## map the constant-decision sentinels into the ordered domain
      ## [-1, Inf] so they join the monotone fill: never true
      ## (-Inf -> -1), always true (Inf, no finite upper count).
      map_root <- function(root) {
        if (is.finite(root)) {
          root
        } else if (root < 0) {
          -1
        } else {
          Inf
        }
      }
      set_val <- function(y2, val) {
        boundary[y2 - lim2[1] + 1] <<- val
      }

      if (getOption("RBesT.decision2S_boundary", "adaptive") == "adaptive") {
        ## Adaptive solve replacement. Because the critical value is
        ## monotone in y2, a run of outcomes bracketed by two equal
        ## critical values is provably constant and filled without any
        ## further (expensive pmixdiff) uniroot_int calls; otherwise the
        ## midpoint is solved, bracketed by the two neighbour roots, and
        ## the halves recurse.
        bracket_for <- function(ra, rb) {
          lo <- max(min(ra, rb), lim1[1])
          hi <- min(max(ra, rb), lim1[2])
          if (lo >= hi) lim1 else c(lo, hi)
        }
        fill <- function(a, ra, b, rb) {
          if (ra == rb) {
            for (yy in a:b) {
              set_val(yy, ra)
            }
            return(invisible())
          }
          if (b - a <= 1L) {
            set_val(a, ra)
            set_val(b, rb)
            return(invisible())
          }
          m <- (a + b) %/% 2L
          rm <- map_root(solve_raw(m, bracket_for(ra, rb)))
          fill(a, ra, m, rm)
          fill(m, rm, b, rb)
        }
        ra <- map_root(solve_raw(lim2[1], lim1))
        rb <- map_root(solve_raw(lim2[2], lim1))
        fill(lim2[1], ra, lim2[2], rb)
      } else {
        ## Legacy sequential sweep: warm-start each grid point's integer
        ## search with a tight bracket around the previous boundary.
        ## Note: the loop variable must not shadow the `y2` argument of
        ## decision_boundary, which is read again by the lookup below.
        prev <- NA
        for (yy in lim2[1]:lim2[2]) {
          bracket <- if (is.na(prev)) {
            lim1
          } else {
            c(max(prev - 2, lim1[1]), min(prev + 2, lim1[2]))
          }
          root <- solve_raw(yy, bracket)
          set_val(yy, map_root(root))
          prev <- if (is.finite(root)) root else NA
        }
      }

      if (lower.tail) {
        ## if lower.tail==TRUE, then the condition becomes
        ## true when going from large to small values, hence
        ## we need to integrate from 0 to the boundary
        boundary <<- pmax(boundary - 1, -1)
      }

      ## save limits of new grid
      clim1 <<- lim1
      clim2 <<- lim2
    }

    assert_numeric(y2, lower = 0, finite = TRUE, any.missing = FALSE)
    crit <- boundary[y2 - clim2[1] + 1]
    ## in case the lower boundary of the searched grid is not
    ## zero, then we cannot say anything about cases when the
    ## decision is always negative
    if (!lower.tail) {
      ## in this case the decision changes from negative to
      ## positive when going from small to large values. Hence,
      ## if the decision is always negative, then we can be sure
      ## that the decision changes past lim1[2]. We set the
      ## boundary to lim1[2]+1 if lim1 has been given; otherwise
      ## to NA.
      ## Should the decision be negative at all values,
      ## it can change at larger values.
      if (missing(lim1)) {
        crit[crit == -1] <- NA
        crit[crit == Inf] <- NA
      } else {
        crit[crit == -1] <- lim1[2] + 1
        crit[crit == Inf] <- lim1[1] - 1
      }
    } else {
      ## now the decision changes from positive to negative
      ## when going from small to large => should the
      ## decision not change in the clim1 domain then we do
      ## not know if it happens later
      if (clim1[1] > 0) {
        crit[crit == -1] <- NA
      }
      ## however, if crit==Inf then we can be sure that the
      ## decision is indeed always positive
    }
    return(crit)
  }
  decision_boundary
}
