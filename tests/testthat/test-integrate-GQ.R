## Tests for Gauss-Hermite quadrature integration

# --- Analytical integral tests -----------------------------------------------

test_that("GH integrates constant g(x) = 1 to 1", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  result <- RBesT:::integrate_density(mix, function(x) rep(1, length(x)))
  expect_equal(result, 1, tolerance = 1e-12)
})

test_that("GH integrates g(x) = x to mixture mean", {
  mix <- mixnorm(c(0.3, -1, 0.5), c(0.7, 2, 1.5), sigma = 1)
  expected_mean <- 0.3 * (-1) + 0.7 * 2
  result <- RBesT:::integrate_density(mix, function(x) x)
  expect_equal(result, expected_mean, tolerance = 1e-12)
})

test_that("GH integrates g(x) = x^2 to mixture second moment", {
  mix <- mixnorm(c(0.4, 1, 2), c(0.6, -1, 0.5), sigma = 1)
  # E[X^2] = sum w_k * (mu_k^2 + sigma_k^2)
  expected <- 0.4 * (1^2 + 2^2) + 0.6 * ((-1)^2 + 0.5^2)
  result <- RBesT:::integrate_density(mix, function(x) x^2)
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("GH integrates polynomial exactly for degree < 2N-1", {
  mix <- mixnorm(c(1, 3, 2), sigma = 1)
  # g(x) = x^4: E[X^4] for N(3, 4) = mu^4 + 6*mu^2*sigma^2 + 3*sigma^4
  mu <- 3
  s <- 2
  expected <- mu^4 + 6 * mu^2 * s^2 + 3 * s^4
  result <- RBesT:::integrate_density(mix, function(x) x^4)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("GH integrates Gaussian CDF against single-component normMix", {
  # int Phi((x - a) / b) * phi(x; mu, sigma) dx
  # = Phi((mu - a) / sqrt(b^2 + sigma^2))
  mu <- 1
  sigma <- 2
  a <- 0.5
  b <- 1.5
  mix <- mixnorm(c(1, mu, sigma), sigma = 1)
  expected <- pnorm((mu - a) / sqrt(b^2 + sigma^2))
  result <- RBesT:::integrate_density(
    mix, function(x) pnorm((x - a) / b)
  )
  expect_equal(result, expected, tolerance = 1e-10)
  # For mixture: sum w_k * Phi((mu_k - a) / sqrt(b^2 + sigma_k^2))
  mix <- mixnorm(c(0.3, -1, 0.5), c(0.7, 2, 1.5), sigma = 1)
  a <- 0
  b <- 1
  expected <- 0.3 * pnorm((-1 - a) / sqrt(b^2 + 0.5^2)) +
    0.7 * pnorm((2 - a) / sqrt(b^2 + 1.5^2))
  result <- RBesT:::integrate_density(
    mix, function(x) pnorm((x - a) / b)
  )
  expect_equal(result, expected, tolerance = 1e-9)
})

# --- Log-space variant -------------------------------------------------------

test_that("GH log variant matches natural-scale variant", {
  mix <- mixnorm(c(0.4, 0, 1), c(0.6, 3, 2), sigma = 1)
  g <- function(x) pnorm(x, 1, 2)
  log_g <- function(x) pnorm(x, 1, 2, log.p = TRUE)

  result_natural <- RBesT:::integrate_density(mix, g)
  result_log <- RBesT:::integrate_density_log(mix, log_g)
  expect_equal(result_log, result_natural, tolerance = 1e-12)
})

# --- Agreement with adaptive method -----------------------------------------

test_that("GH agrees with adaptive for oc2S-style integrand", {
  theta1 <- 0
  theta2 <- 0.5
  sem1 <- 1 / sqrt(10)
  sem2 <- 1 / sqrt(20)
  mix <- mixnorm(c(1, theta2, sem2), sigma = sem2)

  # Simulate an oc2S-like integrand: pnorm of a linear boundary
  log_g <- function(x) {
    pnorm(0.3 * x + 0.1, theta1, sem1, lower.tail = FALSE, log.p = TRUE)
  }

  result_gh <- RBesT:::integrate_density_log(mix, log_g)

  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result_adaptive <- RBesT:::integrate_density_log(mix, log_g)
  })

  expect_equal(result_gh, result_adaptive, tolerance = 1e-8)
})

test_that("GH agrees with adaptive for multi-component normMix", {
  mix <- mixnorm(c(0.3, -1, 0.5), c(0.7, 2, 1.5), sigma = 1)

  log_g <- function(x) pnorm(x, 0, 1, log.p = TRUE)

  result_gh <- RBesT:::integrate_density_log(mix, log_g)

  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result_adaptive <- RBesT:::integrate_density_log(mix, log_g)
  })

  expect_equal(result_gh, result_adaptive, tolerance = 1e-8)
})

test_that("GH agrees with adaptive for natural-scale integration", {
  mix <- mixnorm(c(0.5, 0, 1), c(0.5, 3, 0.5), sigma = 1)

  g <- function(x) dnorm(x, 1, 2)

  result_gh <- RBesT:::integrate_density(mix, g)

  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result_adaptive <- RBesT:::integrate_density(mix, g)
  })

  expect_equal(result_gh, result_adaptive, tolerance = 1e-8)
})

# --- Determinism -------------------------------------------------------------

test_that("GH produces bitwise identical results across repeated calls", {
  mix <- mixnorm(c(0.3, -1, 0.5), c(0.7, 2, 1.5), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)

  results <- replicate(100, RBesT:::integrate_density(mix, g))
  expect_true(all(results == results[1]))
})

# --- Edge cases --------------------------------------------------------------

test_that("GH works with single-component mixture", {
  mix <- mixnorm(c(1, 5, 0.1), sigma = 1)
  result <- RBesT:::integrate_density(mix, function(x) rep(1, length(x)))
  expect_equal(result, 1, tolerance = 1e-12)
})

test_that("GH works with very narrow component", {
  mix <- mixnorm(c(1, 0, 0.001), sigma = 1)
  # g(x) = x should give mean = 0
  result <- RBesT:::integrate_density(mix, function(x) x)
  expect_equal(result, 0, tolerance = 1e-12)
})

test_that("GH works with very wide component", {
  mix <- mixnorm(c(1, 0, 100), sigma = 1)
  result <- RBesT:::integrate_density(mix, function(x) rep(1, length(x)))
  expect_equal(result, 1, tolerance = 1e-10)
})

test_that("GH works with small weight component", {
  mix <- mixnorm(c(0.001, -10, 0.1), c(0.999, 0, 1), sigma = 1)
  expected_mean <- 0.001 * (-10) + 0.999 * 0
  result <- RBesT:::integrate_density(mix, function(x) x)
  expect_equal(result, expected_mean, tolerance = 1e-12)
})

# --- Option switch -----------------------------------------------------------

test_that("RBesT.integrate_method = 'adaptive' uses adaptive path", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)

  result_gh <- RBesT:::integrate_density(mix, g)

  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result_adaptive <- RBesT:::integrate_density(mix, g)
  })

  # Both should give very close results
  expect_equal(result_gh, result_adaptive, tolerance = 1e-8)
})

test_that("RBesT.GQ_nodes option controls node count", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)

  # Disable refinement to test the fixed-node behaviour directly.
  # N=10 should still be reasonable for this smooth integrand
  withr::with_options(list(RBesT.GQ_rel_tol = Inf, RBesT.GQ_nodes = 10L), {
    result_10 <- RBesT:::integrate_density(mix, g)
  })

  withr::with_options(list(RBesT.GQ_rel_tol = Inf, RBesT.GQ_nodes = 50L), {
    result_50 <- RBesT:::integrate_density(mix, g)
  })

  # Both should agree closely for smooth integrand
  expect_equal(result_10, result_50, tolerance = 1e-10)
})

# --- Convergence study -------------------------------------------------------

test_that("GH converges with increasing nodes", {
  mix <- mixnorm(c(0.5, -1, 0.5), c(0.5, 1, 2), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)

  # Disable refinement so each node count is used exactly as set.
  # Reference at high N
  withr::with_options(list(RBesT.GQ_rel_tol = Inf, RBesT.GQ_nodes = 50L), {
    ref <- RBesT:::integrate_density(mix, g)
  })

  errors <- vapply(c(10L, 20L, 30L, 40L), function(n) {
    withr::with_options(list(RBesT.GQ_rel_tol = Inf, RBesT.GQ_nodes = n), {
      abs(RBesT:::integrate_density(mix, g) - ref)
    })
  }, numeric(1))

  # Overall trend: error at N=30 much smaller than at N=10
  expect_lt(errors[3], errors[1])
  # N=30 should be very close to N=50 reference
  expect_lt(errors[3], 1e-6)
})

# --- betaMix (Gauss-Jacobi) tests --------------------------------------------

test_that("GJ integrates exp(log_g) to 1 for betaMix with log_g = 0", {
  mix <- mixbeta(c(1, 3, 7))
  result <- RBesT:::integrate_density_log(mix, function(x) rep(0, length(x)))
  expect_equal(result, 1, tolerance = 1e-12)
})

test_that("GJ integrates log_g = log(x) to betaMix mean", {
  mix <- mixbeta(c(0.4, 3, 7), c(0.6, 10, 5))
  expected_mean <- 0.4 * 3 / 10 + 0.6 * 10 / 15
  result <- RBesT:::integrate_density_log(mix, function(x) log(x))
  expect_equal(result, expected_mean, tolerance = 1e-12)
})

test_that("GJ integrates log_g = log(x^2) to betaMix second moment", {
  mix <- mixbeta(c(1, 5, 3))
  a <- 5
  b <- 3
  # E[X^2] = a*(a+1) / ((a+b)*(a+b+1))
  expected <- a * (a + 1) / ((a + b) * (a + b + 1))
  result <- RBesT:::integrate_density_log(mix, function(x) log(x^2))
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("GJ agrees with adaptive for betaMix CDF integrand", {
  mix <- mixbeta(c(0.3, 2, 8), c(0.7, 15, 3))
  log_g <- function(x) pnorm(x, 0.4, 0.2, log.p = TRUE)

  result_gq <- RBesT:::integrate_density_log(mix, log_g)
  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result_adaptive <- RBesT:::integrate_density_log(mix, log_g)
  })
  expect_equal(result_gq, result_adaptive, tolerance = 1e-8)
})

test_that("GJ works with single-component betaMix", {
  mix <- mixbeta(c(1, 1, 1))
  # Uniform[0,1]: E[X] = 0.5
  result <- RBesT:::integrate_density_log(mix, function(x) log(x))
  expect_equal(result, 0.5, tolerance = 1e-12)
})

test_that("GJ option switch falls through to adaptive", {
  mix <- mixbeta(c(1, 5, 5))
  log_g <- function(x) log(x^3)
  # E[X^3] for Beta(5,5) = 5*6*7/(10*11*12)
  expected <- 5 * 6 * 7 / (10 * 11 * 12)

  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result <- RBesT:::integrate_density_log(mix, log_g)
  })
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("GJ falls through to adaptive for non-identity dlink", {
  mix <- mixbeta(c(1, 5, 5))
  RBesT:::dlink(mix) <- RBesT:::logit_dlink
  log_g <- function(x) pnorm(x, 0, 1, log.p = TRUE)

  # Should not error — falls through to adaptive
  result <- RBesT:::integrate_density_log(mix, log_g)
  expect_true(is.finite(result))
})

# --- gammaMix (Gauss-Laguerre) tests -----------------------------------------

test_that("GL integrates exp(log_g) to 1 for gammaMix with log_g = 0", {
  mix <- mixgamma(c(1, 20, 4))
  result <- RBesT:::integrate_density_log(mix, function(x) rep(0, length(x)))
  expect_equal(result, 1, tolerance = 1e-10)
})

test_that("GL integrates log_g = log(x) to gammaMix mean", {
  mix <- mixgamma(c(0.3, 20, 4), c(0.7, 50, 10))
  expected_mean <- 0.3 * 20 / 4 + 0.7 * 50 / 10
  result <- RBesT:::integrate_density_log(mix, function(x) log(x))
  expect_equal(result, expected_mean, tolerance = 1e-10)
})

test_that("GL integrates log_g = log(x^2) to gammaMix second moment", {
  mix <- mixgamma(c(1, 10, 2))
  shape <- 10
  rate <- 2
  # E[X^2] = shape*(shape+1) / rate^2
  expected <- shape * (shape + 1) / rate^2
  result <- RBesT:::integrate_density_log(mix, function(x) log(x^2))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("GL agrees with adaptive for gammaMix CDF integrand", {
  mix <- mixgamma(c(0.4, 15, 3), c(0.6, 40, 8))
  log_g <- function(x) pnorm(x, 5, 1, log.p = TRUE)

  result_gq <- RBesT:::integrate_density_log(mix, log_g)
  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result_adaptive <- RBesT:::integrate_density_log(mix, log_g)
  })
  expect_equal(result_gq, result_adaptive, tolerance = 1e-6)
})

test_that("GL works with single-component gammaMix", {
  mix <- mixgamma(c(1, 5, 2))
  # E[X] = shape/rate = 5/2
  result <- RBesT:::integrate_density_log(mix, function(x) log(x))
  expect_equal(result, 2.5, tolerance = 1e-10)
})

test_that("GL option switch falls through to adaptive", {
  mix <- mixgamma(c(1, 10, 2))
  log_g <- function(x) log(x)
  expected <- 10 / 2

  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    result <- RBesT:::integrate_density_log(mix, log_g)
  })
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("GL handles large shape parameters without overflow", {
  mix <- mixgamma(c(1, 172, 102))
  log_g <- function(x) pgamma(x, 2, 102, log.p = TRUE)

  result <- RBesT:::integrate_density_log(mix, log_g)
  expect_true(is.finite(result))
  expect_equal(result, 1, tolerance = 1e-8)
})

test_that("GL falls through to adaptive for non-identity dlink", {
  mix <- mixgamma(c(1, 20, 4))
  RBesT:::dlink(mix) <- RBesT:::log_dlink
  log_g <- function(x) pnorm(x, 1, 1, log.p = TRUE)

  # Should not error — falls through to adaptive
  result <- RBesT:::integrate_density_log(mix, log_g)
  expect_true(is.finite(result))
})

# --- Tolerance control (refinement driver) -----------------------------------

test_that("legacy single-shot path (rel_tol = Inf) matches fixed-node GH", {
  mix <- mixnorm(c(0.5, -1, 0.5), c(0.5, 1, 2), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)

  withr::with_options(list(RBesT.GQ_rel_tol = Inf, RBesT.GQ_nodes = 30L), {
    via_driver <- RBesT:::integrate_density(mix, g)
  })

  gh <- statmod::gauss.quad(30L, kind = "hermite")
  total <- 0
  for (k in 1:2) {
    x_k <- mix[2, k] + sqrt(2) * mix[3, k] * gh$nodes
    total <- total + mix[1, k] / sqrt(pi) * sum(gh$weights * g(x_k))
  }
  expect_identical(via_driver, total)
})

test_that("refinement reaches requested tolerance vs high-accuracy reference", {
  mix <- mixnorm(c(0.5, -1, 0.5), c(0.5, 1, 2), sigma = 1)
  # Moderately stiff integrand: a sigmoid narrower than the components but
  # still resolvable by refinement before the node cap.
  g <- function(x) pnorm(x, 0.3, 0.25)

  withr::with_options(list(RBesT.GQ_rel_tol = Inf, RBesT.GQ_nodes = 400L), {
    ref <- RBesT:::integrate_density(mix, g)
  })

  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-4, RBesT.GQ_abs_tol = 1e-6,
         RBesT.GQ_nodes = 15L, RBesT.GQ_max_nodes = 400L),
    {
      val <- RBesT:::integrate_density(mix, g)
    }
  )
  expect_lt(abs(val - ref), max(1e-6, 1e-4 * abs(ref)))
})

test_that("smooth integrand converges within the first refinement step", {
  calls <- 0L
  trace_g <- function(x) {
    calls <<- calls + 1L
    pnorm(x, 0, 1)
  }
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-4, RBesT.GQ_abs_tol = 1e-6,
         RBesT.GQ_nodes = 30L, RBesT.GQ_node_growth = 2),
    {
      RBesT:::integrate_density(mix, trace_g)
    }
  )
  # n0 evaluation + a single doubling step => 2 integrand evaluations
  expect_equal(calls, 2L)
})

test_that("on_nonconvergence = 'warn' warns and returns best estimate", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)
  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-16, RBesT.GQ_abs_tol = 0,
         RBesT.GQ_nodes = 10L, RBesT.GQ_max_nodes = 20L,
         RBesT.GQ_on_nonconvergence = "warn"),
    {
      expect_warning(
        val <- RBesT:::integrate_density(mix, g),
        "did not reach the requested tolerance"
      )
    }
  )
  expect_true(is.finite(val))
})

test_that("on_nonconvergence = 'error' throws at the node cap", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)
  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-16, RBesT.GQ_abs_tol = 0,
         RBesT.GQ_nodes = 10L, RBesT.GQ_max_nodes = 20L,
         RBesT.GQ_on_nonconvergence = "error"),
    {
      expect_error(
        RBesT:::integrate_density(mix, g),
        "did not reach the requested tolerance"
      )
    }
  )
})

test_that("on_nonconvergence = 'silent' returns best estimate quietly", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)
  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-16, RBesT.GQ_abs_tol = 0,
         RBesT.GQ_nodes = 10L, RBesT.GQ_max_nodes = 20L,
         RBesT.GQ_on_nonconvergence = "silent"),
    {
      expect_silent(val <- RBesT:::integrate_density(mix, g))
    }
  )
  expect_true(is.finite(val))
})

test_that("on_nonconvergence = 'adaptive' falls back to adaptive (betaMix)", {
  mix <- mixbeta(c(0.6, 3, 7), c(0.4, 8, 4))
  log_g <- function(x) pnorm(x, 0.4, 0.05, log.p = TRUE)

  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    adaptive_val <- RBesT:::integrate_density_log(mix, log_g)
  })

  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-16, RBesT.GQ_abs_tol = 0,
         RBesT.GQ_nodes = 10L, RBesT.GQ_max_nodes = 20L,
         RBesT.GQ_on_nonconvergence = "adaptive"),
    {
      val <- RBesT:::integrate_density_log(mix, log_g)
    }
  )
  expect_equal(val, adaptive_val, tolerance = 1e-8)
})

test_that("refinement is deterministic across repeated calls", {
  mix <- mixnorm(c(0.4, 0.2, 0.7), c(0.6, -0.5, 1.3), sigma = 1)
  g <- function(x) pnorm(x, 0.1, 0.3)
  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-4, RBesT.GQ_abs_tol = 1e-6),
    {
      a <- RBesT:::integrate_density(mix, g)
      b <- RBesT:::integrate_density(mix, g)
    }
  )
  expect_identical(a, b)
})

test_that("abs_tol floor handles a zero-valued integral", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  g <- function(x) rep(0, length(x))
  withr::with_options(
    list(RBesT.GQ_rel_tol = 1e-4, RBesT.GQ_abs_tol = 1e-6),
    {
      val <- RBesT:::integrate_density(mix, g)
    }
  )
  expect_equal(val, 0, tolerance = 1e-12)
})

test_that("betaMix and gammaMix agree with adaptive under refinement", {
  bmix <- mixbeta(c(0.5, 10, 30), c(0.5, 20, 10))
  log_gb <- function(x) pnorm(x, 0.3, 0.1, log.p = TRUE)
  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    b_ad <- RBesT:::integrate_density_log(bmix, log_gb)
  })
  withr::with_options(list(RBesT.GQ_rel_tol = 1e-4, RBesT.GQ_abs_tol = 1e-6), {
    b_gq <- RBesT:::integrate_density_log(bmix, log_gb)
  })
  expect_equal(b_gq, b_ad, tolerance = 1e-4)

  gmix <- mixgamma(c(0.7, 30, 2), c(0.3, 20, 1))
  log_gg <- function(x) pgamma(x, 5, 1, log.p = TRUE)
  withr::with_options(list(RBesT.integrate_method = "adaptive"), {
    g_ad <- RBesT:::integrate_density_log(gmix, log_gg)
  })
  withr::with_options(list(RBesT.GQ_rel_tol = 1e-4, RBesT.GQ_abs_tol = 1e-6), {
    g_gq <- RBesT:::integrate_density_log(gmix, log_gg)
  })
  expect_equal(g_gq, g_ad, tolerance = 1e-4)
})

# --- Argument validation (checkmate) -----------------------------------------

test_that("GQ refinement validates user-provided options", {
  mix <- mixnorm(c(1, 0, 1), sigma = 1)
  g <- function(x) pnorm(x, 0, 1)

  # n0 must be >= 10
  withr::with_options(list(RBesT.GQ_nodes = 9L), {
    expect_error(RBesT:::integrate_density(mix, g))
  })
  # n0 must be < max_n
  withr::with_options(list(RBesT.GQ_nodes = 50L, RBesT.GQ_max_nodes = 40L), {
    expect_error(RBesT:::integrate_density(mix, g))
  })
  # rel_tol must be >= 0
  withr::with_options(list(RBesT.GQ_rel_tol = -1), {
    expect_error(RBesT:::integrate_density(mix, g))
  })
  # abs_tol must be >= 0
  withr::with_options(list(RBesT.GQ_abs_tol = -1e-6), {
    expect_error(RBesT:::integrate_density(mix, g))
  })
  # growth must be > 1
  withr::with_options(list(RBesT.GQ_node_growth = 1), {
    expect_error(RBesT:::integrate_density(mix, g))
  })
  # on_fail must be a known choice
  withr::with_options(list(RBesT.GQ_on_nonconvergence = "bogus"), {
    expect_error(RBesT:::integrate_density(mix, g))
  })
})
