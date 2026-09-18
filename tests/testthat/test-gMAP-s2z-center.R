## Partial centering under RBesT.MC.ncp = 3.
##
## The checks here are target-level and MC-free: they hold the Stan program to
## the identities that make the option a *reparametrization* rather than a new
## model. Sampled behaviour is covered by the existing s2z fixtures.
##
## The decisive one is the last: partial centering must leave the marginal
## posterior of tau untouched. Getting the log-Jacobian wrong by the
## (J-1) log(tau) that the prior and the determinant have in common is the
## single most likely way to break this, and it is invisible to the boundary
## identities because it vanishes at both c = 0 and c = 1.

skip_on_cran()

#' Re-instantiate a scenario's Stan program with a different centering vector.
#'
#' `re_param` optionally overrides the resolved model parametrization, which is
#' what the boundary identities need: the two shipped parametrizations are
#' only reachable through it.
s2z_center_fit <- function(name, center, re_param = NULL) {
  d <- s2z_test_skeleton(name)$fit.data
  d$re_center <- as.array(rep_len(center, d$n_groups))
  if (is.null(re_param)) {
    if (any(d$re_center > 0)) {
      d$re_param <- 2L
    }
  } else {
    d$re_param <- as.integer(re_param)
  }
  list(fit = s2z_log_prob_fit(d), data = d)
}

#' Differences of the target across a grid, which drop the additive constant
#' Stan is free to leave out.
s2z_center_profile <- function(fit, grid) {
  lp <- vapply(grid, function(u) rstan::log_prob(fit, u), numeric(1))
  lp - lp[1]
}

#' The same, for the conventional R reference target.
conv_center_profile <- function(d, grid, drop = "none") {
  lp <- vapply(
    grid,
    function(u) gmap_conv_reference_log_prob(d, u, drop = drop),
    numeric(1)
  )
  lp - lp[1]
}

#' Laplace-exact log marginal of the target over the group coordinates.
#'
#' At fixed `tau` and `beta_raw` the gaussian-endpoint target is exactly
#' quadratic in the group coordinates, so
#'
#'   log int exp(f) = f(u*) + (n/2) log(2 pi) - 0.5 log det(-H)
#'
#' is exact rather than an approximation. The resulting tau profile is a
#' property of the *model*, so it must not depend on the parametrization. A
#' log-Jacobian that is wrong by a tau-dependent factor -- the single most
#' likely mistake, and one that cancels at both boundary values -- shows up
#' here and nowhere else.
#'
#' @param fit sampler-free stanfit
#' @param idx unconstrained slots holding the group coordinates
#' @param tau_slots unconstrained slots holding `tau_raw`
s2z_center_tau_profile <- function(fit, idx, tau_slots) {
  vapply(seq(-2, 2, length.out = 5), function(tr) {
    u0 <- rep(0, rstan::get_num_upars(fit))
    u0[tau_slots] <- tr
    f <- function(v) {
      u <- u0
      u[idx] <- v
      rstan::log_prob(fit, u)
    }
    g <- function(v) {
      u <- u0
      u[idx] <- v
      rstan::grad_log_prob(fit, u)[idx]
    }
    opt <- stats::optim(
      rep(0, length(idx)),
      function(v) -f(v),
      function(v) -g(v),
      method = "BFGS",
      control = list(reltol = 1e-14)
    )
    ## Hessian by central differences of the analytic gradient
    h <- 1e-4
    H <- vapply(seq_along(idx), function(k) {
      e <- rep(0, length(idx))
      e[k] <- h
      (g(opt$par + e) - g(opt$par - e)) / (2 * h)
    }, numeric(length(idx)))
    H <- (H + t(H)) / 2
    -opt$value +
      0.5 * length(idx) * log(2 * pi) -
      0.5 * determinant(-H, logarithm = TRUE)$modulus
  }, numeric(1))
}

## the conventional (S2Z-off) scenarios the parametrization has to serve:
## binomial with J = 8 groups, normal with J = 6, and the Student-t random
## effect -- the parametrization is a location-scale statement, so it must
## serve both re_dist values
conv_center_scenarios <- c(
  "optout", "optout_normal", "fallback_student_t",
  "fallback_tau_strata", "fallback_no_intercept"
)
conv_legacy_center_scenarios <- setdiff(
  conv_center_scenarios,
  "fallback_no_intercept"
)

test_that("ncp accepts only integer sampling policies from zero through three", {
  make_fit <- function(ncp) {
    withr::with_options(
      list(RBesT.MC.ncp = ncp),
      suppressMessages(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = RBesT::AS,
      family = binomial,
      tau.dist = "HalfNormal",
      tau.prior = 0.5,
      beta.prior = 2,
      chains = 0
      ))
    )
  }

  for (ncp in 0:3) {
    expect_s3_class(make_fit(ncp), "gMAP")
  }
  for (ncp in list(-1, 4, 1.5, NA_real_)) {
    expect_error(make_fit(ncp))
  }
})

test_that("the automatic rule returns per-group fractions in the unit interval", {
  fit <- withr::with_options(
    list(RBesT.MC.ncp = 3),
    suppressMessages(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = RBesT::AS,
      family = binomial,
      tau.dist = "HalfNormal",
      tau.prior = 0.5,
      beta.prior = 2,
      chains = 0
    ))
  )
  center <- as.numeric(fit$fit.data$re_center)

  expect_length(center, nlevels(factor(RBesT::AS$study)))
  expect_true(all(center > 0 & center < 1))
  expect_identical(fit$fit.data$re_param, 2L)
})

test_that("fixed tau gives the Fisher centering fraction directly", {
  tau <- 0.5
  fit <- withr::with_options(
    list(RBesT.MC.ncp = 3),
    suppressMessages(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = RBesT::AS,
      family = binomial,
      tau.dist = "Fixed",
      tau.prior = tau,
      beta.prior = 2,
      chains = 0
    ))
  )
  p <- (RBesT::AS$r + 0.5) / (RBesT::AS$n + 1)
  information <- RBesT::AS$n * p * (1 - p)
  expected <- tau^2 * information / (1 + tau^2 * information)

  expect_equal(as.numeric(fit$fit.data$re_center), expected)
})

test_that("working information follows each likelihood approximation", {
  gaussian <- RBesT:::.gmap_working_data(
    "gaussian", c(1, 2), c(0.5, 0.25), numeric(2), rep(1, 2),
    numeric(2), numeric(2), 1:2, 2
  )
  expect_equal(gaussian$information, c(4, 16))

  r <- c(0, 7)
  n <- c(10, 20)
  p <- (r + 0.5) / (n + 1)
  binomial <- RBesT:::.gmap_working_data(
    "binomial", numeric(2), numeric(2), r, n, numeric(2),
    numeric(2), 1:2, 2
  )
  expect_equal(binomial$information, n * p * (1 - p))

  count <- c(0, 4)
  poisson <- RBesT:::.gmap_working_data(
    "poisson", numeric(2), numeric(2), numeric(2), rep(1, 2),
    count, log(c(2, 4)), 1:2, 2
  )
  expect_equal(poisson$information, count + 0.5)
  expect_equal(poisson$response, log(count + 0.5) - log(c(2, 4)))
})

test_that("partial initialization supports every tau prior", {
  priors <- list(
    HalfNormal = 0.5,
    TruncNormal = c(0.2, 0.5),
    Uniform = c(0, 1),
    Gamma = c(2, 4),
    InvGamma = c(3, 1),
    LogNormal = c(-1, 0.5),
    TruncCauchy = c(0, 0.5),
    Exp = 2,
    Fixed = 0.5
  )

  fits <- lapply(names(priors), function(tau.dist) {
    withr::with_options(
      list(RBesT.MC.ncp = 3),
      suppressMessages(gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = RBesT::AS,
        family = binomial,
        tau.dist = tau.dist,
        tau.prior = priors[[tau.dist]],
        beta.prior = 2,
        chains = 0
      ))
    )
  })

  for (fit in fits) {
    expect_identical(fit$fit.data$re_param, 2L)
    expect_true(all(fit$fit.data$re_center > 0))
    expect_true(all(fit$fit.data$re_center < 1))
  }
})

test_that("partial initialization covers conventional model variants", {
  fit_scenario <- function(name) {
    scenario <- s2z_test_scenarios()[[name]]
    withr::with_options(
      utils::modifyList(scenario$opts, list(RBesT.MC.ncp = 3)),
      suppressMessages(suppressWarnings(do.call(
        gMAP, c(scenario$args(), list(chains = 0))
      )))
    )
  }

  for (name in c(
    "fallback_student_t", "fallback_tau_strata",
    "fallback_no_intercept", "normal_ncp", "poisson_ncp"
  )) {
    fit <- fit_scenario(name)
    expect_identical(fit$fit.data$re_param, 2L, label = name)
    expect_true(all(fit$fit.data$re_center >= 0), label = name)
    expect_true(any(fit$fit.data$re_center > 0), label = name)
    expect_true(all(fit$fit.data$re_center < 1), label = name)
  }
})

test_that("partial initialization is deterministic and prior-only safe", {
  make_fit <- function(prior_PD = FALSE) {
    withr::with_options(
      list(RBesT.MC.ncp = 3),
      suppressMessages(gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = RBesT::AS,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 0.5,
        beta.prior = 2,
        prior_PD = prior_PD,
        chains = 0
      ))
    )
  }

  first <- make_fit()
  second <- make_fit()
  expect_identical(first$fit.data$re_center, second$fit.data$re_center)

  prior <- make_fit(TRUE)
  expect_identical(prior$fit.data$re_param, 1L)
  expect_true(all(prior$fit.data$re_center == 0))
})

test_that("unused group levels do not alter observed centering fractions", {
  fit_data <- function(data) {
    withr::with_options(
      list(RBesT.MC.ncp = 3),
      suppressMessages(gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = data,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 0.5,
        beta.prior = 2,
        chains = 0
      ))
    )
  }

  baseline <- fit_data(RBesT::AS)
  expanded <- RBesT::AS
  expanded$study <- factor(
    expanded$study,
    levels = c(unique(expanded$study), "Unobserved")
  )
  with_unused <- fit_data(expanded)

  expect_equal(
    head(with_unused$fit.data$re_center, -1),
    baseline$fit.data$re_center
  )
  expect_identical(as.numeric(tail(with_unused$fit.data$re_center, 1)), 0)
})

test_that("initializer failure warns and falls back to whole-fit NCP", {
  data <- transform(
    RBesT::AS,
    stratum = factor(rep(LETTERS[1:4], each = 2))
  )
  expect_warning(
    fit <- withr::with_options(
      list(RBesT.MC.ncp = 3),
      suppressMessages(gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = data,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = rep(0.5, 4),
        tau.strata = stratum,
        beta.prior = 2,
        chains = 0
      ))
    ),
    "more than three observed non-fixed tau strata"
  )
  expect_identical(fit$fit.data$re_param, 1L)
  expect_true(all(fit$fit.data$re_center == 0))
})

test_that("centered no-intercept models use the exact partial endpoint", {
  scenario <- s2z_test_scenarios()[["fallback_no_intercept"]]
  fit <- withr::with_options(
    list(RBesT.MC.ncp = 0),
    suppressMessages(suppressWarnings(do.call(
      gMAP, c(scenario$args(), list(chains = 0))
    )))
  )

  expect_identical(fit$fit.data$re_param, 2L)
  expect_true(all(fit$fit.data$re_center == 1))
})

test_that("center = 0 reproduces the shipped non-centered parametrization", {
  ## c = 0 never consults the numeraire g (partial_center_scale() returns 1
  ## unconditionally, partial_center_loc()/_effect() drop the g term too), so
  ## this remains a pointwise match against the dedicated fast path regardless
  ## of group_scale_guess.
  for (name in c("binomial_ncp", "normal_ncp")) {
    grid <- s2z_upars_grid(name)
    ref_d <- s2z_test_skeleton(name)$fit.data
    ref_d$re_center <- as.array(rep(0, ref_d$n_groups))
    ref_d$re_param <- 1L
    ref <- s2z_center_profile(s2z_log_prob_fit(ref_d), grid)

    prop <- s2z_center_profile(s2z_center_fit(name, 0)$fit, grid)

    expect_equal(
      prop,
      ref,
      tolerance = 1e-8,
      label = paste0(name, " center = 0")
    )
  }
})

test_that("center = 1 reproduces its own reference parametrization", {
  ## Unlike c = 0, c = 1 uses group_scale_guess as the numeraire g, which is a
  ## per-group quadrature-derived vector and generally differs from the
  ## scalar beta_raw_guess[2,1] the dedicated re_param = 0 fast path uses.
  ## partial_center_scale()'s doc comment shows both are *exact*
  ## reparametrizations of the same physical N(location, tau) model for any
  ## positive g (shared or per-group), but they are no longer pointwise
  ## identical in xi_eta-space once the numeraire differs -- this checks the
  ## general path against its own reference formula (which is aware of
  ## group_scale_guess) instead of the dedicated fast path.
  for (name in c("binomial_ncp", "normal_ncp")) {
    grid <- s2z_upars_grid(name)
    f <- s2z_center_fit(name, 1)
    stan_lp <- vapply(grid, function(u) rstan::log_prob(f$fit, u), numeric(1))
    ref_lp <- vapply(
      grid,
      function(u) s2z_reference_log_prob(f$data, u),
      numeric(1)
    )
    expect_equal(
      stan_lp - stan_lp[1],
      ref_lp - ref_lp[1],
      tolerance = 1e-8,
      label = name
    )
  }
})

test_that("neutral group guesses recover the shipped centered S2Z path", {
  for (name in c("binomial_ncp", "normal_ncp")) {
    grid <- s2z_upars_grid(name)
    legacy <- s2z_test_skeleton(name)$fit.data
    legacy$re_param <- 0L
    legacy$re_center[] <- 0

    generalized <- legacy
    generalized$re_param <- 2L
    generalized$re_center[] <- 1
    generalized$group_scale_guess[] <- legacy$beta_raw_guess[2, 1]
    generalized$group_location_guess[] <- 0

    expect_equal(
      s2z_center_profile(s2z_log_prob_fit(generalized), grid),
      s2z_center_profile(s2z_log_prob_fit(legacy), grid),
      tolerance = 1e-8,
      label = name
    )
  }
})

test_that("the Stan target matches the reference target under partial centering", {
  for (name in c("binomial_ncp", "normal_ncp")) {
    grid <- s2z_upars_grid(name)
    J <- s2z_test_skeleton(name)$fit.data$n_groups
    ## a deliberately heterogeneous vector: a constant one cannot detect a
    ## missing re-projection onto the zero-sum subspace
    center <- seq(0.1, 0.9, length.out = J)

    f <- s2z_center_fit(name, center)
    stan_lp <- vapply(grid, function(u) rstan::log_prob(f$fit, u), numeric(1))
    ref_lp <- vapply(
      grid,
      function(u) s2z_reference_log_prob(f$data, u),
      numeric(1)
    )

    ## Stan drops parameter-independent constants, so only differences are
    ## comparable
    expect_equal(
      stan_lp - stan_lp[1],
      ref_lp - ref_lp[1],
      tolerance = 1e-8,
      label = name
    )
  }
})

test_that("the group effects stay exactly zero-sum for a heterogeneous fraction", {
  ## dividing a zero-sum vector by a non-constant scale leaves the subspace,
  ## so the re-projection in the Stan code is load bearing. theta minus its
  ## fixed-effect part must still sum to zero across groups.
  name <- "normal_ncp"
  d <- s2z_test_skeleton(name)$fit.data
  J <- d$n_groups
  f <- s2z_center_fit(name, seq(0.05, 0.95, length.out = J))

  for (u in s2z_upars_grid(name)) {
    pars <- rstan::constrain_pars(f$fit, u)
    ## theta = X beta_param + re; the intercept part is common, so the group
    ## effects are recovered by centering theta within the design
    re <- pars$theta - as.numeric(d$X %*% pars$beta_param)
    re_group <- tapply(re, d$group_index, mean)
    expect_lt(abs(sum(re_group - mean(re_group))), 1e-10)
  }
})

test_that("partial centering leaves the tau marginal invariant", {
  ## The single most valuable check in this file, see
  ## s2z_center_tau_profile(). A log-Jacobian that is wrong by (J-1) log(tau)
  ## -- the term the prior and the determinant have in common -- shows up here
  ## and nowhere else.
  name <- "normal_ncp"
  d <- s2z_test_skeleton(name)$fit.data
  J <- d$n_groups
  idx <- d$mX + d$n_tau_strata + seq_len(J - 1L)
  tau_slots <- d$mX + seq_len(d$n_tau_strata)

  profile <- function(fit) s2z_center_tau_profile(fit, idx, tau_slots)

  ref <- profile(s2z_center_fit(name, 0)$fit)
  for (center in list(0.5, 1, seq(0.1, 0.9, length.out = J))) {
    prop <- profile(s2z_center_fit(name, center)$fit)
    ## a constant offset is allowed (Stan drops constants); a tau-dependent
    ## one is not
    expect_equal(
      prop - prop[1],
      ref - ref[1],
      tolerance = 1e-5,
      label = paste0("center = ", paste(round(center, 2), collapse = ", "))
    )
  }
})

## ---------------------------------------------------------------------------
## Conventional (S2Z-off) partial centering
##
## The same option drives the ordinary parametrization, where the parameter
## mapping is the per-group affine map
##
##   sc_j = (1 - c_j) + c_j tau_j / g,   m_j = c_j (beta[1] - anchor) / g,
##   eps_j = tau_j (xi_j - m_j) / sc_j,   xi_j ~ D(m_j, sc_j).
##
## Here there is no zero-sum projection and hence no determinant correction:
## the location-scale density on the sampled coordinate *is* the target,
## Jacobian included. The checks below are therefore aimed at the three
## ingredients that can individually be wrong -- location, scale and the
## implied Jacobian -- and each is shown to be detectable.
## ---------------------------------------------------------------------------

test_that("the conventional layout carries J group coordinates, not J - 1", {
  ## the parametrization is indexed by group, so a test that assumed the s2z
  ## layout would write into the wrong unconstrained slots
  for (name in conv_center_scenarios) {
    skeleton <- s2z_test_skeleton(name)
    d <- skeleton$fit.data

    expect_false(s2z_active(skeleton), label = name)

    n_up <- rstan::get_num_upars(s2z_stanfit(name))
    expect_identical(
      n_up,
      as.integer(d$mX + d$n_tau_strata + d$n_groups),
      label = name
    )
    expect_length(s2z_constrained(name, rep(0, n_up))$xi_eta, d$n_groups)
  }
})

test_that("conventional c = 0 reproduces the non-centered parametrization", {
  ## non-vacuous: the all-zero vector takes the shipped fast path, while the
  ## reference always evaluates the parametrization, so this pins the fast
  ## path to c = 0 rather than comparing the fast path with itself
  for (name in conv_center_scenarios) {
    grid <- s2z_upars_grid(name)
    f <- s2z_center_fit(name, 0, re_param = 1L)

    expect_equal(
      s2z_center_profile(f$fit, grid),
      conv_center_profile(f$data, grid),
      tolerance = 1e-8,
      label = name
    )
  }
})

test_that("conventional c = 1 reproduces its own reference parametrization", {
  ## As with the s2z endpoint above, c = 1 uses group_scale_guess (per-group)
  ## rather than the scalar beta_raw_guess[2,1] the legacy re_param = 0 fast
  ## path uses, so this checks the general path against
  ## gmap_conv_reference_log_prob() (which reads group_scale_guess) instead of
  ## the dedicated fast path.
  for (name in conv_legacy_center_scenarios) {
    grid <- s2z_upars_grid(name)
    prop <- s2z_center_fit(name, 1, re_param = 2L)
    stan_lp <- vapply(
      grid, function(u) rstan::log_prob(prop$fit, u), numeric(1)
    )
    ref_lp <- vapply(
      grid, function(u) gmap_conv_reference_log_prob(prop$data, u), numeric(1)
    )
    expect_equal(
      stan_lp - stan_lp[1],
      ref_lp - ref_lp[1],
      tolerance = 1e-8,
      label = name
    )
  }
})

test_that("neutral group guesses recover the shipped centered conventional path", {
  for (name in conv_legacy_center_scenarios) {
    grid <- s2z_upars_grid(name)
    legacy <- s2z_test_skeleton(name)$fit.data
    legacy$re_param <- 0L
    legacy$re_center[] <- 0

    generalized <- legacy
    generalized$re_param <- 2L
    generalized$re_center[] <- 1
    generalized$group_scale_guess[] <- legacy$beta_raw_guess[2, 1]
    generalized$group_location_guess[] <- 0

    expect_equal(
      s2z_center_profile(s2z_log_prob_fit(generalized), grid),
      s2z_center_profile(s2z_log_prob_fit(legacy), grid),
      tolerance = 1e-8,
      label = name
    )
  }
})

test_that("heterogeneous conventional fractions match the reference target", {
  for (name in conv_center_scenarios) {
    grid <- s2z_upars_grid(name)
    J <- s2z_test_skeleton(name)$fit.data$n_groups
    ## a constant vector cannot separate the per-group scale from a global
    ## one; the leading exact zero additionally exercises a group that stays
    ## non-centered *inside* the partial branch
    center <- seq(0, 0.95, length.out = J)

    f <- s2z_center_fit(name, center)
    stan <- s2z_center_profile(f$fit, grid)

    expect_equal(
      stan,
      conv_center_profile(f$data, grid),
      tolerance = 1e-8,
      label = name
    )

    ## the comparison has teeth: dropping any single ingredient of the
    ## parametrization from the reference must break the agreement
    drops <- c("scale", "jacobian")
    if (f$data$has_intercept == 1) {
      drops <- c("loc", drops)
    }
    for (drop in drops) {
      expect_gt(
        max(abs(stan - conv_center_profile(f$data, grid, drop = drop))),
        1e-3,
        label = paste0(name, " drop = ", drop)
      )
    }
  }
})

test_that("resolved parametrizations ignore inactive centering fractions", {
  for (name in conv_center_scenarios) {
    grid <- s2z_upars_grid(name)
    J <- s2z_test_skeleton(name)$fit.data$n_groups
    center <- seq(0.1, 0.9, length.out = J)

    expect_equal(
      s2z_center_profile(
        s2z_center_fit(name, center, re_param = 1L)$fit,
        grid
      ),
      s2z_center_profile(
        s2z_center_fit(name, 0, re_param = 1L)$fit,
        grid
      ),
      tolerance = 1e-8,
      label = name
    )
  }
})

test_that("conventional partial centering leaves the tau marginal invariant", {
  ## the conventional counterpart of the s2z invariance check above; the
  ## marginal over the J group coordinates at fixed tau is a property of the
  ## model and must not move with the parametrization
  name <- "optout_normal"
  d <- s2z_test_skeleton(name)$fit.data
  J <- d$n_groups
  idx <- d$mX + d$n_tau_strata + seq_len(J)
  tau_slots <- d$mX + seq_len(d$n_tau_strata)

  profile <- function(f) s2z_center_tau_profile(f$fit, idx, tau_slots)

  ref <- profile(s2z_center_fit(name, 0, re_param = 1L))
  specs <- list(
    list(center = 0.5, re_param = 2L),
    list(center = 1, re_param = 2L),
    list(center = seq(0, 0.9, length.out = J), re_param = 2L),
    ## the legacy centered program, which is the same model in yet another
    ## parametrization
    list(center = 0, re_param = 0L)
  )
  for (spec in specs) {
    prop <- profile(
      s2z_center_fit(name, spec$center, re_param = spec$re_param)
    )
    expect_equal(
      prop - prop[1],
      ref - ref[1],
      tolerance = 1e-5,
      label = paste0(
        "center = ",
        paste(round(spec$center, 2), collapse = ", "),
        ", re_param = ",
        spec$re_param
      )
    )
  }
})

test_that("automatic centering is the default with s2z on and off", {
  gmap_AS <- function(opts) {
    withr::with_options(
      opts,
      suppressMessages(gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = RBesT::AS,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 0.5,
        beta.prior = 2,
        chains = 0
      ))
    )
  }

  off_auto <- gmap_AS(list(RBesT.MC.s2z = FALSE, RBesT.MC.ncp = 3))
  center <- as.numeric(off_auto$fit.data$re_center)

  expect_false(s2z_active(off_auto))
  expect_identical(off_auto$fit.data$use_s2z, 0L)
  expect_length(center, nlevels(factor(RBesT::AS$study)))
  expect_true(all(center > 0 & center < 1))
  ## the studies differ in size, so the fractions must too -- a constant
  ## vector would mean the rule never saw the per-group information
  expect_gt(diff(range(center)), 0.05)

  ## the rule is a statement about the model, not about the parametrization,
  ## so it must return what the s2z path derives from the same data
  on_auto <- gmap_AS(list(RBesT.MC.ncp = 3))
  expect_equal(center, as.numeric(on_auto$fit.data$re_center))
  expect_true(all(as.numeric(on_auto$fit.data$re_center) > 0))

  defaults <- gmap_AS(list(RBesT.MC.s2z = NULL, RBesT.MC.ncp = NULL))
  expect_true(s2z_active(defaults))
  expect_equal(center, as.numeric(defaults$fit.data$re_center))

  ## Both switches recover the pre-1.12 non-centered sampling scheme.
  expect_true(all(
    as.numeric(gmap_AS(list(
      RBesT.MC.s2z = FALSE,
      RBesT.MC.ncp = 1
    ))$fit.data$re_center) == 0
  ))

  ## and the data list this path produces really drives the conventional
  ## partial parametrization -- not a hand-patched one, and not the fast path
  d <- off_auto$fit.data
  fit <- s2z_log_prob_fit(d)
  grid <- gmap_upars_grid(d, fit)
  expect_equal(
    s2z_center_profile(fit, grid),
    conv_center_profile(d, grid),
    tolerance = 1e-8
  )
})
