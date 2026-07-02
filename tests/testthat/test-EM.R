## test that the EM algorithms recover reliably test distributions;
## test criterium is a "sufficiently" small KL divergence

## the test considers the uni-variate normal, beta & gamma case in
## three flavours each:

## - single-component
## - two component mixture with heavy tails
## - three component mixture with bi-modal density and heavy tails

EM_settings <- function() {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    list(Nsim = 1e4, verbose = FALSE, KLthresh = 1e-2)
  } else {
    list(Nsim = 1e3, verbose = FALSE, KLthresh = 1e-1)
  }
}

EM_covariance <- function() {
  p <- 4
  Rho <- diag(p)
  Rho[lower.tri(Rho)] <- c(0.3, 0.8, -0.2, 0.1, 0.5, -0.4)
  Rho[upper.tri(Rho)] <- t(Rho)[upper.tri(Rho)]
  s <- c(1, 2, 3, 4)
  diag(s, p) %*% Rho %*% diag(s, p)
}

EM_reference <- function(name) {
  S <- EM_covariance()
  zero <- rep(0, nrow(S))
  switch(
    name,
    norm_single = mixnorm(c(1, 1, 5), param = "mn", sigma = 1),
    norm_heavy = mixnorm(
      c(0.5, 0, 0.25),
      c(0.5, 1, 5),
      param = "mn",
      sigma = 1
    ),
    norm_bi = mixnorm(
      c(0.5, 0, 0.5),
      c(0.25, 1, 5),
      c(0.25, -1, 2),
      param = "mn",
      sigma = 1
    ),
    mvnorm_single = mixmvnorm(c(1, zero, 5), param = "mn", sigma = S),
    mvnorm_heavy = mixmvnorm(
      c(0.5, zero, 0.25),
      c(0.5, zero + 1, 5),
      param = "mn",
      sigma = S
    ),
    mvnorm_bi = mixmvnorm(
      c(0.5, zero, 0.5),
      c(0.25, zero + 1, 5),
      c(0.25, zero - 1, 2),
      param = "mn",
      sigma = S
    ),
    mvnorm_bi_1D = mixmvnorm(
      c(0.5, 0, 0.5),
      c(0.25, 1, 5),
      c(0.25, -1, 2),
      param = "mn",
      sigma = S[1, 1, drop = FALSE]
    ),
    beta_single = mixbeta(c(1, 0.3, 10), param = "mn"),
    ## This density is challenging for constrained beta EM and leads to a large KLdiv.
    beta_single_alt = mixbeta(c(1, 0.2, 3)),
    beta_heavy = mixbeta(c(0.8, 0.3, 10), c(0.2, 0.5, 2.5), param = "mn"),
    beta_bi = mixbeta(
      c(0.3, 0.3, 20),
      c(0.2, 0.5, 2),
      c(0.5, 0.7, 10),
      param = "mn"
    ),
    gamma_single = mixgamma(c(1, 7.5, 5), param = "mn", likelihood = "poisson"),
    gamma_heavy = mixgamma(
      c(0.5, 7.5, 0.5),
      c(0.5, 5, 10),
      param = "mn",
      likelihood = "poisson"
    ),
    gamma_bi = mixgamma(
      c(0.5, 7.5, 1),
      c(0.25, 15, 15),
      c(0.25, 5, 10),
      param = "mn",
      likelihood = "poisson"
    )
  )
}

EM_test <- function(mixTest, seed, Nsim = 1e4, verbose = FALSE, ...) {
  KLthresh <- EM_settings()$KLthresh
  samp <- withr::with_seed(seed, rmix(mixTest, Nsim))
  withr::local_seed(seed)
  EMmix1 <- mixfit(
    samp,
    type = switch(
      class(mixTest)[1],
      gammaMix = "gamma",
      normMix = "norm",
      betaMix = "beta",
      mvnormMix = "mvnorm"
    ),
    thin = 1,
    eps = 2,
    Nc = ncol(mixTest),
    verbose = verbose,
    ...
  )
  kl1 <- abs(KLdivmix(mixTest, EMmix1))
  expect_true(kl1 < KLthresh)
  ## results must not depend on the seed, but only on the order of
  ## the input sample
  EMmix2 <- withr::with_seed(
    seed + 657858,
    mixfit(
      samp,
      type = switch(
        class(mixTest)[1],
        gammaMix = "gamma",
        normMix = "norm",
        betaMix = "beta",
        mvnormMix = "mvnorm"
      ),
      thin = 1,
      eps = 2,
      Nc = ncol(mixTest),
      verbose = verbose,
      ...
    )
  )
  expect_true(
    all(EMmix1 == EMmix2),
    info = "Result of EM is independent of random seed."
  )
}

EM_mvn_test <- function(mixTest, seed, Nsim = 1e4, verbose = FALSE, ...) {
  samp <- withr::with_seed(seed, rmix(mixTest, Nsim))
  withr::local_seed(seed)
  EMmix1 <- mixfit(
    samp,
    type = "mvnorm",
    thin = 1,
    eps = 2,
    Nc = ncol(mixTest),
    verbose = verbose,
    ...
  )
  expect_equal(summary(mixTest)$mean, summary(EMmix1)$mean, tolerance = 0.1)
  expect_equal(summary(mixTest)$cov, summary(EMmix1)$cov, tolerance = 0.1)
  expect_equal(likelihood(EMmix1), likelihood(mixTest))
  EMmix2 <- withr::with_seed(
    seed + 476767,
    mixfit(
      samp,
      type = "mvnorm",
      thin = 1,
      eps = 2,
      Nc = ncol(mixTest),
      verbose = verbose,
      ...
    )
  )
  expect_true(
    all(EMmix1 == EMmix2),
    info = "Result of EM is independent of random seed."
  )
}

test_that("Normal EM fits single component", {
  settings <- EM_settings()
  EM_test(EM_reference("norm_single"), 3453563, settings$Nsim, settings$verbose)
})
test_that("Normal EM fits heavy-tailed mixture", {
  settings <- EM_settings()
  EM_test(EM_reference("norm_heavy"), 9275624, settings$Nsim, settings$verbose)
})
test_that("Normal EM fits bi-modal mixture", {
  settings <- EM_settings()
  EM_test(EM_reference("norm_bi"), 9345726, settings$Nsim, settings$verbose)
})

test_that("Multivariate Normal EM fits single component", {
  settings <- EM_settings()
  EM_mvn_test(EM_reference("mvnorm_single"), 3453563, max(1E4, settings$Nsim), settings$verbose)
})
test_that("Multivariate Normal EM fits heavy-tailed mixture", {
  settings <- EM_settings()
  EM_mvn_test(EM_reference("mvnorm_heavy"), 9275624, max(1E4, settings$Nsim), settings$verbose)
})
test_that("Multivariate Normal EM fits bi-modal mixture", {
  settings <- EM_settings()
  EM_mvn_test(EM_reference("mvnorm_bi"), 9345726, max(1E4, settings$Nsim), settings$verbose)
})
test_that("Multivariate Normal EM fits bi-modal mixture 1D", {
  settings <- EM_settings()
  EM_mvn_test(EM_reference("mvnorm_bi_1D"), 9345726, max(1E4, settings$Nsim), settings$verbose)
})

test_that("Gamma EM fits single component", {
  settings <- EM_settings()
  EM_test(EM_reference("gamma_single"), 9345835, settings$Nsim, settings$verbose)
})
test_that("Gamma EM fits heavy-tailed mixture", {
  settings <- EM_settings()
  EM_test(EM_reference("gamma_heavy"), 5629389, settings$Nsim, settings$verbose)
})
test_that("Gamma EM fits bi-modal mixture", {
  settings <- EM_settings()
  EM_test(EM_reference("gamma_bi"), 9373515, settings$Nsim, settings$verbose)
})

test_that("Beta EM fits single component", {
  settings <- EM_settings()
  EM_test(EM_reference("beta_single"), 7265355, settings$Nsim, settings$verbose)
})
test_that("Beta EM fits single component with mass at boundary", {
  settings <- EM_settings()
  EM_test(EM_reference("beta_single_alt"), 7265355, settings$Nsim, settings$verbose, constrain_gt1 = FALSE)
})
test_that("Beta EM fits heavy-tailed mixture", {
  settings <- EM_settings()
  EM_test(EM_reference("beta_heavy"), 2946562, settings$Nsim, settings$verbose)
})
test_that("Beta EM fits bi-modal mixture", {
  settings <- EM_settings()
  EM_test(EM_reference("beta_bi"), 9460370, settings$Nsim, settings$verbose)
})

test_that("Constrained Beta EM respects a>1 & b>1", {
  settings <- EM_settings()
  unconstrained <- mixbeta(c(0.6, 2.8, 64), c(0.25, 0.5, 0.92), c(0.15, 3, 15))
  samp <- withr::with_seed(45747, rmix(unconstrained, settings$Nsim))
  constrained <- mixfit(samp, type = "beta", Nc = 3, constrain_gt1 = TRUE)
  expect_numeric(constrained[2, ], lower = 1, any.missing = FALSE, len = 3)
  expect_numeric(constrained[3, ], lower = 1, any.missing = FALSE, len = 3)
})

test_that("EM plotting functions do not use deprecated ggplot2 functionality", {
  settings <- EM_settings()
  unconstrained <- mixbeta(c(0.6, 2.8, 64), c(0.25, 0.5, 0.92), c(0.15, 3, 15))
  samp <- withr::with_seed(45747, rmix(unconstrained, settings$Nsim))
  constrained <- mixfit(samp, type = "beta", Nc = 3, constrain_gt1 = TRUE)
  withr::with_options(
    list(lifecycle_verbosity = "error", RBesT.verbose = TRUE),
    {
      expect_no_error(plot(constrained))
    }
  )
})
