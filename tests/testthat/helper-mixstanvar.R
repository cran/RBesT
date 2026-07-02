mixstanvar_tolerance <- function() 1E-1

mixstanvar_univariate_mixes <- function() {
  list(
    beta = mixbeta(c(1, 11, 4)),
    betaMix = mixbeta(c(0.8, 11, 4), c(0.2, 1, 1)),
    gamma = mixgamma(c(1, 5, 10), param = "mn"),
    gammaMix = mixgamma(
      rob = c(0.25, 8, 0.5),
      inf = c(0.75, 8, 10),
      param = "mn"
    ),
    norm = mixnorm(c(1, 0, sqrt(2)), sigma = 1),
    normMix = mixnorm(c(0.2, 0, 1), c(0.8, 1, 1), sigma = 1)
  )
}

mixstanvar_multivariate_mixes <- function(dimension) {
  checkmate::assert_choice(dimension, c(2, 4))

  p <- 4
  Rho <- diag(p)
  Rho[lower.tri(Rho)] <- c(0.3, 0.8, -0.2, 0.1, 0.5, -0.4)
  Rho[upper.tri(Rho)] <- t(Rho)[upper.tri(Rho)]
  s <- c(1, 2, 3, 4)
  S <- diag(s, p) %*% Rho %*% diag(s, p)
  m1 <- 0:3
  m2 <- 1:4

  if (dimension == 4) {
    return(list(
      single = mixmvnorm(c(1, m1, 5), param = "mn", sigma = S),
      heavy = mixmvnorm(
        c(0.5, m1, 0.25),
        c(0.5, m2, 5),
        param = "mn",
        sigma = S
      )
    ))
  }

  list(
    single = mixmvnorm(
      c(1, m1[1:2], 5),
      param = "mn",
      sigma = S[1:2, 1:2]
    ),
    heavy = mixmvnorm(
      c(0.5, m1[1:2], 0.25),
      c(0.5, m2[1:2], 5),
      param = "mn",
      sigma = S[1:2, 1:2]
    )
  )
}

brms_sampling_args <- function() list(chains = 1, iter = 20000, init = 0.25)

brms_test_options <- function() {
  opts <- list()
  if (getOption("brms.backend", "not_set") == "not_set") {
    brms_backend <- Sys.getenv("BRMS_BACKEND", "not_set")
    if (brms_backend != "not_set") {
      opts$brms.backend <- brms_backend
    }
  }
  if (getOption("cmdstanr_write_stan_file_dir", "not_set") == "not_set") {
    brms_cache_dir <- Sys.getenv("BRMS_CACHE_DIR", "not_set")
    if (brms_cache_dir != "not_set") {
      opts$cmdstanr_write_stan_file_dir <- brms_cache_dir
    }
  }
  opts
}

brms_prior_model <- function(mix, brms_args, empty = FALSE) {
  withr::with_options(brms_test_options(), {
    model_args <- c(
      brms_args,
      list(
        seed = 1423545,
        refresh = 0,
        sample_prior = "only",
        stanvars = mixstanvar(prior = mix)
      )
    )
    if (empty) {
      model_args$empty <- TRUE
    }
    do.call(brms::brm, model_args)
  })
}

brms_beta_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        r | trials(n) ~ 1,
        family = brms::brmsfamily("binomial", link = "identity"),
        center = FALSE
      ),
      data = data.frame(r = 0, n = 0),
      prior = brms::prior(mixbeta(prior_w, prior_a, prior_b), coef = Intercept)
    )
  )
}

brms_beta_trunc_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        r | trials(n) ~ 1,
        family = brms::brmsfamily("binomial", link = "identity"),
        center = FALSE
      ),
      data = data.frame(r = 0, n = 0),
      prior = brms::prior(
        mixbeta(prior_w, prior_a, prior_b),
        class = b,
        lb = 0.1,
        ub = 0.9
      )
    )
  )
}

brms_normal_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        y ~ 1,
        family = brms::brmsfamily("gaussian", link = "identity"),
        center = FALSE
      ),
      data = data.frame(y = 0),
      prior = brms::prior(mixnorm(prior_w, prior_m, prior_s), coef = Intercept) +
        brms::prior(constant(1), class = sigma)
    )
  )
}

brms_normal_trunc_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        y ~ 1,
        family = brms::brmsfamily("gaussian", link = "identity"),
        center = FALSE
      ),
      data = data.frame(y = 0),
      prior = brms::prior(
        mixnorm(prior_w, prior_m, prior_s),
        class = b,
        lb = -5,
        ub = 5
      ) +
        brms::prior(constant(1), class = sigma)
    )
  )
}

brms_gamma_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        y ~ 1,
        family = brms::brmsfamily("gaussian", link = "identity"),
        center = FALSE
      ),
      data = data.frame(y = 1),
      prior = brms::prior(mixgamma(prior_w, prior_a, prior_b), coef = Intercept) +
        brms::prior(constant(1), class = sigma)
    )
  )
}

brms_gamma_trunc_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        y ~ 1,
        family = brms::brmsfamily("gaussian", link = "identity"),
        center = FALSE
      ),
      data = data.frame(y = 1),
      prior = brms::prior(
        mixgamma(prior_w, prior_a, prior_b),
        class = b,
        lb = 0.1,
        ub = 10
      ) +
        brms::prior(constant(1), class = sigma)
    )
  )
}

brms_mvn_4_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        y ~ 1 + l1 + l2 + l3,
        family = brms::brmsfamily("gaussian", link = "identity"),
        center = FALSE
      ),
      data = data.frame(y = 1, l1 = 0, l2 = 0, l3 = 0),
      prior = brms::prior(mixmvnorm(prior_w, prior_m, prior_sigma_L), class = b) +
        brms::prior(constant(1), class = sigma)
    )
  )
}

brms_mvn_2_args <- function() {
  skip_if_not_installed("brms")

  modifyList(
    brms_sampling_args(),
    list(
      formula = brms::bf(
        y ~ 1 + l1,
        family = brms::brmsfamily("gaussian", link = "identity"),
        center = FALSE
      ),
      data = data.frame(y = 1, l1 = 0),
      prior = brms::prior(mixmvnorm(prior_w, prior_m, prior_sigma_L), class = b) +
        brms::prior(constant(1), class = sigma)
    )
  )
}

expect_mixstanvar_declared <- function(mix, brms_args) {
  skip_if_not_installed("brms")

  brms_prior_empty <- brms_prior_model(mix, brms_args, empty = TRUE)
  mix_class <- gsub("Mix$", "", class(mix)[1])
  stan_dist_lpdf <- paste0("mix", mix_class, "_lpdf")
  stan_dist_lcdf <- paste0("mix", mix_class, "_lcdf")
  stan_dist_lccdf <- paste0("mix", mix_class, "_lccdf")
  stan_dist_cdf <- paste0("mix", mix_class, "_cdf")

  stan_code <- brms::stancode(brms_prior_empty)
  expect_true(
    grep(stan_dist_lpdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density pdf in generated brms Stan code."
  )
  expect_true(
    grep(stan_dist_lcdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density cdf in generated brms Stan code."
  )
  expect_true(
    grep(stan_dist_lccdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density ccdf in generated brms Stan code."
  )
  expect_true(
    grep(stan_dist_cdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density natural scale cdf in generated brms Stan code."
  )

  stan_data <- brms::standata(brms_prior_empty)
  for (i in 1:3) {
    param <- paste0("prior_", rownames(mix)[i])
    expect_true(all(unname(mix[i, ]) == unname(stan_data[[param]])))
  }
}

expect_mixstanvar_mvnorm_declared <- function(mix, brms_args) {
  skip_if_not_installed("brms")

  brms_prior_empty <- brms_prior_model(mix, brms_args, empty = TRUE)
  stan_dist <- paste0("mix", gsub("Mix$", "", class(mix)[1]), "_lpdf")
  expect_true(
    grep(stan_dist, brms::stancode(brms_prior_empty)) == 1,
    info = "Looking for declared Stan mixture density in generated brms Stan code."
  )

  stan_data <- brms::standata(brms_prior_empty)
  Nc <- ncol(mix)
  expect_equal(stan_data$prior_Nc, Nc)
  p <- length(summary(mix)$mean)
  expect_equal(stan_data$prior_p, p)
  expect_equal(stan_data$prior_w, array(unname(mix[1, ])))
  expect_equal(unname(stan_data$prior_m), unname(t(mix[2:(p + 1), ])))
  for (i in seq_len(Nc)) {
    S_c <- matrix(stan_data$prior_sigma[i, , , drop = FALSE], p, p)
    expect_equal(sqrt(diag(S_c)), unname(mix[(p + 2):(1 + 2 * p), i]))
    Rho_c <- cov2cor(S_c)
    expect_equal(
      Rho_c[lower.tri(Rho_c)],
      unname(mix[(1 + 2 * p + 1):nrow(mix), i])
    )
  }
}

expect_mixstanvar_sampled_prior <- function(
  mix,
  brms_args,
  eps,
  qtest,
  ptest = seq(0.2, 0.8, by = 0.2)
) {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  capture.output(brms_prior <- brms_prior_model(mix, brms_args))
  samp <- as.numeric(brms::as_draws_matrix(
    brms_prior,
    variable = "b_Intercept"
  )[, 1])
  qtest_samp <- quantile(samp, ptest)
  qref_qmix <- qmix(mix, ptest)
  expect_equal(unname(qtest_samp), unname(qref_qmix), tolerance = eps)
  ptest_samp <- vapply(qtest, function(q) mean(samp < q), c(0.1))
  pref_pmix <- pmix(mix, qtest)
  expect_equal(unname(ptest_samp), unname(pref_pmix), tolerance = eps)
  invisible(TRUE)
}

KLdiv_mvnorm <- function(m_1, sigma_1, m_2, sigma_2) {
  m_delta <- (m_2 - m_1)
  inv_sigma_2 <- solve(sigma_2)
  p <- length(m_1)
  0.5 *
    (t(m_delta) %*%
      inv_sigma_2 %*%
      m_delta +
      sum(diag(inv_sigma_2 %*% sigma_1)) -
      log(det(sigma_1)) +
      log(det(sigma_2)) -
      p)
}

expect_mixstanvar_mvnorm_sampled_prior <- function(mvmix, brms_args, eps) {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  capture.output(brms_prior <- brms_prior_model(mvmix, brms_args))
  samp <- brms::as_draws_matrix(brms_prior, variable = "^b_", regex = TRUE)
  samp_m <- colMeans(samp)
  samp_sigma <- cov(samp)
  mix_m <- summary(mvmix)$mean
  mix_sigma <- summary(mvmix)$cov
  kl <- KLdiv_mvnorm(samp_m, samp_sigma, mix_m, mix_sigma)
  expect_true(abs(kl) < eps)
}
