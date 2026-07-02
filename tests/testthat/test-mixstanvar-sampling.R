test_that("Beta quantiles are correct for brms sampled prior", {
  mix <- mixstanvar_univariate_mixes()$beta
  expect_mixstanvar_sampled_prior(
    mix,
    brms_beta_args(),
    mixstanvar_tolerance(),
    c(0.1, 0.9)
  )
})

test_that("Beta mixture quantiles are correct for brms sampled prior", {
  mix <- mixstanvar_univariate_mixes()$betaMix
  expect_mixstanvar_sampled_prior(
    mix,
    brms_beta_args(),
    mixstanvar_tolerance(),
    c(0.1, 0.9)
  )
})

test_that("Normal quantiles are correct for brms sampled prior", {
  mix <- mixstanvar_univariate_mixes()$norm
  expect_mixstanvar_sampled_prior(
    mix,
    brms_normal_args(),
    mixstanvar_tolerance(),
    c(-1, 0)
  )
})

test_that("Normal mixture quantiles are correct for brms sampled prior", {
  mix <- mixstanvar_univariate_mixes()$normMix
  expect_mixstanvar_sampled_prior(
    mix,
    brms_normal_args(),
    mixstanvar_tolerance(),
    c(-0.5, 1),
    ptest = c(0.4, 0.5, 0.6)
  )
})

test_that("Gamma quantiles are correct for brms sampled prior", {
  mix <- mixstanvar_univariate_mixes()$gamma
  expect_mixstanvar_sampled_prior(
    mix,
    brms_gamma_args(),
    mixstanvar_tolerance(),
    c(2, 7)
  )
})

test_that("Gamma mixture quantile function is correct for brms sampled prior", {
  mix <- mixstanvar_univariate_mixes()$gammaMix
  expect_mixstanvar_sampled_prior(
    mix,
    brms_gamma_args(),
    mixstanvar_tolerance(),
    c(2, 7),
    ptest = seq(0.2, 0.8, by = 0.2)
  )
})

test_that("Multivariate normal (4D) is correct for brms sampled prior", {
  mix <- mixstanvar_multivariate_mixes(4)$single
  expect_mixstanvar_mvnorm_sampled_prior(
    mix,
    brms_mvn_4_args(),
    mixstanvar_tolerance()
  )
})

test_that("Multivariate normal with heavy (4D) tails is correct for brms sampled prior", {
  mix <- mixstanvar_multivariate_mixes(4)$heavy
  expect_mixstanvar_mvnorm_sampled_prior(
    mix,
    brms_mvn_4_args(),
    mixstanvar_tolerance()
  )
})

test_that("Multivariate normal (2D) is correct for brms sampled prior", {
  mix <- mixstanvar_multivariate_mixes(2)$single
  expect_mixstanvar_mvnorm_sampled_prior(
    mix,
    brms_mvn_2_args(),
    mixstanvar_tolerance()
  )
})

test_that("Multivariate normal with heavy (2D) tails is correct for brms sampled prior", {
  mix <- mixstanvar_multivariate_mixes(2)$heavy
  expect_mixstanvar_mvnorm_sampled_prior(
    mix,
    brms_mvn_2_args(),
    mixstanvar_tolerance()
  )
})
