test_that("Beta prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$beta
  expect_mixstanvar_declared(mix, brms_beta_args())
})

test_that("Beta mixture prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$betaMix
  expect_mixstanvar_declared(mix, brms_beta_args())
})

test_that("Beta (truncated) prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$beta
  expect_mixstanvar_declared(mix, brms_beta_trunc_args())
})

test_that("Beta mixture (truncated) prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$betaMix
  expect_mixstanvar_declared(mix, brms_beta_trunc_args())
})

test_that("Normal prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$norm
  expect_mixstanvar_declared(mix, brms_normal_args())
})

test_that("Normal mixture prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$normMix
  expect_mixstanvar_declared(mix, brms_normal_args())
})

test_that("Normal (truncated) prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$norm
  expect_mixstanvar_declared(mix, brms_normal_trunc_args())
})

test_that("Normal mixture (truncated) prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$normMix
  expect_mixstanvar_declared(mix, brms_normal_trunc_args())
})

test_that("Gamma prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$gamma
  expect_mixstanvar_declared(mix, brms_gamma_args())
})

test_that("Gamma mixture prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$gammaMix
  expect_mixstanvar_declared(mix, brms_gamma_args())
})

test_that("Gamma (truncated) prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$gamma
  expect_mixstanvar_declared(mix, brms_gamma_trunc_args())
})

test_that("Gamma mixture (truncated) prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_univariate_mixes()$gammaMix
  expect_mixstanvar_declared(mix, brms_gamma_trunc_args())
})

test_that("Multivariate normal (4D) prior is declared correctly in brms generated model and data", {
  mix <- mixstanvar_multivariate_mixes(4)$single
  expect_mixstanvar_mvnorm_declared(mix, brms_mvn_4_args())
})

test_that("Multivariate normal with heavy (4D) tails is declared correctly in brms generated model and data", {
  mix <- mixstanvar_multivariate_mixes(4)$heavy
  expect_mixstanvar_mvnorm_declared(mix, brms_mvn_4_args())
})

test_that("Multivariate normal (2D) is declared correctly in brms generated model and data", {
  mix <- mixstanvar_multivariate_mixes(2)$single
  expect_mixstanvar_mvnorm_declared(mix, brms_mvn_2_args())
})

test_that("Multivariate normal with heavy (2D) is declared correctly in brms generated model and data", {
  mix <- mixstanvar_multivariate_mixes(2)$heavy
  expect_mixstanvar_mvnorm_declared(mix, brms_mvn_2_args())
})
