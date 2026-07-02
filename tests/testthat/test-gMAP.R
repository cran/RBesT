## test gMAP results using SBC and cached gMAP fixtures

test_that("gMAP meets SBC requirements wrt to a Chi-Square statistic.", {
  require(dplyr)
  require(tidyr)
  sbc_chisq_test <- RBesT:::calibration_data %>%
    gather(count.mu, count.tau, key = "param", value = "count") %>%
    group_by(data_scenario, family, sd_tau, param) %>%
    do(as.data.frame(chisq.test(.$count)[c("statistic", "p.value")]))
  num_tests <- nrow(sbc_chisq_test)
  num_failed <- sum(sbc_chisq_test$p.value < 0.05)
  pv <- pbinom(num_failed, num_tests, 0.05)
  expect_true(pv > 0.025 & pv < 0.975)
})

test_that("gMAP meets SBC requirements per bin.", {
  require(dplyr)
  require(tidyr)
  B <- RBesT:::calibration_meta$B
  S <- RBesT:::calibration_meta$S
  alpha <- 0.2
  ptrue <- 1 / B
  crit_low <- qbinom(alpha / 2, S, ptrue)
  crit_high <- qbinom(1 - alpha / 2, S, ptrue)
  sbc_binom_test <- RBesT:::calibration_data %>%
    gather(count.mu, count.tau, key = "param", value = "count") %>%
    group_by(data_scenario, family, sd_tau, param) %>%
    summarise(crit = sum(count < crit_low | count > crit_high)) %>%
    mutate(
      pvalue = pbinom(crit, B, alpha),
      extreme = pvalue < 0.025 | pvalue > 0.975
    )
  num_tests <- nrow(sbc_binom_test)
  num_failed <- sum(sbc_binom_test$extreme)
  pv <- pbinom(num_failed, num_tests, 0.05)
  expect_true(pv > 0.025 & pv < 0.975)
})

test_that("SBC data was up to date at package creation.", {
  calibration_datum <- RBesT:::calibration_meta$created
  package_datum <- RBesT:::pkg_create_date
  delta <- difftime(package_datum, calibration_datum, units = "weeks")
  expect_true(delta < 52. / 2.)
})

## match against respective rstanarm model
test_that("gMAP matches RStanArm binomial family", {
  skip("RStanArm has issues loading since 2024-01-02 on CI/CD systems.")
})

## add test case with a single data
test_that("gMAP processes single trial case", {
  map1 <- load_gmap_fixture("gmap_binomial_single_trial", type = "compact")
  expect_true(nrow(fitted(map1)) == 1)
})

test_that("gMAP processes not continuously labeled studies", {
  map1 <- load_gmap_fixture("gmap_binomial_noncontinuous_studies", type = "compact")
  expect_true(nrow(fitted(map1)) == nrow(AS) - 1)
})

test_that("gMAP reports divergences", {
  mcmc_div <- load_gmap_fixture("gmap_binomial_divergent", type = "mcmc")

  div_draws <- posterior::subset_draws(
    mcmc_div$draws_diag,
    variable = "divergent__"
  )
  expect_true(sum(as.array(div_draws)) > 0)
})

test_that("gMAP handles extreme response rates", {
  map1 <- load_gmap_fixture("gmap_binomial_extreme_all_response", type = "compact")
  expect_true(nrow(fitted(map1)) == 4)

  map2 <- load_gmap_fixture("gmap_binomial_extreme_no_response", type = "compact")
  expect_true(nrow(fitted(map2)) == 4)

  map3 <- load_gmap_fixture("gmap_binomial_extreme_two_studies", type = "compact")
  expect_true(nrow(fitted(map3)) == 4)
})

test_that("gMAP handles fixed tau case", {
  map1 <- load_gmap_fixture("gmap_binomial_fixed_tau", type = "compact")
  expect_true(map1$Rhat.max >= 1)
})

test_that("gMAP labels data rows correctly when using covariates", {
  map_covs <- load_gmap_fixture("gmap_binomial_covariate_rows", type = "compact")
  data_covs <- map_covs$data
  expect_true(all(
    rownames(fitted(map_covs)) ==
      paste(data_covs$study, data_covs$stratum, sep = "/")
  ))

  map_tau_strata <- load_gmap_fixture("gmap_binomial_tau_strata_rows", type = "compact")
  expect_true(all(
    rownames(fitted(map_tau_strata)) == as.character(map_tau_strata$data$id)
  ))
})


test_that("plot.gMAP and forest_plot does not use deprecated ggplot2 size aesthetic", {
  map1 <- load_gmap_fixture("gmap_binomial_fixed_tau", type = "compact")

  withr::with_options(
    list(lifecycle_verbosity = "error", RBesT.verbose = TRUE),
    {
      expect_no_error(plot(map1))
      expect_no_error(forest_plot(map1))
    }
  )
})
