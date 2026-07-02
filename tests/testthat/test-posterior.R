# test S3 methods in alphabetical order
test_that("as_draws and friends have resonable outputs", {
  map <- load_gmap_fixture("gmap_binomial_fixed_tau", type = "compact")

  n_iter <- map$metadata_mcmc$iter
  n_warmup <- map$metadata_mcmc$warmup
  n_chains <- map$metadata_mcmc$chains

  draws <- as_draws(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- suppressMessages(as_draws_matrix(
    map,
    variable = "theta_resp_pred"
  ))
  expect_s3_class(draws, "draws_matrix")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_array(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_array")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_df(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_df")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_list(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_rvars(map)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::nvariables(draws) > 0)
  expect_equal(posterior::ndraws(draws), nsamples(map))

  expect_equal(
    posterior::niterations(map$draws),
    n_iter - n_warmup
  )
  expect_equal(
    posterior::nchains(map$draws),
    n_chains
  )
  expect_equal(posterior::ndraws(draws), (n_iter - n_warmup) * n_chains)
  expect_equal(nsamples(map), (n_iter - n_warmup) * n_chains)
})

test_that("as_draws_rvars can include warmup draws", {
  map_full <- load_gmap_fixture("gmap_binomial_fixed_tau_warmup", type = "mcmc")

  n_iter_full <- map_full$metadata_mcmc$iter
  n_warmup_full <- map_full$metadata_mcmc$warmup
  n_chains_full <- map_full$metadata_mcmc$chains

  draws_full <- as_draws_rvars(map_full, inc_warmup = TRUE)
  expect_s3_class(draws_full, "draws_rvars")
  expect_true(posterior::nvariables(draws_full) > 0)
  expect_equal(
    posterior::ndraws(draws_full),
    nsamples(map_full) + n_warmup_full * n_chains_full
  )

  expect_equal(
    posterior::niterations(map_full$draws),
    n_iter_full - n_warmup_full
  )
  expect_equal(
    posterior::niterations(map_full$draws_warmup),
    n_warmup_full
  )
  expect_equal(
    posterior::nchains(map_full$draws),
    n_chains_full
  )
  expect_equal(posterior::ndraws(draws_full), n_iter_full * n_chains_full)
  expect_equal(nsamples(map_full), (n_iter_full - n_warmup_full) * n_chains_full)
})

test_that("as_draws methods warn for draw-free gMAP skeletons", {
  map <- suppressMessages(gMAP(
    cbind(r, n - r) ~ 1 | study,
    family = binomial,
    data = AS,
    tau.dist = "Fixed",
    tau.prior = 0.5,
    beta.prior = 2,
    chains = 0
  ))

  expect_warning(
    as_draws(map),
    "gMAP object \"map\" does not contain any samples"
  )
  expect_warning(
    posterior::as_draws_matrix(map),
    "gMAP object \"map\" does not contain any samples"
  )
  expect_warning(
    posterior::as_draws_array(map),
    "gMAP object \"map\" does not contain any samples"
  )
  expect_warning(
    posterior::as_draws_df(map),
    "gMAP object \"map\" does not contain any samples"
  )
  expect_warning(
    posterior::as_draws_list(map),
    "gMAP object \"map\" does not contain any samples"
  )
  expect_warning(
    posterior::as_draws_rvars(map),
    "gMAP object \"map\" does not contain any samples"
  )
})
