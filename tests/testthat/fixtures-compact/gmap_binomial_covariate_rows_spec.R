data_covs <- 
structure(list(n = c(10, 10, 10), r = c(3, 3, 3), study = c(1, 
2, 2), stratum = structure(c(1L, 1L, 2L), levels = c("A", "B"
), class = "factor"), group = c("1/A", "2/A", "2/B"), id = 1:3), row.names = c(NA, 
-3L), class = "data.frame")

gmap_binomial_covariate_rows_compact_spec <- 
structure(list(name = "gmap_binomial_covariate_rows", seed = 873461L, nc = 1L, 
    has_intercept = TRUE, tau_fixed = TRUE, tau_values = c(`tau[1]` = 0.25), 
    variables_model = c("beta[1]", "beta[2]", "lp__"), variables_theta = c("theta[1]", 
    "theta[2]", "theta[3]", "theta_pred"), variables = c("theta[1]", 
    "theta[2]", "theta[3]", "beta[1]", "beta[2]", "tau[1]", "theta_pred", 
    "theta_resp_pred", "lp__"), draws_dim = c(1000L, 4L, 9L), 
    diag_variables = c("accept_stat__", "stepsize__", "treedepth__", 
    "n_leapfrog__", "divergent__", "energy__"), metadata_mcmc = list(
        iter = 2000, warmup = 1000, warmup_saved = 0L, chains = 4, 
        n_save_per_chain = 1000, post_warmup_saved = 1000, thin_input = 1, 
        thin_post = 1L, save_warmup = FALSE), generation_config = list(
        recipe = "tests/testthat/fixtures-compact-src/gmap_binomial_covariate_rows_fixture.R", 
        nc = 1L, digits = 4L, seed = 873461L, data_objects = "data_covs")), class = "compact_gmap_fixture_spec")

gmap_binomial_covariate_rows_compact_spec$mvn_model <- read_mix_json(
  testthat::test_path("fixtures-compact", "gmap_binomial_covariate_rows_mvn_model.json"),
  rescale = TRUE
)
gmap_binomial_covariate_rows_compact_spec$mvn_theta <- read_mix_json(
  testthat::test_path("fixtures-compact", "gmap_binomial_covariate_rows_mvn_theta.json"),
  rescale = TRUE
)
gmap_binomial_covariate_rows_compact_spec$builder <- function() {
  withr::with_options(
    list(RBesT.MC.save_warmup = FALSE),
    suppressWarnings(suppressMessages(eval(
      gMAP(cbind(r, n - r) ~ 1 + stratum | study, family = binomial, 
          data = data_covs, tau.dist = "Fixed", tau.prior = 0.25, beta.prior = 2, 
          warmup = 1000, iter = 2000, chains = 0L, thin = 1)
    )))
  )
}

