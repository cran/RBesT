gmap_binomial_single_trial_compact_spec <- 
structure(list(name = "gmap_binomial_single_trial", seed = 873461L, nc = 2L, 
    has_intercept = TRUE, tau_fixed = FALSE, tau_values = NULL, 
    variables_model = c("beta[1]", "tau[1]", "lp__"), variables_theta = c("theta[1]", 
    "theta_pred"), variables = c("theta[1]", "beta[1]", "tau[1]", 
    "theta_pred", "theta_resp_pred", "lp__"), draws_dim = c(1000L, 
    4L, 6L), diag_variables = c("accept_stat__", "stepsize__", 
    "treedepth__", "n_leapfrog__", "divergent__", "energy__"), 
    metadata_mcmc = list(iter = 2000, warmup = 1000, warmup_saved = 0L, 
        chains = 4, n_save_per_chain = 1000, post_warmup_saved = 1000, 
        thin_input = 1, thin_post = 1L, save_warmup = FALSE), 
    generation_config = list(recipe = "tests/testthat/fixtures-compact-src/gmap_binomial_single_trial_fixture.R", 
        nc = 2L, digits = 4L, seed = 873461L, data_objects = character(0))), class = "compact_gmap_fixture_spec")

gmap_binomial_single_trial_compact_spec$mvn_model <- read_mix_json(
  testthat::test_path("fixtures-compact", "gmap_binomial_single_trial_mvn_model.json"),
  rescale = TRUE
)
gmap_binomial_single_trial_compact_spec$mvn_theta <- read_mix_json(
  testthat::test_path("fixtures-compact", "gmap_binomial_single_trial_mvn_theta.json"),
  rescale = TRUE
)
gmap_binomial_single_trial_compact_spec$builder <- function() {
  withr::with_options(
    list(RBesT.MC.save_warmup = FALSE),
    suppressWarnings(suppressMessages(eval(
      gMAP(cbind(r, n - r) ~ 1, data = colitis[1, ], family = binomial, 
          tau.dist = "HalfNormal", tau.prior = c(0.5), beta.prior = cbind(0, 
              1), warmup = 1000, iter = 2000, chains = 0L, thin = 1)
    )))
  )
}

