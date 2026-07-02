trans_cov <- 
structure(list(study = c("Study 1", "Study 2", "Study 3", "Study 4", 
"Study 5", "Study 6", "Study 7", "Study 8", "Study 9", "Study 10", 
"Study 11"), n = c(33, 45, 74, 103, 140, 49, 83, 59, 22, 109, 
213), r = c(6, 8, 17, 28, 26, 8, 22, 8, 6, 16, 53), country = structure(c(1L, 
1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L), levels = c("CH", "US", 
"DE"), class = "factor")), class = "data.frame", row.names = c(NA, 
-11L))

gmap_binomial_covariate_compact_spec <- 
structure(list(name = "gmap_binomial_covariate", seed = 873461L, nc = 1L, 
    has_intercept = TRUE, tau_fixed = FALSE, tau_values = NULL, 
    variables_model = c("beta[1]", "beta[2]", "beta[3]", "tau[1]", 
    "lp__"), variables_theta = c("theta[1]", "theta[2]", "theta[3]", 
    "theta[4]", "theta[5]", "theta[6]", "theta[7]", "theta[8]", 
    "theta[9]", "theta[10]", "theta[11]", "theta_pred"), variables = c("theta[1]", 
    "theta[2]", "theta[3]", "theta[4]", "theta[5]", "theta[6]", 
    "theta[7]", "theta[8]", "theta[9]", "theta[10]", "theta[11]", 
    "beta[1]", "beta[2]", "beta[3]", "tau[1]", "theta_pred", 
    "theta_resp_pred", "lp__"), draws_dim = c(1000L, 4L, 18L), 
    diag_variables = c("accept_stat__", "stepsize__", "treedepth__", 
    "n_leapfrog__", "divergent__", "energy__"), metadata_mcmc = list(
        iter = 2000, warmup = 1000, warmup_saved = 0L, chains = 4, 
        n_save_per_chain = 1000, post_warmup_saved = 1000, thin_input = 1, 
        thin_post = 1L, save_warmup = FALSE), generation_config = list(
        recipe = "tests/testthat/fixtures-compact-src/gmap_binomial_covariate_fixture.R", 
        nc = 1L, digits = 4L, seed = 873461L, data_objects = "trans_cov")), class = "compact_gmap_fixture_spec")

gmap_binomial_covariate_compact_spec$mvn_model <- read_mix_json(
  testthat::test_path("fixtures-compact", "gmap_binomial_covariate_mvn_model.json"),
  rescale = TRUE
)
gmap_binomial_covariate_compact_spec$mvn_theta <- read_mix_json(
  testthat::test_path("fixtures-compact", "gmap_binomial_covariate_mvn_theta.json"),
  rescale = TRUE
)
gmap_binomial_covariate_compact_spec$builder <- function() {
  withr::with_options(
    list(RBesT.MC.save_warmup = FALSE),
    suppressWarnings(suppressMessages(eval(
      gMAP(cbind(r, n - r) ~ 1 + country | study, data = trans_cov, 
          tau.dist = "HalfNormal", tau.prior = 1, beta.prior = rbind(c(0, 
              2), c(0, 1), c(0, 1)), family = binomial, warmup = 1000, 
          iter = 2000, thin = 1, chains = 0L)
    )))
  )
}

