## Verbose S2Z fit used by the summary/layout tests.
##
## `RBesT.verbose = TRUE` keeps the sampled parameters (`beta_raw`, `tau_raw`,
## `xi_eta`, `beta_param`) in the draws, which is what makes it observable that
## the downstream selectors do not pick `beta_param` up as `beta`.
fixture <- withr::with_options(
  list(
    RBesT.verbose = TRUE,
    RBesT.MC.ncp = 1,
    RBesT.MC.save_warmup = FALSE
  ),
  {
    set.seed(46711)
    suppressMessages(suppressWarnings(
      gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = AS,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 0.5,
        beta.prior = 2,
        warmup = 1000,
        iter = 2000,
        chains = 4,
        thin = 1
      )
    ))
  }
)
