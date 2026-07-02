fixture <- withr::with_seed(
  23434,
  withr::with_options(
    list(
      RBesT.MC.chains = 1,
      RBesT.MC.control = list(adapt_delta = 0.75, stepsize = 1),
      RBesT.MC.ncp = 0,
      RBesT.MC.save_warmup = FALSE
    ),
    suppressMessages(suppressWarnings(
      gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = AS[1, , drop = FALSE],
        family = binomial,
        tau.dist = "Uniform",
        tau.prior = cbind(0, 1000),
        beta.prior = cbind(0, 1E5),
        iter = 1000,
        warmup = 0,
        chains = 1,
        thin = 1,
        init = 10
      )
    ))
  )
)
