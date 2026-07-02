fixture <- withr::with_seed(
  34563,
  withr::with_options(
    list(RBesT.MC.save_warmup = TRUE),
    gMAP(
      cbind(r, n - r) ~ 1 | study,
      family = binomial,
      data = AS,
      tau.dist = "Fixed",
      tau.prior = 0.5,
      beta.prior = 2,
      warmup = 100,
      iter = 200,
      chains = 2,
      thin = 1
    )
  )
)
