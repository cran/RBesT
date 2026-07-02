fixture_call <- quote(
  gMAP(
    cbind(r, n - r) ~ 1 | study,
    family = binomial,
    data = AS,
    tau.dist = "Fixed",
    tau.prior = 0.5,
    beta.prior = 2,
    warmup = 1000,
    iter = 2000,
    chains = 4,
    thin = 1
  )
)

fixture_seed <- 34563L
compact_seed <- 873461L
compact_nc <- 1L
compact_digits <- 4L

fixture <- compact_gmap_fixture_sample(fixture_call, seed = fixture_seed)
