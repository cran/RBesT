fixture_call <- quote(
  gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = AS[-1, ],
    family = binomial,
    tau.dist = "HalfNormal",
    tau.prior = 0.5,
    warmup = 1000,
    iter = 2000,
    chains = 4,
    thin = 1
  )
)

fixture_seed <- 34251L
compact_seed <- 873461L
compact_nc <- 1L
compact_digits <- 4L

fixture <- suppressWarnings(
  compact_gmap_fixture_sample(fixture_call, seed = fixture_seed)
)
