data_extreme <- data.frame(
  n = c(5, 5, 5, 5),
  r = c(5, 5, 5, 5),
  study = c(1, 1, 2, 2)
)

fixture_call <- quote(
  gMAP(
    cbind(r, n - r) ~ 1 | study,
    family = binomial,
    data = data_extreme,
    tau.dist = "HalfNormal",
    tau.prior = 2.0,
    beta.prior = 2,
    warmup = 1000,
    iter = 2000,
    chains = 4,
    thin = 1
  )
)

fixture_seed <- 34250L
compact_seed <- 873461L
compact_nc <- 1L
compact_digits <- 4L
compact_data_objects <- "data_extreme"

fixture <- suppressWarnings(
  compact_gmap_fixture_sample(fixture_call, seed = fixture_seed)
)
