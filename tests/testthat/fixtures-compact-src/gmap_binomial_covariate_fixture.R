trans_cov <- transform(
  transplant,
  country = cut(1:11, c(0, 5, 8, Inf), c("CH", "US", "DE"))
)

fixture_call <- quote(
  gMAP(
    cbind(r, n - r) ~ 1 + country | study,
    data = trans_cov,
    tau.dist = "HalfNormal",
    tau.prior = 1,
    beta.prior = rbind(c(0, 2), c(0, 1), c(0, 1)),
    family = binomial,
    warmup = 1000,
    iter = 2000,
    thin = 1,
    chains = 4
  )
)

fixture_seed <- 34246L
compact_seed <- 873461L
compact_nc <- 1L
compact_digits <- 4L
compact_data_objects <- "trans_cov"

fixture <- compact_gmap_fixture_sample(fixture_call, seed = fixture_seed)
