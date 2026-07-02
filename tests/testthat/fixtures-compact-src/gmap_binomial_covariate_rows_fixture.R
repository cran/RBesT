data_covs <- data.frame(
  n = 10,
  r = 3,
  study = c(1, 2, 2),
  stratum = factor(c("A", "A", "B"))
)
data_covs$group <- paste(data_covs$study, data_covs$stratum, sep = "/")
data_covs$id <- as.integer(factor(data_covs$group))

fixture_call <- quote(
  gMAP(
    cbind(r, n - r) ~ 1 + stratum | study,
    family = binomial,
    data = data_covs,
    tau.dist = "Fixed",
    tau.prior = 0.25,
    beta.prior = 2,
    warmup = 1000,
    iter = 2000,
    chains = 4,
    thin = 1
  )
)

fixture_seed <- 34247L
compact_seed <- 873461L
compact_nc <- 1L
compact_digits <- 4L
compact_data_objects <- "data_covs"

fixture <- suppressWarnings(
  compact_gmap_fixture_sample(fixture_call, seed = fixture_seed)
)
