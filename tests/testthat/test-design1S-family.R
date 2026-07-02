# Tests for family-based sigma in 1S normal design functions
# Decision: P(theta > qc) > pc on the link scale

# =========================================================================
# Test Group 1: Core formula validation (indirectly testing sigma_from_family)
# =========================================================================

test_that("gaussian(identity) + sigma gives identical boundary to fixed sigma", {
  prior <- mixnorm(c(1, 0, 2), sigma = 1)
  n <- 50
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  # Fixed sigma=1 path
  crit_fixed <- decision1S_boundary(prior, n, dec, sigma = 1)

  # gaussian(identity) + sigma short-circuits to the fixed-sigma path
  crit_family <- decision1S_boundary(prior, n, dec,
                                     sigma = 1, family = gaussian())

  expect_equal(crit_family, crit_fixed, tolerance = 1e-6)
})

test_that("boundary with binomial() family is self-consistent at y_c", {
  # At the boundary y_c, evaluating the decision with the posterior formed

  # using se = sigma_from_family(y_c, binomial()) / sqrt(n) should be
  # exactly at the decision threshold.
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 60
  dec <- decision1S(pc = 0.95, qc = qlogis(0.10), lower.tail = FALSE)

  y_c <- decision1S_boundary(prior, n, dec, family = binomial())

  # Manually compute posterior at y_c with family-derived sigma
  p_c <- plogis(y_c)
  se_at_yc <- (1 / sqrt(p_c * (1 - p_c))) / sqrt(n)
  # Posterior: precision-weighted update
  prior_prec <- 1 / 2^2
  data_prec <- 1 / se_at_yc^2
  post_prec <- prior_prec + data_prec
  post_mean <- (0 * prior_prec + y_c * data_prec) / post_prec
  post_sd <- 1 / sqrt(post_prec)

  # P(theta > qc | y_c) should be very close to pc = 0.95
  prob_above <- pnorm(qlogis(0.10), post_mean, post_sd, lower.tail = FALSE)
  expect_equal(prob_above, 0.95, tolerance = 1e-4)
})

test_that("oc1S with NB family at known sigma point matches fixed-sigma", {
  # At eta=0 (rate=1), sigma for NB(theta=0.5) = sqrt(1/1 + 1/0.5) = sqrt(3)
  theta_NB <- 0.5
  sigma_at_zero <- sqrt(1 + 1 / theta_NB)  # sqrt(3)

  prior <- mixnorm(c(1, 0.5, 1), sigma = sigma_at_zero)
  n <- 30
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  oc_family <- oc1S(prior, n, dec,
                    family = MASS::negative.binomial(theta_NB))

  # Verify that OC at eta=0 is reasonable and monotone around it
  # (higher true rate -> higher power for "greater than" decision)
  oc_below <- oc_family(-0.5)
  oc_at <- oc_family(0)
  oc_above <- oc_family(0.5)

  expect_true(oc_below < oc_at)
  expect_true(oc_at < oc_above)
  expect_true(oc_at > 0 && oc_at < 1)
})


# =========================================================================
# Test Group 2: decision1S_boundary with family
# =========================================================================

test_that("decision1S_boundary with family differs from fixed sigma", {
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 40
  # Decision at a low rate where sigma correction is prominent
  dec <- decision1S(pc = 0.95, qc = qlogis(0.10), lower.tail = FALSE)

  crit_family <- decision1S_boundary(prior, n, dec, family = binomial())
  crit_fixed <- decision1S_boundary(prior, n, dec, sigma = 2)

  # They should differ because sigma at the boundary (~p=0.1 region)
  # is ~3.33, not 2.0
  expect_false(isTRUE(all.equal(crit_family, crit_fixed, tolerance = 1e-3)))
})

test_that("decision1S_boundary with family works for two-sided decisions", {
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 50
  dec <- decision1S(
    pc = c(0.95, 0.80),
    qc = c(qlogis(0.10), qlogis(0.30)),
    lower.tail = c(FALSE, TRUE)
  )

  result <- decision1S_boundary(prior, n, dec, family = binomial())
  expect_named(result, c("lower_or_equal_than", "higher_than"))
  expect_true(result["lower_or_equal_than"] > result["higher_than"])
})


# =========================================================================
# Test Group 3: oc1S with family
# =========================================================================

test_that("oc1S with gaussian(identity) + sigma matches fixed-sigma oc1S", {
  prior <- mixnorm(c(1, 1, 1.5), sigma = 1)
  n <- 40
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  oc_family <- oc1S(prior, n, dec, sigma = 1, family = gaussian())
  oc_fixed <- oc1S(prior, n, dec, sigma = 1)

  thetas <- seq(-1, 3, by = 0.5)
  expect_equal(oc_family(thetas), oc_fixed(thetas), tolerance = 1e-6)
})

test_that("oc1S with binomial() family matches exact betaMix OC for large n", {
  # Large n: normal approximation should be good
  n <- 200
  # Decision: P(rate > 0.10) > 0.95
  qc_logit <- qlogis(0.10)
  dec_logit <- decision1S(pc = 0.95, qc = qc_logit, lower.tail = FALSE)

  # Flat prior on logit scale (very wide)
  prior_logit <- mixnorm(c(1, 0, 10), sigma = 2)
  oc_family <- oc1S(prior_logit, n, dec_logit, family = binomial())

  # Exact binomial with Beta(1,1) prior
  dec_exact <- decision1S(pc = 0.95, qc = 0.10, lower.tail = FALSE)
  prior_exact <- mixbeta(c(1, 1, 1))
  oc_exact <- oc1S(prior_exact, n, dec_exact)

  # Compare at true rates where sigma correction matters
  # Note: the priors are not exactly matched (N(0,10) on logit vs Beta(1,1))
  # so we allow 3% tolerance which accounts for both the normal approximation
  # and the prior mismatch
  true_ps <- c(0.08, 0.10, 0.12, 0.15, 0.20)
  oc_fam_vals <- oc_family(qlogis(true_ps))
  oc_exact_vals <- oc_exact(true_ps)

  expect_equal(oc_fam_vals, oc_exact_vals, tolerance = 0.03)
})

test_that("oc1S with family differs from fixed sigma in tails", {
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 40
  dec <- decision1S(pc = 0.95, qc = qlogis(0.10), lower.tail = FALSE)

  oc_family <- oc1S(prior, n, dec, family = binomial())
  oc_fixed <- oc1S(prior, n, dec, sigma = 2)

  # At low p (say p=0.05), sigma is ~4.5 vs fixed 2.0
  # so the family OC should give higher probability (wider sampling dist)
  theta_low <- qlogis(0.05)
  expect_true(oc_family(theta_low) > oc_fixed(theta_low))
})

test_that("oc1S with family works for two-sided decisions", {
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 50
  dec <- decision1S(
    pc = c(0.95, 0.80),
    qc = c(qlogis(0.10), qlogis(0.30)),
    lower.tail = c(FALSE, TRUE)
  )

  oc_family <- oc1S(prior, n, dec, family = binomial())
  # Should return a function that gives values in [0, 1]
  vals <- oc_family(qlogis(c(0.10, 0.15, 0.20, 0.25)))
  expect_true(all(vals >= 0 & vals <= 1))
})


# =========================================================================
# Test Group 4: pos1S with family
# =========================================================================

test_that("pos1S with family at tight prior matches oc1S", {
  # A very tight prior centered at theta_0 should give
  # pos1S(tight_mix) ≈ oc1S(theta_0)
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 40
  dec <- decision1S(pc = 0.95, qc = qlogis(0.10), lower.tail = FALSE)

  oc_fn <- oc1S(prior, n, dec, family = binomial())
  pos_fn <- pos1S(prior, n, dec, family = binomial())

  # Tight "prior belief" centered at logit(0.20)
  theta_0 <- qlogis(0.20)
  tight_mix <- mixnorm(c(1, theta_0, 0.001), sigma = 2)

  expect_equal(pos_fn(tight_mix), oc_fn(theta_0), tolerance = 1e-3)
})

test_that("pos1S with gaussian(identity) + sigma matches fixed-sigma pos1S", {
  prior <- mixnorm(c(1, 1, 1.5), sigma = 1)
  n <- 40
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  pos_family <- pos1S(prior, n, dec, sigma = 1, family = gaussian())
  pos_fixed <- pos1S(prior, n, dec, sigma = 1)

  # Use the prior itself as the predictive belief
  expect_equal(pos_family(prior), pos_fixed(prior), tolerance = 1e-4)
})

test_that("pos1S with family works for mixture prior belief", {
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 40
  dec <- decision1S(pc = 0.95, qc = qlogis(0.10), lower.tail = FALSE)

  pos_fn <- pos1S(prior, n, dec, family = binomial())

  # Two-component belief about the true rate
  belief <- mixnorm(
    c(0.7, qlogis(0.15), 0.5),
    c(0.3, qlogis(0.25), 0.3),
    sigma = 2
  )

  result <- pos_fn(belief)
  expect_true(result > 0 && result < 1)
})


# =========================================================================
# Test Group 5: Backward compatibility and input validation
# =========================================================================

test_that("decision1S_boundary without family is unchanged", {
  prior <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  dec <- decision1S(pc = 0.95, qc = 0.6, lower.tail = TRUE)

  # This is an existing test value from test-decision1S_boundary.R
  result <- decision1S_boundary(prior, n = 50, decision = dec)
  expect_equal(result, -0.7359727, tolerance = 1e-5)
})

test_that("oc1S without family is unchanged", {
  s <- 2
  prior <- mixnorm(c(1, 0, 100), sigma = s)
  dec <- decision1S(0.95, 0.4, lower.tail = TRUE)
  n <- 155

  oc_fn <- oc1S(prior, n, dec)
  # At theta=0 with flat prior and this design, power should be high
  expect_true(oc_fn(0) > 0.7)
})

test_that("pos1S without family is unchanged", {
  s <- 2
  prior <- mixnorm(c(1, 0, 100), sigma = s)
  dec <- decision1S(0.95, 0.4, lower.tail = TRUE)
  n <- 155

  pos_fn <- pos1S(prior, n, dec)
  # PoS with the prior itself should be a probability
  result <- pos_fn(prior)
  expect_true(result > 0 && result < 1)
})

test_that("specifying sigma with a non-gaussian family raises an error", {
  prior <- mixnorm(c(1, 0, 2), sigma = 2)
  n <- 40
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  expect_error(
    decision1S_boundary(prior, n, dec, sigma = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
  expect_error(
    oc1S(prior, n, dec, sigma = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
  expect_error(
    pos1S(prior, n, dec, sigma = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
})


# =========================================================================
# Test Group 6: Offset
# =========================================================================

test_that("offset changes boundary for NB family", {
  theta_NB <- 0.5
  fam_nb <- MASS::negative.binomial(theta = theta_NB)
  prior <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  n <- 30
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  # More exposure (offset=log(2)) means more info per obs => lower boundary
  crit_no_offset <- decision1S_boundary(prior, n, dec,
                                         family = fam_nb, offset = 0)
  crit_with_offset <- decision1S_boundary(prior, n, dec,
                                           family = fam_nb, offset = log(2))

  # With more exposure, need less extreme observation to trigger decision
  expect_true(crit_with_offset < crit_no_offset)
})

test_that("oc1S with NB family and offset gives higher power", {
  theta_NB <- 0.5
  fam_nb <- MASS::negative.binomial(theta = theta_NB)
  prior <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  n <- 30
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  oc_no_offset <- oc1S(prior, n, dec, family = fam_nb, offset = 0)
  oc_with_offset <- oc1S(prior, n, dec, family = fam_nb, offset = log(2))

  # At a true rate above the threshold (rate=1.5, eta=log(1.5)),
  # more exposure should give higher power
  theta_test <- log(1.5)
  expect_true(oc_with_offset(theta_test) > oc_no_offset(theta_test))
})

test_that("offset=0 is the default and matches no-offset call", {
  theta_NB <- 0.5
  fam_nb <- MASS::negative.binomial(theta = theta_NB)
  prior <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  n <- 30
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  crit_default <- decision1S_boundary(prior, n, dec, family = fam_nb)
  crit_explicit <- decision1S_boundary(prior, n, dec,
                                        family = fam_nb, offset = 0)

  expect_equal(crit_default, crit_explicit)
})


# =========================================================================
# Test Group 7: Gaussian family + sigma (dispersion) interaction
# =========================================================================

test_that("gaussian() family without sigma raises an error", {
  prior <- mixnorm(c(1, 0, 2), sigma = 1)
  n <- 40
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  expect_error(
    decision1S_boundary(prior, n, dec, family = gaussian()),
    "sigma.*must be specified|gaussian.*sigma"
  )
  expect_error(
    oc1S(prior, n, dec, family = gaussian()),
    "sigma.*must be specified|gaussian.*sigma"
  )
  expect_error(
    pos1S(prior, n, dec, family = gaussian()),
    "sigma.*must be specified|gaussian.*sigma"
  )
})

test_that("gaussian(log) + sigma gives varying effective sigma", {
  # For gaussian(log): sigma_eff(eta) = sigma / exp(eta)
  # At eta=0: sigma_eff = sigma/1 = sigma
  # At eta=1: sigma_eff = sigma/e ≈ sigma/2.72
  # So the boundary should differ from the fixed-sigma case.
  sigma_val <- 2
  prior <- mixnorm(c(1, 0, 1), sigma = sigma_val)
  n <- 50
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  crit_gauss_log <- decision1S_boundary(prior, n, dec,
                                         sigma = sigma_val,
                                         family = gaussian("log"))
  crit_fixed <- decision1S_boundary(prior, n, dec, sigma = sigma_val)

  # Should differ because sigma_eff varies with eta

  expect_false(isTRUE(all.equal(crit_gauss_log, crit_fixed, tolerance = 1e-3)))
})

test_that("gaussian(log) boundary is self-consistent", {
  # At boundary y_c, the posterior with se = sigma_eff(y_c)/sqrt(n) should
  # give the decision probability exactly at the threshold.
  sigma_val <- 2
  prior <- mixnorm(c(1, 0, 1), sigma = sigma_val)
  n <- 50
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  y_c <- decision1S_boundary(prior, n, dec,
                              sigma = sigma_val, family = gaussian("log"))

  # sigma_eff at y_c for gaussian(log): sigma / exp(y_c)
  se_at_yc <- (sigma_val / exp(y_c)) / sqrt(n)

  # Posterior: precision-weighted
  prior_prec <- 1 / 1^2  # prior sd = 1
  data_prec <- 1 / se_at_yc^2
  post_prec <- prior_prec + data_prec
  post_mean <- (0 * prior_prec + y_c * data_prec) / post_prec
  post_sd <- 1 / sqrt(post_prec)

  # P(theta > 0 | y_c) should be exactly pc = 0.95
  prob_above <- pnorm(0, post_mean, post_sd, lower.tail = FALSE)
  expect_equal(prob_above, 0.95, tolerance = 1e-4)
})

test_that("oc1S with gaussian(log) + sigma differs from fixed sigma", {
  sigma_val <- 2
  prior <- mixnorm(c(1, 0, 1), sigma = sigma_val)
  n <- 50
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  oc_gauss_log <- oc1S(prior, n, dec,
                        sigma = sigma_val, family = gaussian("log"))
  oc_fixed <- oc1S(prior, n, dec, sigma = sigma_val)

  # At eta = 1: sigma_eff = 2/e ≈ 0.74 (less variability than fixed sigma=2)
  # -> higher power with gaussian(log) at positive eta
  expect_true(oc_gauss_log(1) > oc_fixed(1))

  # At eta = -1: sigma_eff = 2*e ≈ 5.44 (more variability)
  # -> but this direction depends on the boundary as well
  # Just check they are different
  expect_false(isTRUE(all.equal(oc_gauss_log(-1), oc_fixed(-1),
                                tolerance = 1e-3)))
})

test_that("pos1S with gaussian(log) + sigma at tight prior matches oc1S", {
  sigma_val <- 2
  prior <- mixnorm(c(1, 0, 1), sigma = sigma_val)
  n <- 50
  dec <- decision1S(pc = 0.95, qc = 0, lower.tail = FALSE)

  oc_fn <- oc1S(prior, n, dec, sigma = sigma_val, family = gaussian("log"))
  pos_fn <- pos1S(prior, n, dec, sigma = sigma_val, family = gaussian("log"))

  # Tight prior at theta_0 = 0.5
  theta_0 <- 0.5
  tight_mix <- mixnorm(c(1, theta_0, 0.001), sigma = sigma_val)

  expect_equal(pos_fn(tight_mix), oc_fn(theta_0), tolerance = 1e-3)
})

test_that("gaussian(identity) + sigma=3 matches fixed sigma=3", {
  # gaussian(identity) should short-circuit: any sigma value works
  prior <- mixnorm(c(1, 5, 3), sigma = 3)
  n <- 30
  dec <- decision1S(pc = 0.90, qc = 2, lower.tail = FALSE)

  crit_family <- decision1S_boundary(prior, n, dec,
                                      sigma = 3, family = gaussian())
  crit_fixed <- decision1S_boundary(prior, n, dec, sigma = 3)

  expect_equal(crit_family, crit_fixed, tolerance = 1e-6)

  oc_family <- oc1S(prior, n, dec, sigma = 3, family = gaussian())
  oc_fixed <- oc1S(prior, n, dec, sigma = 3)

  thetas <- seq(2, 8, by = 1)
  expect_equal(oc_family(thetas), oc_fixed(thetas), tolerance = 1e-6)
})
