# Tests for family-based sigma in 2S normal design functions

# =========================================================================
# Test Group 1: Backward compatibility
# =========================================================================

test_that("decision2S_boundary without family is unchanged", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  crit_fn <- decision2S_boundary(prior1, prior2, 30, 30, dec, sigma1 = 2, sigma2 = 2)
  # Should return a function
  expect_true(is.function(crit_fn))
  # Boundary over a range should give finite values
  y2_grid <- seq(-1, 1, by = 0.5)
  vals <- crit_fn(y2_grid)
  expect_true(all(is.finite(vals)))
})

test_that("oc2S without family is unchanged", {
  prior <- mixnorm(c(1, 0, 100), sigma = 2)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  oc_fn <- oc2S(prior, prior, n1, n2, dec)
  # At equal thetas (type I error), should be < 0.05
  expect_true(oc_fn(0, 0) < 0.05)
  # With a true difference, power should be moderate
  expect_true(oc_fn(1, 0) > 0.1)
})

test_that("pos2S without family is unchanged", {
  prior <- mixnorm(c(1, 0, 100), sigma = 2)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  pos_fn <- pos2S(prior, prior, n1, n2, dec)
  # PoS with the priors themselves
  result <- pos_fn(prior, prior)
  expect_true(result > 0 && result < 1)
})


# =========================================================================
# Test Group 2: Gaussian(identity) + sigma short-circuits
# =========================================================================

test_that("decision2S_boundary with gaussian(identity) matches fixed sigma", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  crit_fixed <- decision2S_boundary(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)
  crit_family <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                      sigma1 = 2, sigma2 = 2,
                                      family = gaussian())

  # Evaluate at a few y2 points
  y2_pts <- seq(-2, 2, by = 1)
  expect_equal(crit_fixed(y2_pts), crit_family(y2_pts), tolerance = 1e-4)
})

test_that("oc2S with gaussian(identity) matches fixed sigma", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  oc_fixed <- oc2S(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)
  oc_family <- oc2S(prior1, prior2, n1, n2, dec,
                     sigma1 = 2, sigma2 = 2, family = gaussian())

  expect_equal(oc_fixed(1, 0), oc_family(1, 0), tolerance = 1e-4)
  expect_equal(oc_fixed(0, 0), oc_family(0, 0), tolerance = 1e-4)
})


# =========================================================================
# Test Group 3: Binomial family boundary
# =========================================================================

test_that("decision2S_boundary with binomial() family differs from fixed sigma", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, qlogis(0.10), lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  y2_pts <- qlogis(c(0.05, 0.10, 0.15, 0.20))
  crit_family <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                      family = binomial())
  crit_fixed <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                     sigma1 = 2, sigma2 = 2)

  # Boundaries should differ at some points
  expect_false(isTRUE(all.equal(crit_family(y2_pts), crit_fixed(y2_pts),
                                tolerance = 1e-2)))
})

test_that("decision2S_boundary with binomial() returns a function", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, qlogis(0.10), lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  crit_fn <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                  family = binomial())
  expect_true(is.function(crit_fn))

  # Should give finite values across a range of y2
  y2_pts <- qlogis(c(0.05, 0.10, 0.15, 0.20))
  vals <- crit_fn(y2_pts)
  expect_true(all(is.finite(vals)))
})


# =========================================================================
# Test Group 4: OC with family
# =========================================================================

test_that("oc2S with binomial() family gives reasonable type I error", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 60
  n2 <- 60

  oc_fn <- oc2S(prior1, prior2, n1, n2, dec, family = binomial())

  # Under null (both arms at same rate, logit(0.30)):
  theta_null <- qlogis(0.30)
  type_I <- oc_fn(theta_null, theta_null)
  expect_true(type_I < 0.10)

  # Under alternative (arm 1 higher than arm 2):
  theta_alt1 <- qlogis(0.50)
  theta_alt2 <- qlogis(0.20)
  power <- oc_fn(theta_alt1, theta_alt2)
  expect_true(power > type_I)
})

test_that("oc2S with binomial() family differs from fixed sigma in tails", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  oc_family <- oc2S(prior1, prior2, n1, n2, dec, family = binomial())
  oc_fixed <- oc2S(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)

  # At extreme rates where sigma differs most from 2
  theta1 <- qlogis(0.05)
  theta2 <- qlogis(0.05)
  expect_false(isTRUE(all.equal(oc_family(theta1, theta2),
                                oc_fixed(theta1, theta2),
                                tolerance = 1e-2)))
})


# =========================================================================
# Test Group 5: Offset
# =========================================================================

test_that("offset changes boundary for NB family in 2S case", {
  theta_NB <- 0.5
  fam_nb <- MASS::negative.binomial(theta = theta_NB)
  prior1 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  prior2 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  n1 <- 30
  n2 <- 30
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  crit_no_offset <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                         family = fam_nb)
  crit_with_offset <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                           family = fam_nb,
                                           offset1 = log(2), offset2 = log(2))

  # Boundaries should differ with different exposure
  y2_pts <- seq(0, 1, by = 0.25)
  expect_false(isTRUE(all.equal(crit_no_offset(y2_pts),
                                crit_with_offset(y2_pts),
                                tolerance = 1e-2)))
})

test_that("offset2 defaults to offset1", {
  theta_NB <- 0.5
  fam_nb <- MASS::negative.binomial(theta = theta_NB)
  prior1 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  prior2 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  n1 <- 30
  n2 <- 30
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  # offset1 = log(2), offset2 should default to log(2)
  crit_default <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                       family = fam_nb, offset1 = log(2))
  crit_explicit <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                        family = fam_nb,
                                        offset1 = log(2), offset2 = log(2))

  y2_pts <- seq(0, 1, by = 0.25)
  expect_equal(crit_default(y2_pts), crit_explicit(y2_pts), tolerance = 1e-6)
})


# =========================================================================
# Test Group 6: Error conditions
# =========================================================================

test_that("specifying sigma with non-gaussian family raises error in 2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  expect_error(
    decision2S_boundary(prior1, prior2, 30, 30, dec,
                        sigma1 = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
  expect_error(
    oc2S(prior1, prior2, 30, 30, dec,
         sigma1 = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
  expect_error(
    pos2S(prior1, prior2, 30, 30, dec,
          sigma1 = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
})

test_that("gaussian() family without sigma raises error in 2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  expect_error(
    decision2S_boundary(prior1, prior2, 30, 30, dec, family = gaussian()),
    "sigma.*must be specified|gaussian.*sigma"
  )
})

test_that("invalid family object raises error in 2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  expect_error(
    decision2S_boundary(prior1, prior2, 30, 30, dec, family = "binomial"),
    "family.*must be"
  )
})


# =========================================================================
# Test Group 7: PoS consistency
# =========================================================================

test_that("pos2S with family at tight priors matches oc2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, qlogis(0.10), lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  oc_fn <- oc2S(prior1, prior2, n1, n2, dec, family = binomial())
  pos_fn <- pos2S(prior1, prior2, n1, n2, dec, family = binomial())

  # Tight priors centered at specific values
  theta1 <- qlogis(0.20)
  theta2 <- qlogis(0.10)
  tight1 <- mixnorm(c(1, theta1, 0.001), sigma = 2)
  tight2 <- mixnorm(c(1, theta2, 0.001), sigma = 2)

  expect_equal(pos_fn(tight1, tight2), oc_fn(theta1, theta2), tolerance = 0.01)
})

test_that("pos2S with gaussian(identity) matches fixed-sigma pos2S", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  pos_fixed <- pos2S(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)
  pos_family <- pos2S(prior1, prior2, n1, n2, dec,
                       sigma1 = 2, sigma2 = 2, family = gaussian())

  expect_equal(pos_fixed(prior1, prior2), pos_family(prior1, prior2),
               tolerance = 1e-3)
})
